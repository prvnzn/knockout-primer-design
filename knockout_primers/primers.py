"""Primer design and quality scoring, with primer3 underneath.

* :func:`score` — the thermodynamics of one oligo (Tm, GC, hairpin, self-dimer)
  plus a check against every configured criterion.
* :func:`make_pair` — scores two oligos together, adding the Tm match and the
  cross-dimer check.
* :func:`design_inner` — a pair *inside* the ORF, spanning as much of it as the
  criteria allow (the presence/absence probe).
* :func:`design_outer` — a pair in the constant flanks either side of the ORF
  (the size-shift pair).

Every threshold is read from the config. If primer3 cannot satisfy them it
returns nothing rather than quietly relaxing anything, and the caller reports
that as a failure for the gene — see :mod:`knockout_primers.design`.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import primer3

from .reference import Locus


@dataclass
class Primer:
    seq: str
    tm: float
    gc: float
    length: int
    hairpin_tm: float
    homodimer_tm: float
    warnings: list[str] = field(default_factory=list)
    genome_hits: int | None = None      # filled in by the specificity check

    @property
    def ok(self) -> bool:
        return not self.warnings


@dataclass
class Pair:
    fwd: Primer
    rev: Primer
    heterodimer_tm: float
    tm_diff: float
    warnings: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.warnings and self.fwd.ok and self.rev.ok


def _thermo_args(cfg: dict) -> dict:
    t = cfg["primers"]["thermo"]
    return dict(mv_conc=t["salt_monovalent_mM"], dna_conc=t["dna_conc_nM"])


# --------------------------------------------------------------------------- #
# Scoring
# --------------------------------------------------------------------------- #
def score(seq: str, cfg: dict) -> Primer:
    """Compute one oligo's properties and flag every criterion it misses."""
    seq = "".join(seq.split()).upper()
    P, T = cfg["primers"], cfg["primers"]["thermo"]
    ta = _thermo_args(cfg)

    gc = 100.0 * (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0
    tm = primer3.calc_tm(seq, **ta)
    hairpin = primer3.calc_hairpin(seq, **ta).tm
    homodimer = primer3.calc_homodimer(seq, **ta).tm

    w: list[str] = []
    if not (P["length"]["min"] <= len(seq) <= P["length"]["max"]):
        w.append(f"length {len(seq)} outside {P['length']['min']}-{P['length']['max']} nt")
    if not (P["gc"]["min"] <= gc <= P["gc"]["max"]):
        w.append(f"GC {gc:.0f}% outside {P['gc']['min']}-{P['gc']['max']}%")
    if not (P["tm"]["accept_min"] <= tm <= P["tm"]["accept_max"]):
        w.append(f"Tm {tm:.1f} outside {P['tm']['accept_min']}-{P['tm']['accept_max']} C")
    if hairpin > T["max_hairpin_tm"]:
        w.append(f"hairpin Tm {hairpin:.1f} > {T['max_hairpin_tm']} C")
    if homodimer > T["max_self_any_tm"]:
        w.append(f"self-dimer Tm {homodimer:.1f} > {T['max_self_any_tm']} C")
    if _max_poly_run(seq) > P["max_poly_x"]:
        w.append(f"run of {_max_poly_run(seq)} identical bases > {P['max_poly_x']}")

    return Primer(seq, round(tm, 1), round(gc, 1), len(seq),
                  round(hairpin, 1), round(homodimer, 1), w)


def _max_poly_run(seq: str) -> int:
    best = run = 0
    for i, b in enumerate(seq):
        run = run + 1 if i and b == seq[i - 1] else 1
        best = max(best, run)
    return best


def make_pair(fwd_seq: str, rev_seq: str, cfg: dict) -> Pair:
    """Score a forward/reverse pair together: Tm match and cross-dimer."""
    f, r = score(fwd_seq, cfg), score(rev_seq, cfg)
    het = primer3.calc_heterodimer(f.seq, r.seq, **_thermo_args(cfg)).tm
    tm_diff = abs(f.tm - r.tm)

    w: list[str] = []
    if tm_diff > cfg["primers"]["tm"]["max_pair_diff"]:
        w.append(f"pair Tm differs by {tm_diff:.1f} > {cfg['primers']['tm']['max_pair_diff']} C")
    if het > cfg["primers"]["thermo"]["max_dimer_tm"]:
        w.append(f"cross-dimer Tm {het:.1f} > {cfg['primers']['thermo']['max_dimer_tm']} C")
    return Pair(f, r, round(het, 1), round(tm_diff, 1), w)


# --------------------------------------------------------------------------- #
# primer3 settings, built from the config
# --------------------------------------------------------------------------- #
def _global_args(cfg: dict) -> dict:
    P, T = cfg["primers"], cfg["primers"]["thermo"]
    return {
        "PRIMER_OPT_SIZE": P["length"]["opt"],
        "PRIMER_MIN_SIZE": P["length"]["min"],
        "PRIMER_MAX_SIZE": P["length"]["max"],
        "PRIMER_OPT_TM": P["tm"]["opt"],
        "PRIMER_MIN_TM": P["tm"]["min"],
        "PRIMER_MAX_TM": P["tm"]["max"],
        "PRIMER_MIN_GC": P["gc"]["min"],
        "PRIMER_MAX_GC": P["gc"]["max"],
        "PRIMER_PAIR_MAX_DIFF_TM": P["tm"]["max_pair_diff"],
        "PRIMER_MAX_POLY_X": P["max_poly_x"],
        "PRIMER_SALT_MONOVALENT": T["salt_monovalent_mM"],
        "PRIMER_DNA_CONC": T["dna_conc_nM"],
        "PRIMER_MAX_SELF_ANY_TH": T["max_self_any_tm"],
        "PRIMER_MAX_SELF_END_TH": T["max_self_end_tm"],
        "PRIMER_MAX_HAIRPIN_TH": T["max_hairpin_tm"],
        # We only ever take the best pair, and the thermodynamic checks are the
        # expensive part of the search — so don't pay for runners-up.
        "PRIMER_NUM_RETURN": 1,
    }


def _first_pair(res: dict) -> tuple[str, str] | None:
    if res.get("PRIMER_LEFT_NUM_RETURNED", 0) == 0:
        return None
    return res["PRIMER_LEFT_0_SEQUENCE"], res["PRIMER_RIGHT_0_SEQUENCE"]


# --------------------------------------------------------------------------- #
# The two pairs
# --------------------------------------------------------------------------- #
def design_inner(locus: Locus, cfg: dict) -> Pair | None:
    """Design the inner pair: near-terminal primers spanning the ORF.

    Both primers must land inside the ORF, because the whole point is that the
    knockout — where the ORF is gone — gives no product.
    """
    orf, n = locus.orf, locus.orf_len
    window = cfg["inner"]["max_end_window_bp"]
    g = _global_args(cfg)
    g["PRIMER_PRODUCT_SIZE_RANGE"] = [[max(80, n - window), n]]
    g["PRIMER_PRODUCT_OPT_SIZE"] = n
    g["PRIMER_PAIR_WT_PRODUCT_SIZE_LT"] = 1.0      # reward near-full-length products
    res = primer3.design_primers({"SEQUENCE_TEMPLATE": orf}, g)
    fr = _first_pair(res)
    return make_pair(*fr, cfg) if fr else None


def design_outer(locus: Locus, cfg: dict) -> Pair | None:
    """Design the outer pair: primers in the flanks, kept close to the ORF.

    ``SEQUENCE_TARGET`` covers the ORF, which forces primer3 to put both primers
    outside it — in sequence the knockout still has, so the pair amplifies in
    both strains.
    """
    n = locus.orf_len
    m = cfg["outer"]["product_margin"]
    g = _global_args(cfg)
    g["PRIMER_PRODUCT_SIZE_RANGE"] = [[n + m["min"], n + m["max"]]]
    g["PRIMER_PRODUCT_OPT_SIZE"] = n + m["opt"]
    g["PRIMER_PAIR_WT_PRODUCT_SIZE_GT"] = 0.5
    g["PRIMER_PAIR_WT_PRODUCT_SIZE_LT"] = 0.1

    # The largest product we accept is ORF + max margin, so no primer can sit
    # further than that from the ORF however the margin is split between the two
    # sides. Searching beyond it only costs time, so hand primer3 just the window
    # that can actually contribute. (The flank fetched may be much larger.)
    start = max(0, locus.orf_start - m["max"])
    end = min(len(locus.full), locus.orf_start + n + m["max"])
    seq_args = {
        "SEQUENCE_TEMPLATE": locus.full,
        "SEQUENCE_INCLUDED_REGION": [start, end - start],
        "SEQUENCE_TARGET": [locus.orf_start, n],      # forces primers outside the ORF
    }
    res = primer3.design_primers(seq_args, g)
    fr = _first_pair(res)
    return make_pair(*fr, cfg) if fr else None
