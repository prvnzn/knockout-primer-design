"""Design a full verification primer set for one gene.

The assay, in one block:

    inner pair   both primers inside the ORF
                 wild type -> band (M_inner)        knockout -> no product
    outer pair   both primers in the flanks
                 wild type -> M_wt                  knockout -> M_ko
    control      cassette-internal primer + the outer primer in the flank it faces
                 wild type -> no product            knockout -> band

The outer shift is fixed and worth internalising:

    M_ko - M_wt = cassette_length - ORF_length

The flanks are identical in both strains, so they cancel and the placement of
the outer primers does not affect the shift at all. A short ORF therefore shifts
*up* strongly, and an ORF close to the cassette length barely moves — those genes
are flagged, and for them the inner pair is what you trust.

Two entry points, both returning a :class:`GeneResult`:

* :func:`design_gene` — design fresh pairs with primer3.
* :func:`evaluate_primers` — score and size primers you already have, so an
  existing, trusted set is never silently replaced.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Optional

from .config import load_config
from .control import ControlResult, design_control
from .primers import Pair, design_inner, design_outer, make_pair
from .reference import (Locus, count_genome_hits, find_binding_site, load_cassette,
                        locate)


@dataclass
class GeneResult:
    """Everything the tool knows about one gene's verification set."""
    gene: str
    ok: bool
    source: str                            # 'designed' or 'existing'
    strand: str
    orf_len: int
    product: str                           # NCBI description of the gene

    inner: Optional[Pair]
    inner_product_wt: Optional[int]        # band in wild type; none in the knockout
    inner_coverage_frac: Optional[float]   # fraction of the ORF the inner band spans

    outer: Optional[Pair]
    outer_product_wt: Optional[int]        # M_wt
    outer_product_ko: Optional[int]        # M_ko
    outer_shift: Optional[int]             # M_ko - M_wt

    control: Optional[ControlResult]
    warnings: list[str]

    def as_dict(self) -> dict:
        d = asdict(self)
        d["control"] = asdict(self.control) if self.control else None
        return d


# --------------------------------------------------------------------------- #
# Amplicon geometry
# --------------------------------------------------------------------------- #
def _span(locus: Locus, pair: Pair):
    """Where a pair's product sits on the locus: ((lo, hi), starts, ends) or None.

    Orientation-agnostic: a working pair is convergent, so the product runs from
    the leftmost primer start to the rightmost primer end whichever primer is
    labelled forward.
    """
    hf = find_binding_site(locus.full, pair.fwd.seq)
    hr = find_binding_site(locus.full, pair.rev.seq)
    if not hf or not hr:
        return None
    starts = [hf[0], hr[0]]
    ends = [hf[0] + hf[2], hr[0] + hr[2]]
    return (min(starts), max(ends)), starts, ends


def _product_wt(locus: Locus, pair: Pair) -> Optional[int]:
    span = _span(locus, pair)
    return None if span is None else span[0][1] - span[0][0]


def _outer_products(locus: Locus, pair: Pair, cassette_len: int):
    """``(M_wt, M_ko, shift, flag)`` for the outer pair."""
    span = _span(locus, pair)
    if span is None:
        return None, None, None, "unmapped"
    (lo, hi), starts, ends = span
    wt = hi - lo
    orf_lo, orf_hi = locus.orf_start, locus.orf_start + locus.orf_len
    # A primer that overlaps the ORF is deleted in the knockout -> no product there.
    overlaps = any(not (e <= orf_lo or s >= orf_hi) for s, e in zip(starts, ends))
    shift = cassette_len - locus.orf_len
    return wt, wt + shift, shift, ("orf_overlap" if overlaps else None)


def _inner_inside_orf(locus: Locus, pair: Pair) -> bool:
    """True when both inner primers really do sit within the ORF."""
    span = _span(locus, pair)
    if span is None:
        return False
    (lo, hi), _, _ = span
    return lo >= locus.orf_start and hi <= locus.orf_start + locus.orf_len


# --------------------------------------------------------------------------- #
# Assembly
# --------------------------------------------------------------------------- #
def _assemble(gene, locus, inner, outer, cfg, source, genome=None) -> GeneResult:
    """Turn a pair of pairs into sizes, a pass/fail verdict, and warnings."""
    cassette = load_cassette(cfg)
    cassette_len = cfg["cassette"]["length_bp"]
    warnings: list[str] = []

    inner_wt = coverage = None
    if inner is None:
        warnings.append("inner: no pair meets the criteria")
    else:
        inner_wt = _product_wt(locus, inner)
        if inner_wt is None:
            warnings.append("inner: the primers do not both bind this locus")
        else:
            coverage = round(inner_wt / locus.orf_len, 3)
            if coverage < cfg["inner"]["min_coverage_frac"]:
                warnings.append(f"inner band covers only {coverage:.0%} of the ORF")
            if not _inner_inside_orf(locus, inner):
                warnings.append("inner: a primer lies outside the ORF — it would still "
                                "amplify in the knockout")
        warnings += [f"inner-F: {w}" for w in inner.fwd.warnings]
        warnings += [f"inner-R: {w}" for w in inner.rev.warnings]
        warnings += [f"inner pair: {w}" for w in inner.warnings]

    outer_wt = outer_ko = shift = None
    if outer is None:
        warnings.append("outer: no pair meets the criteria")
    else:
        outer_wt, outer_ko, shift, flag = _outer_products(locus, outer, cassette_len)
        if flag == "unmapped":
            warnings.append("outer: the primers do not both bind this locus")
        elif flag == "orf_overlap":
            warnings.append("outer: a primer overlaps the ORF, so it is deleted in the "
                            "knockout — move it further into the flank")
        if shift is not None and abs(shift) < cfg["outer"]["min_shift_bp"]:
            warnings.append(
                f"outer bands differ by only {shift:+d} bp (ORF {locus.orf_len} bp vs "
                f"cassette {cassette_len} bp) — hard to tell apart on a gel, rely on "
                "the inner pair")
        warnings += [f"outer-F: {w}" for w in outer.fwd.warnings]
        warnings += [f"outer-R: {w}" for w in outer.rev.warnings]
        warnings += [f"outer pair: {w}" for w in outer.warnings]

    control = design_control(locus, outer, cassette, cfg)
    if control is not None and control.note:
        warnings.append(f"control: {control.note}")

    # Specificity: an exact match anywhere else in the genome risks extra bands.
    if genome is not None:
        for tag, pair in (("inner", inner), ("outer", outer)):
            for role, primer in (("F", pair.fwd), ("R", pair.rev)) if pair else ():
                primer.genome_hits = count_genome_hits(genome, primer.seq)
                if primer.genome_hits > 1:
                    warnings.append(f"{tag}-{role}: binds {primer.genome_hits} places in "
                                    "the genome — not specific")

    ok = bool(
        inner and outer and inner.ok and outer.ok
        and inner_wt and outer_wt
        and not any(w.startswith(("inner:", "outer:")) for w in warnings)
        and not any("differ by only" in w or "not specific" in w for w in warnings)
    )
    return GeneResult(
        gene=gene, ok=ok, source=source, strand=locus.strand, orf_len=locus.orf_len,
        product=locus.product,
        inner=inner, inner_product_wt=inner_wt, inner_coverage_frac=coverage,
        outer=outer, outer_product_wt=outer_wt, outer_product_ko=outer_ko,
        outer_shift=shift, control=control, warnings=warnings,
    )


# --------------------------------------------------------------------------- #
# Forward / reverse convention
# --------------------------------------------------------------------------- #
def _swap_roles(pair: Pair | None) -> Pair | None:
    """Swap a pair's forward/reverse labels. The sequences are untouched.

    primer3 designs on the gene-sense template, so for a gene on the minus strand
    its 'left' primer is the one that reads leftwards along the chromosome. We
    report forward/reverse relative to the reference (+) strand — the usual
    convention — so those labels have to be swapped back.
    """
    if pair is None:
        return None
    return Pair(fwd=pair.rev, rev=pair.fwd,
                heterodimer_tm=pair.heterodimer_tm, tm_diff=pair.tm_diff,
                warnings=pair.warnings)


# --------------------------------------------------------------------------- #
# Entry points
# --------------------------------------------------------------------------- #
def design_pairs(locus: Locus, cfg: dict) -> tuple[Pair | None, Pair | None]:
    """Run primer3 for one located gene and return ``(inner, outer)``.

    Pure computation: no network, no files. Kept separate from :func:`design_gene`
    so a batch can run this part in parallel across processes — it is by far the
    most expensive step.
    """
    inner, outer = design_inner(locus, cfg), design_outer(locus, cfg)
    if locus.strand == "-":
        inner, outer = _swap_roles(inner), _swap_roles(outer)
    return inner, outer


def assemble(gene: str, locus: Locus, inner: Pair | None, outer: Pair | None,
             cfg: dict, source: str = "designed", genome: str | None = None) -> GeneResult:
    """Turn located-and-designed pairs into the finished :class:`GeneResult`."""
    return _assemble(gene, locus, inner, outer, cfg, source, genome)


def design_gene(gene: str, cds: str | None = None, cfg: dict | None = None,
                genome: str | None = None) -> GeneResult:
    """Design a fresh verification set for one gene.

    >>> from knockout_primers.design import design_gene
    >>> r = design_gene("hslV")
    >>> r.outer_product_wt, r.outer_product_ko
    """
    cfg = cfg or load_config()
    locus = locate(gene, cds, cfg, genome=genome)
    inner, outer = design_pairs(locus, cfg)
    return _assemble(gene, locus, inner, outer, cfg, "designed", genome)


def evaluate_primers(gene: str, primers: dict, cds: str | None = None,
                     cfg: dict | None = None, genome: str | None = None) -> GeneResult:
    """Score and size primers you already have — nothing is redesigned.

    ``primers`` holds any of ``I_F``, ``I_R``, ``O_F``, ``O_R`` (5'->3').
    """
    cfg = cfg or load_config()
    locus = locate(gene, cds, cfg, genome=genome)
    have = lambda a, b: primers.get(a) and primers.get(b)
    inner = make_pair(primers["I_F"], primers["I_R"], cfg) if have("I_F", "I_R") else None
    outer = make_pair(primers["O_F"], primers["O_R"], cfg) if have("O_F", "O_R") else None
    return _assemble(gene, locus, inner, outer, cfg, "existing", genome)


def design_or_keep(gene: str, cds: str | None = None, existing: dict | None = None,
                   cfg: dict | None = None, genome: str | None = None) -> GeneResult:
    """Keep a complete existing set; design one for any gene that lacks it."""
    cfg = cfg or load_config()
    if existing and all(existing.get(k) for k in ("I_F", "I_R", "O_F", "O_R")):
        return evaluate_primers(gene, existing, cds, cfg, genome)
    return design_gene(gene, cds, cfg, genome)
