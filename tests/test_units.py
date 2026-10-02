"""Offline tests: config handling, input parsing, and the assay geometry.

Nothing here touches the network. The geometry tests build a synthetic locus, so
the band arithmetic and the cassette-control geometry are checked without needing
a genome at all.
"""
from __future__ import annotations

import random

import pytest

from knockout_primers.config import apply_overrides, load_config
from knockout_primers.control import design_control
from knockout_primers.design import _outer_products, _span, design_gene
from knockout_primers.io import read_genes
from knockout_primers.primers import design_inner, design_outer, score
from knockout_primers.reference import Locus, count_genome_hits, find_binding_site


# --------------------------------------------------------------------------- #
# Config
# --------------------------------------------------------------------------- #
def test_defaults_are_self_consistent():
    cfg = load_config()
    P = cfg["primers"]
    assert P["length"]["min"] <= P["length"]["opt"] <= P["length"]["max"]
    assert P["tm"]["min"] <= P["tm"]["opt"] <= P["tm"]["max"]
    assert P["tm"]["accept_min"] <= P["tm"]["min"]
    assert P["tm"]["accept_max"] >= P["tm"]["max"]
    assert P["gc"]["min"] < P["gc"]["max"]
    assert cfg["cassette"]["length_bp"] > 0


def test_overrides_change_values_and_coerce_types():
    cfg = load_config()
    out = apply_overrides(cfg, ["primers.tm.opt=62", "control.enabled=false",
                                "primers.max_poly_x=5"])
    assert out["primers"]["tm"]["opt"] == 62
    assert out["control"]["enabled"] is False
    assert out["primers"]["max_poly_x"] == 5
    assert cfg["primers"]["tm"]["opt"] == 60, "the original config must not be mutated"


def test_unknown_override_is_rejected():
    cfg = load_config()
    with pytest.raises(KeyError):
        apply_overrides(cfg, ["primers.tm.optimum=62"])
    with pytest.raises(KeyError):
        apply_overrides(cfg, ["nonsense.key=1"])
    with pytest.raises(ValueError):
        apply_overrides(cfg, ["primers.tm.opt"])


# --------------------------------------------------------------------------- #
# Input parsing
# --------------------------------------------------------------------------- #
def test_plain_gene_list_with_comments(tmp_path):
    f = tmp_path / "genes.txt"
    f.write_text("# my targets\nhslV\n\nclpP\n  lon  \n", encoding="utf-8")
    assert [e.gene for e in read_genes(f)] == ["hslV", "clpP", "lon"]


def test_csv_with_cds_and_existing_primers(tmp_path):
    f = tmp_path / "genes.csv"
    f.write_text(
        "gene,cds,outer_f,outer_r,inner_f,inner_r\n"
        "hslV,atgc gtac,AAAA,CCCC,GGGG,TTTT\n"
        "clpP,,,,,\n", encoding="utf-8")
    a, b = read_genes(f)
    assert a.cds == "ATGCGTAC", "whitespace stripped and upper-cased"
    assert a.existing == {"O_F": "AAAA", "O_R": "CCCC", "I_F": "GGGG", "I_R": "TTTT"}
    assert a.has_complete_set and not b.has_complete_set


def test_csv_needs_a_gene_column(tmp_path):
    f = tmp_path / "bad.csv"
    f.write_text("gene_name,cds\nhslV,ATGC\n", encoding="utf-8")
    with pytest.raises(ValueError, match="gene"):
        read_genes(f)


# --------------------------------------------------------------------------- #
# Sequence helpers
# --------------------------------------------------------------------------- #
def test_find_binding_site_both_strands():
    template = "AAAAGGGTTTCCCAAAA"
    assert find_binding_site(template, "GGGTTTCCC") == (4, "+", 9)
    assert find_binding_site(template, "GGGAAACCC") == (4, "-", 9)   # reverse complement
    assert find_binding_site(template, "TTTTTTTTTTTT") is None


def test_count_genome_hits_counts_both_strands():
    # One match as written, one as the reverse complement.
    genome = "TTT" + "ACGTACGTAC" + "TTT" + "GTACGTACGT" + "TTT"
    assert count_genome_hits(genome, "ACGTACGTAC") == 2


# --------------------------------------------------------------------------- #
# Assay geometry, on a synthetic locus
# --------------------------------------------------------------------------- #
def _synthetic_locus(orf_len: int, flank: int = 700, seed: int = 42) -> Locus:
    """A made-up locus with balanced base composition, so primer3 can work on it."""
    rng = random.Random(seed)
    seq = "".join(rng.choice("ACGT") for _ in range(orf_len + 2 * flank))
    return Locus(gene="synth", accession="SYNTH", strand="+",
                 lo=flank, hi=flank + orf_len - 1, flank=flank,
                 up=seq[:flank], orf=seq[flank:flank + orf_len],
                 down=seq[flank + orf_len:])


@pytest.mark.parametrize("orf_len", [300, 900, 1500])
def test_outer_shift_is_cassette_minus_orf(orf_len):
    """The shift must be exactly cassette - ORF, whatever the ORF length."""
    cfg = load_config()
    locus = _synthetic_locus(orf_len)
    outer = design_outer(locus, cfg)
    assert outer is not None
    wt, ko, shift, flag = _outer_products(locus, outer, cfg["cassette"]["length_bp"])
    assert flag is None, "outer primers must sit outside the ORF"
    assert shift == cfg["cassette"]["length_bp"] - orf_len
    assert ko - wt == shift


def test_outer_primers_lie_in_the_flanks():
    cfg = load_config()
    locus = _synthetic_locus(900)
    outer = design_outer(locus, cfg)
    (lo, hi), _, _ = _span(locus, outer)
    assert lo < locus.orf_start, "forward primer must be upstream of the ORF"
    assert hi > locus.orf_start + locus.orf_len, "reverse primer must be downstream"


def test_inner_primers_lie_inside_the_orf_and_span_most_of_it():
    cfg = load_config()
    locus = _synthetic_locus(900)
    inner = design_inner(locus, cfg)
    assert inner is not None
    (lo, hi), _, _ = _span(locus, inner)
    assert lo >= locus.orf_start
    assert hi <= locus.orf_start + locus.orf_len
    assert (hi - lo) / locus.orf_len >= cfg["inner"]["min_coverage_frac"]


def test_control_pairs_convergently_with_one_outer_primer():
    """The control primer must face exactly one outer primer, with a sane size."""
    cfg = load_config()
    from knockout_primers.reference import load_cassette
    cassette = load_cassette(cfg)
    locus = _synthetic_locus(900)
    outer = design_outer(locus, cfg)

    ctrl = design_control(locus, outer, cassette, cfg)
    assert ctrl is not None and ctrl.ok, ctrl.note if ctrl else "no control"
    assert ctrl.partner in ("O_F", "O_R")
    # The product runs from that outer primer, across the junction, into the
    # cassette — so it is shorter than flank + cassette and longer than nothing.
    assert 0 < ctrl.product_bp < locus.flank + cfg["cassette"]["length_bp"]

    # And it must genuinely amplify on the knockout allele but not on wild type.
    allele = locus.keio_allele(cassette)
    assert find_binding_site(allele, ctrl.primer.seq) is not None
    assert find_binding_site(locus.full, ctrl.primer.seq) is None


def test_control_can_be_switched_off():
    cfg = apply_overrides(load_config(), ["control.enabled=false"])
    from knockout_primers.reference import load_cassette
    locus = _synthetic_locus(900)
    outer = design_outer(locus, cfg)
    assert design_control(locus, outer, load_cassette(cfg), cfg) is None


def test_scoring_flags_a_bad_primer():
    cfg = load_config()
    bad = score("AAAAAAAAAAAAAAAAAAAA", cfg)          # low GC, low Tm, poly-A run
    assert not bad.ok
    assert any("GC" in w for w in bad.warnings)
    assert any("identical bases" in w for w in bad.warnings)

    good = score("GTAGCTTAATTCTTCGATGGTGTG", cfg)     # a real, validated primer
    assert good.ok, good.warnings
    assert 55 <= good.tm <= 68
