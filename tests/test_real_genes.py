"""End-to-end tests against real E. coli genes.

These need NCBI access the first time (the chromosome is cached afterwards), so
they are marked ``network``:

    pytest -m "not network"     # skip them
    pytest                      # run everything

The expected values are not guesses: they were produced by this assay's original
implementation and checked at the bench, so a change in any of them means the
geometry has regressed.
"""
from __future__ import annotations

import os

import pytest

from knockout_primers.config import apply_overrides, load_config
from knockout_primers.design import design_gene
from knockout_primers.reference import fetch_genome

pytestmark = pytest.mark.network

EMAIL = os.environ.get("NCBI_EMAIL", "keio.primers.tests@example.com")


@pytest.fixture(scope="module")
def cfg():
    return apply_overrides(load_config(), [f"reference.entrez_email={EMAIL}"])


@pytest.fixture(scope="module")
def genome(cfg):
    return fetch_genome(cfg)


# gene,  strand, ORF length, expected outer shift, expected control band
REAL_GENES = [
    ("hslV", "-", 531, 1323 - 531, 637),
    ("clpP", "+", 624, 1323 - 624, None),
    ("lon", "+", 2355, 1323 - 2355, None),
    ("pepQ", "+", 1332, 1323 - 1332, 568),
]


@pytest.mark.parametrize("gene,strand,orf_len,shift,control_bp", REAL_GENES)
def test_gene_geometry(gene, strand, orf_len, shift, control_bp, cfg, genome):
    res = design_gene(gene, cfg=cfg, genome=genome)

    assert res.strand == strand
    assert res.orf_len == orf_len
    assert res.outer_shift == shift
    assert res.outer_product_ko - res.outer_product_wt == shift
    assert res.inner_product_wt and res.inner_product_wt <= orf_len
    assert res.inner_coverage_frac >= cfg["inner"]["min_coverage_frac"]
    if control_bp is not None:
        assert res.control is not None and res.control.product_bp == control_bp


def test_hslV_primers_are_stable(cfg, genome):
    """The exact oligos for a minus-strand gene, including forward/reverse labels.

    This pins the convention: forward and reverse are reported relative to the
    reference (+) strand, not to the gene's own direction.
    """
    res = design_gene("hslV", cfg=cfg, genome=genome)
    assert res.outer.fwd.seq == "TCATTTCAGACATGGGAGATCC"
    assert res.outer.rev.seq == "GGTCATACAAACTGCATTCGG"
    assert res.control.partner == "O_R", "the control faces the primer in the 5' flank"


def test_small_shift_is_flagged(cfg, genome):
    """pepQ is 1332 bp against a 1323 bp cassette: the bands are 9 bp apart."""
    res = design_gene("pepQ", cfg=cfg, genome=genome)
    assert abs(res.outer_shift) < cfg["outer"]["min_shift_bp"]
    assert not res.ok
    assert any("differ by only" in w for w in res.warnings)


def test_designed_primers_are_specific(cfg, genome):
    res = design_gene("hslV", cfg=cfg, genome=genome)
    for pair in (res.inner, res.outer):
        for primer in (pair.fwd, pair.rev):
            assert primer.genome_hits == 1, f"{primer.seq} is not unique"


def test_unknown_gene_name_fails_clearly(cfg, genome):
    with pytest.raises(ValueError, match="no NCBI gene record"):
        design_gene("notagene999", cfg=cfg, genome=genome)


def test_supplied_cds_locates_the_same_gene(cfg, genome):
    """Locating by sequence must agree with locating by name."""
    from knockout_primers.reference import locate_by_name
    by_name = design_gene("hslV", cfg=cfg, genome=genome)
    cds = locate_by_name("hslV", cfg, genome=genome).orf
    by_cds = design_gene("hslV", cds=cds, cfg=cfg, genome=genome)
    assert by_cds.strand == by_name.strand
    assert by_cds.orf_len == by_name.orf_len
    assert by_cds.outer.fwd.seq == by_name.outer.fwd.seq
