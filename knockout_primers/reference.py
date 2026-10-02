"""Reference-genome layer: find a gene and cut out its ORF plus flanks.

A gene can be located two ways, both giving the same :class:`Locus`:

* **by name** (:func:`locate_by_name`) — ask NCBI Entrez for the annotated
  coordinates. This is all most users need: a list of gene names is enough.
* **by coding sequence** (:func:`locate_by_cds`) — align a CDS you supply to the
  chromosome. Name-independent and exact; useful for a gene whose name is
  ambiguous, or to pin down a specific allele.

The chromosome is downloaded from NCBI once and cached in ``cache/``, so a whole
batch costs a single download and every later run is offline.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from Bio import Entrez, SeqIO
from Bio.Seq import Seq

from .config import load_config, resolve_path

# Where the downloaded chromosome is kept. Relative to wherever the tool is run,
# so a project keeps its own cache; override with --cache-dir.
CACHE_DIR = Path("cache")


@dataclass
class Locus:
    """A located gene, written in the gene's own 5'->3' orientation.

    ``up``/``orf``/``down`` are the upstream flank, the ORF, and the downstream
    flank. In the Keio strain ``orf`` is what gets replaced by the cassette; the
    flanks are identical in both strains, which is why the outer pair works.
    """
    gene: str
    accession: str
    strand: str          # '+' or '-' relative to the reference chromosome
    lo: int              # 0-based inclusive start on the (+) strand
    hi: int              # 0-based inclusive end on the (+) strand
    flank: int           # bp of flanking context captured each side
    up: str
    orf: str
    down: str
    product: str = ""    # NCBI product description, when a name lookup found one

    @property
    def full(self) -> str:
        """The whole window, 5'->3' in gene orientation: up + orf + down."""
        return self.up + self.orf + self.down

    @property
    def orf_len(self) -> int:
        return len(self.orf)

    @property
    def orf_start(self) -> int:
        """Offset of the ORF within :attr:`full` (= the upstream flank actually captured)."""
        return len(self.up)

    def keio_allele(self, cassette: str) -> str:
        """The modelled knockout locus: flanks unchanged, ORF replaced by the cassette.

        The Keio allele is modelled as a clean, complete-ORF replacement (Baba
        et al. 2006). Real junctions carry a few scar bases, so sizes computed
        from this are accurate to within a handful of bp.
        """
        return self.up + cassette + self.down


# --------------------------------------------------------------------------- #
# Genome download + cache
# --------------------------------------------------------------------------- #
def _set_email(cfg: dict) -> None:
    email = (cfg.get("reference", {}) or {}).get("entrez_email") or ""
    if not email.strip():
        raise ValueError(
            "NCBI requires an e-mail address for Entrez queries. Set "
            "reference.entrez_email in the config, or pass --email you@example.com"
        )
    Entrez.email = email.strip()


def fetch_genome(cfg: dict | None = None, cache_dir: Path | None = None) -> str:
    """Return the reference chromosome as an uppercase string, caching the FASTA."""
    cfg = cfg or load_config()
    acc = cfg["reference"]["accession"]
    cache = (cache_dir or CACHE_DIR) / f"{acc}.fasta"
    if cache.exists():
        return str(SeqIO.read(cache, "fasta").seq).upper()
    _set_email(cfg)
    cache.parent.mkdir(parents=True, exist_ok=True)
    with Entrez.efetch(db="nuccore", id=acc, rettype="fasta", retmode="text") as fh:
        rec = SeqIO.read(fh, "fasta")
    SeqIO.write(rec, cache, "fasta")
    return str(rec.seq).upper()


def load_cassette(cfg: dict) -> str:
    """Return the cassette sequence, written in the gene's 5'->3' orientation (P1->P4)."""
    path = resolve_path(cfg["cassette"]["sequence_fasta"], cfg)
    seq = str(SeqIO.read(path, "fasta").seq).upper()
    declared = cfg["cassette"]["length_bp"]
    if len(seq) != declared:
        raise ValueError(
            f"cassette FASTA is {len(seq)} bp but cassette.length_bp says {declared}; "
            "fix one of them (the band sizes depend on it)"
        )
    return seq


# --------------------------------------------------------------------------- #
# Locating a gene
# --------------------------------------------------------------------------- #
def _build_locus(gene, accession, strand, lo, hi, genome, flank, product="") -> Locus:
    """Cut the window out of the genome and orient it 5'->3' along the gene."""
    win = genome[max(0, lo - flank): hi + 1 + flank]
    orf_len = hi - lo + 1
    left = lo - max(0, lo - flank)          # flank actually captured (genome-edge safe)
    seq = win if strand == "+" else str(Seq(win).reverse_complement())
    if strand == "-":
        left = len(win) - orf_len - left    # the other flank leads after reverse-complementing
    return Locus(gene, accession, strand, lo, hi, flank,
                 up=seq[:left], orf=seq[left:left + orf_len], down=seq[left + orf_len:],
                 product=product)


def locate_by_cds(gene: str, cds: str, cfg: dict | None = None,
                  flank: int | None = None, genome: str | None = None) -> Locus:
    """Locate a gene by exact alignment of a CDS to the chromosome (either strand)."""
    cfg = cfg or load_config()
    flank = cfg["reference"]["flank_bp"] if flank is None else flank
    genome = genome if genome is not None else fetch_genome(cfg)
    acc = cfg["reference"]["accession"]
    cds = "".join(cds.split()).upper()

    i = genome.find(cds)
    if i >= 0:
        return _build_locus(gene, acc, "+", i, i + len(cds) - 1, genome, flank)
    j = genome.find(str(Seq(cds).reverse_complement()))
    if j >= 0:
        return _build_locus(gene, acc, "-", j, j + len(cds) - 1, genome, flank)
    raise ValueError(
        f"the CDS given for {gene!r} is not an exact match anywhere in {acc} "
        "(checked both strands) — check the sequence, or drop it and let the "
        "gene name be looked up instead"
    )


def locate_by_name(gene: str, cfg: dict | None = None, flank: int | None = None,
                   genome: str | None = None) -> Locus:
    """Locate a gene by its NCBI-annotated coordinates."""
    cfg = cfg or load_config()
    flank = cfg["reference"]["flank_bp"] if flank is None else flank
    genome = genome if genome is not None else fetch_genome(cfg)
    org = cfg["reference"]["organism"]
    _set_email(cfg)

    lookup = (cfg["reference"].get("gene_aliases") or {}).get(gene, gene)
    with Entrez.esearch(db="gene", term=f'{lookup}[gene] AND "{org}"[orgn]') as fh:
        res = Entrez.read(fh)
    if not res["IdList"]:
        raise ValueError(
            f"no NCBI gene record for {gene!r} in {org} — check the spelling, add it "
            "to reference.gene_aliases, or supply the gene's sequence instead"
        )
    with Entrez.esummary(db="gene", id=res["IdList"][0]) as fh:
        summ = Entrez.read(fh)
    doc = summ["DocumentSummarySet"]["DocumentSummary"][0]
    info = doc["GenomicInfo"][0]
    start, stop = int(info["ChrStart"]), int(info["ChrStop"])
    strand = "+" if start < stop else "-"
    return _build_locus(gene, cfg["reference"]["accession"], strand,
                        min(start, stop), max(start, stop), genome, flank,
                        product=str(doc.get("Description", "")))


def locate(gene: str, cds: str | None = None, cfg: dict | None = None,
           flank: int | None = None, genome: str | None = None) -> Locus:
    """Locate a gene, preferring a supplied CDS and falling back to its name."""
    cfg = cfg or load_config()
    if cds and cds.strip():
        return locate_by_cds(gene, cds, cfg, flank, genome)
    return locate_by_name(gene, cfg, flank, genome)


# --------------------------------------------------------------------------- #
# Mapping a primer onto a sequence
# --------------------------------------------------------------------------- #
def find_binding_site(template: str, primer: str) -> tuple[int, str, int] | None:
    """Find a primer on ``template`` (given as a 5'->3' top strand).

    Returns ``(start, strand, length)``, where strand is ``'+'`` when the primer
    matches the top strand as written and ``'-'`` when it matches the bottom
    strand. Only exact, full-length matches count. ``None`` if it does not bind.
    """
    p = "".join(primer.split()).upper()
    if not p:
        return None
    i = template.find(p)
    if i >= 0:
        return (i, "+", len(p))
    j = template.find(str(Seq(p).reverse_complement()))
    if j >= 0:
        return (j, "-", len(p))
    return None


def count_genome_hits(genome: str, primer: str) -> int:
    """How many places in the genome a primer matches exactly (both strands).

    More than one means the primer is not specific and may give extra bands.
    """
    p = "".join(primer.split()).upper()
    if not p:
        return 0
    def n(hay: str, needle: str) -> int:
        c = start = 0
        while (k := hay.find(needle, start)) >= 0:
            c += 1
            start = k + 1
        return c
    rc = str(Seq(p).reverse_complement())
    return n(genome, p) + (n(genome, rc) if rc != p else 0)
