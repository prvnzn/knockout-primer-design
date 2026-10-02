"""Reading the gene list and writing the results.

**Input.** The simplest thing that works: a text file with one gene name per
line, or a CSV with a ``gene`` column. Blank lines and lines starting with ``#``
are ignored, so a list can carry comments.

A CSV may also carry optional columns, in any order, with any capitalisation:

    gene        the gene name (required)
    cds         the coding sequence, if you want the gene located by alignment
                rather than by name
    inner_f, inner_r, outer_f, outer_r
                primers you already have. When all four are given, that gene is
                scored and sized but never redesigned.

**Output.** One CSV with four rows per gene — outer forward/reverse, then inner
forward/reverse — carrying each primer's sequence and properties, both outer band
sizes, the inner band size, the cassette-control row, and any warnings. It opens
in Excel and pastes into an ordering sheet as-is.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

from .design import GeneResult

PRIMER_COLUMNS = {
    "inner_f": "I_F", "inner_r": "I_R", "outer_f": "O_F", "outer_r": "O_R",
    "i_f": "I_F", "i_r": "I_R", "o_f": "O_F", "o_r": "O_R",
}


@dataclass
class GeneEntry:
    """One gene as read from the input file."""
    gene: str
    cds: str = ""
    existing: dict[str, str] = field(default_factory=dict)

    @property
    def has_complete_set(self) -> bool:
        return all(self.existing.get(k) for k in ("I_F", "I_R", "O_F", "O_R"))


def _looks_like_csv(first_line: str) -> bool:
    return ("," in first_line or ";" in first_line) and "gene" in first_line.lower()


def read_genes(path: str | Path) -> list[GeneEntry]:
    """Read a gene list: plain one-per-line, or a CSV with a ``gene`` column."""
    text = Path(path).read_text(encoding="utf-8-sig")
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.lstrip().startswith("#")]
    if not lines:
        raise ValueError(f"{path} contains no gene names")

    if not _looks_like_csv(lines[0]):
        # Plain list. Tolerate a stray trailing comma or a one-column CSV.
        return [GeneEntry(ln.split(",")[0].strip()) for ln in lines]

    delim = ";" if lines[0].count(";") > lines[0].count(",") else ","
    reader = csv.DictReader(lines, delimiter=delim)
    norm = {(h or "").strip().lower(): (h or "") for h in (reader.fieldnames or [])}
    if "gene" not in norm:
        raise ValueError(f"{path}: a CSV input needs a 'gene' column "
                         f"(found: {', '.join(reader.fieldnames or []) or 'nothing'})")

    entries: list[GeneEntry] = []
    for row in reader:
        name = (row.get(norm["gene"]) or "").strip()
        if not name:
            continue
        e = GeneEntry(name)
        if "cds" in norm:
            e.cds = "".join((row.get(norm["cds"]) or "").split()).upper()
        for col, key in PRIMER_COLUMNS.items():
            if col in norm:
                seq = "".join((row.get(norm[col]) or "").split()).upper()
                if seq:
                    e.existing[key] = seq
        entries.append(e)
    return entries


# --------------------------------------------------------------------------- #
# Output
# --------------------------------------------------------------------------- #
HEADER = [
    "Gene", "Strand", "ORF (bp)", "Source", "Status",
    "Primer", "Direction", "Sequence (5'-3')", "Length", "GC%", "Tm (C)",
    "Genome hits", "Band: wild type (bp)", "Band: knockout (bp)", "Notes",
]


def _primer_rows(res: GeneResult) -> list[list]:
    """Four rows per gene (outer F/R, inner F/R), plus a control row if present."""
    base = [res.gene, res.strand, res.orf_len, res.source, "ok" if res.ok else "check"]
    notes = "; ".join(res.warnings) if res.warnings else "ok"
    rows: list[list] = []

    def row(kind, direction, primer, wt, ko, note=""):
        return base + [
            kind, direction,
            primer.seq if primer else "",
            primer.length if primer else "",
            primer.gc if primer else "",
            primer.tm if primer else "",
            primer.genome_hits if primer and primer.genome_hits is not None else "",
            wt, ko, note,
        ]

    o_wt = res.outer_product_wt or ""
    o_ko = res.outer_product_ko or ""
    of = res.outer.fwd if res.outer else None
    orv = res.outer.rev if res.outer else None
    shift = f"shift {res.outer_shift:+d} bp" if res.outer_shift is not None else ""
    rows.append(row("Outer", "Forward", of, o_wt, o_ko, shift))
    rows.append(row("Outer", "Reverse", orv, o_wt, o_ko, ""))

    i_wt = res.inner_product_wt or ""
    inf = res.inner.fwd if res.inner else None
    inr = res.inner.rev if res.inner else None
    cov = (f"covers {res.inner_coverage_frac:.0%} of ORF"
           if res.inner_coverage_frac is not None else "")
    rows.append(row("Inner", "Forward", inf, i_wt, "none", cov))
    rows.append(row("Inner", "Reverse", inr, i_wt, "none", ""))

    c = res.control
    if c is not None and c.ok:
        rows.append(row("Control", "Forward", c.primer, "none", c.product_bp,
                        f"pair with {c.partner} ({c.partner_seq})"))

    rows[-1][-1] = (rows[-1][-1] + "; " if rows[-1][-1] else "") + notes
    return rows


def write_results_csv(results: list[GeneResult], path: str | Path) -> None:
    """Write every gene's set to one CSV."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(HEADER)
        for res in results:
            w.writerows(_primer_rows(res))


def write_failures_csv(failures: list[tuple[str, str]], path: str | Path) -> None:
    """Write the genes that could not be designed, with the reason for each."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Gene", "Reason"])
        w.writerows(failures)
