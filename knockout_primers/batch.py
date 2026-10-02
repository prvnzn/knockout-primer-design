"""Run a whole gene list.

The work splits cleanly in two:

1. **Locate** every gene — needs NCBI (unless a CDS was supplied), but is quick.
   Done one gene at a time on purpose: NCBI asks callers not to hammer Entrez,
   and firing off a dozen parallel lookups is how you get rate-limited.
2. **Design** the primers — no network at all, and the slow part (primer3 spends
   a few seconds per gene on thermodynamic checks). This part runs across
   processes when ``jobs > 1``.

Genes that already carry a complete primer set are scored, never redesigned, and
skip step 2 entirely.
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from typing import Callable, Iterable

from .config import load_config
from .design import GeneResult, assemble, design_pairs, evaluate_primers
from .io import GeneEntry
from .primers import Pair
from .reference import Locus, locate


@dataclass
class BatchOutcome:
    results: list[GeneResult]
    failures: list[tuple[str, str]]          # (gene, reason)


def _reason(exc: Exception) -> str:
    return str(exc) if isinstance(exc, ValueError) else f"{type(exc).__name__}: {exc}"


def _design_one(payload: tuple[Locus, dict]) -> tuple[Pair | None, Pair | None]:
    """Worker entry point: must be importable at module level for pickling."""
    locus, cfg = payload
    return design_pairs(locus, cfg)


def run_batch(entries: Iterable[GeneEntry], cfg: dict | None = None,
              genome: str | None = None, jobs: int = 1,
              on_locate: Callable[[str, bool, str], None] | None = None,
              on_result: Callable[[GeneResult], None] | None = None) -> BatchOutcome:
    """Design or score every gene in ``entries``.

    ``on_locate(gene, ok, message)`` and ``on_result(result)`` are optional
    progress callbacks, so a caller can report as it goes instead of waiting for
    the whole batch.
    """
    cfg = cfg or load_config()
    entries = list(entries)
    results: list[GeneResult] = []
    failures: list[tuple[str, str]] = []

    # --- 1. locate, and score anything that already has primers ---------------
    to_design: list[tuple[GeneEntry, Locus]] = []
    for e in entries:
        try:
            if e.has_complete_set:
                res = evaluate_primers(e.gene, e.existing, e.cds or None, cfg, genome)
                results.append(res)
                if on_locate:
                    on_locate(e.gene, True, "kept existing primers")
                if on_result:
                    on_result(res)
                continue
            locus = locate(e.gene, e.cds or None, cfg, genome=genome)
            to_design.append((e, locus))
            if on_locate:
                on_locate(e.gene, True, f"{locus.strand} strand, ORF {locus.orf_len} bp")
        except Exception as exc:
            failures.append((e.gene, _reason(exc)))
            if on_locate:
                on_locate(e.gene, False, _reason(exc))

    # --- 2. design (the expensive, network-free part) -------------------------
    def finish(entry: GeneEntry, locus: Locus, pairs) -> None:
        inner, outer = pairs
        res = assemble(entry.gene, locus, inner, outer, cfg, "designed", genome)
        results.append(res)
        if on_result:
            on_result(res)

    if jobs > 1 and len(to_design) > 1:
        with ProcessPoolExecutor(max_workers=jobs) as pool:
            payloads = [(locus, cfg) for _, locus in to_design]
            for (entry, locus), pairs in zip(to_design, pool.map(_design_one, payloads)):
                finish(entry, locus, pairs)
    else:
        for entry, locus in to_design:
            try:
                finish(entry, locus, design_pairs(locus, cfg))
            except Exception as exc:
                failures.append((entry.gene, _reason(exc)))

    # Keep the output in the order the genes were given, not the order they finished.
    order = {e.gene: i for i, e in enumerate(entries)}
    results.sort(key=lambda r: order.get(r.gene, 0))
    return BatchOutcome(results, failures)
