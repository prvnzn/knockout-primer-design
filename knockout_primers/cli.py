"""Command-line interface.

    knockout-primers genes.txt --email you@example.com -o primers.csv

Or, without installing the package:

    python -m knockout_primers genes.txt --email you@example.com -o primers.csv

Every design criterion can be changed from here (``--tm-opt``, ``--gc``, …), or
in bulk by copying ``config/criteria.yaml`` and passing ``--config``. Anything
without a dedicated flag is reachable with ``--set section.key=value``.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

from . import __version__
from .batch import run_batch
from .config import apply_overrides, load_config, write_default_config
from .io import read_genes, write_failures_csv, write_results_csv
from .reference import fetch_genome


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="knockout-primers",
        description="Design nested colony-PCR primers that verify Keio-collection "
                    "gene deletions in E. coli.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""examples:
  knockout-primers genes.txt --email you@example.com
  knockout-primers genes.csv -o out/primers.csv --tm-opt 62 --gc 45 65
  knockout-primers genes.txt --no-control --set outer.min_shift_bp=100
""")
    p.add_argument("genes", type=Path, nargs="?",
                   help="gene list: one name per line, or a CSV with a 'gene' column")
    p.add_argument("-o", "--out", type=Path, default=Path("primers.csv"),
                   help="output CSV (default: primers.csv)")
    p.add_argument("--report", type=Path,
                   help="also write a JSON report with every number and warning")
    p.add_argument("--config", type=Path, help="criteria YAML to use instead of the default")
    p.add_argument("--write-config", type=Path, metavar="FILE",
                   help="write a copy of the default criteria file for editing, then exit")
    p.add_argument("--cache-dir", type=Path, default=Path("cache"), metavar="DIR",
                   help="where to keep the downloaded genome (default: ./cache)")
    p.add_argument("-j", "--jobs", type=int, metavar="N",
                   help="design this many genes at once (default: one per core). "
                        "Gene lookups stay sequential either way, to be kind to NCBI")
    p.add_argument("--email", help="your e-mail address; NCBI requires one for lookups")
    p.add_argument("--version", action="version", version=f"knockout-primers {__version__}")

    g = p.add_argument_group("design criteria (override the config)")
    g.add_argument("--tm-opt", type=float, metavar="C", help="target melting temperature")
    g.add_argument("--tm", type=float, nargs=2, metavar=("MIN", "MAX"),
                   help="melting-temperature window to design within")
    g.add_argument("--tm-accept", type=float, nargs=2, metavar=("MIN", "MAX"),
                   help="melting temperatures to accept before flagging a primer")
    g.add_argument("--tm-pair-diff", type=float, metavar="C",
                   help="largest Tm difference allowed within a pair")
    g.add_argument("--gc", type=float, nargs=2, metavar=("MIN", "MAX"), help="GC%% window")
    g.add_argument("--length", type=int, nargs=2, metavar=("MIN", "MAX"),
                   help="primer length window (nt)")
    g.add_argument("--cassette-bp", type=int, metavar="BP",
                   help="length of the inserted cassette (default 1323, pKD13 FRT-kan-FRT)")
    g.add_argument("--flank-bp", type=int, metavar="BP",
                   help="genomic context fetched each side of the ORF")
    g.add_argument("--min-shift", type=int, metavar="BP",
                   help="flag a gene when the two outer bands differ by less than this")
    g.add_argument("--no-control", action="store_true",
                   help="skip the cassette (kanamycin) control")
    g.add_argument("--set", dest="overrides", action="append", metavar="KEY=VALUE",
                   help="set any config value, e.g. --set primers.max_poly_x=5 "
                        "(repeatable)")
    return p


def _config_from_args(args) -> dict:
    """Load the config, then layer the command-line criteria on top of it."""
    cfg = load_config(args.config)
    overrides = list(args.overrides or [])

    def add(key, value):
        if value is not None:
            overrides.append(f"{key}={value}")

    add("primers.tm.opt", args.tm_opt)
    if args.tm:
        add("primers.tm.min", args.tm[0]); add("primers.tm.max", args.tm[1])
    if args.tm_accept:
        add("primers.tm.accept_min", args.tm_accept[0])
        add("primers.tm.accept_max", args.tm_accept[1])
    add("primers.tm.max_pair_diff", args.tm_pair_diff)
    if args.gc:
        add("primers.gc.min", args.gc[0]); add("primers.gc.max", args.gc[1])
    if args.length:
        add("primers.length.min", args.length[0]); add("primers.length.max", args.length[1])
    add("cassette.length_bp", args.cassette_bp)
    add("reference.flank_bp", args.flank_bp)
    add("outer.min_shift_bp", args.min_shift)
    if args.no_control:
        add("control.enabled", "false")
    if args.email:
        add("reference.entrez_email", args.email)

    cfg = apply_overrides(cfg, overrides)
    if args.length and cfg["primers"]["length"]["opt"] > cfg["primers"]["length"]["max"]:
        cfg["primers"]["length"]["opt"] = cfg["primers"]["length"]["max"]
    return cfg


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.write_config:
        try:
            out = write_default_config(args.write_config)
        except (OSError, FileExistsError) as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 2
        print(f"wrote {out}")
        print(f"edit it, then run:  knockout-primers genes.txt --config {out}")
        return 0

    if args.genes is None:
        parser.print_usage(sys.stderr)
        print("error: a gene list is required (or use --write-config)", file=sys.stderr)
        return 2

    try:
        cfg = _config_from_args(args)
    except (KeyError, ValueError, FileNotFoundError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    try:
        entries = read_genes(args.genes)
    except (OSError, ValueError) as exc:
        print(f"error reading {args.genes}: {exc}", file=sys.stderr)
        return 2

    print(f"{len(entries)} gene(s) from {args.genes}")
    print(f"criteria: {cfg['_config_path']}")
    print(f"Tm {cfg['primers']['tm']['min']}-{cfg['primers']['tm']['max']} C "
          f"(target {cfg['primers']['tm']['opt']}), "
          f"GC {cfg['primers']['gc']['min']}-{cfg['primers']['gc']['max']}%, "
          f"length {cfg['primers']['length']['min']}-{cfg['primers']['length']['max']} nt, "
          f"cassette {cfg['cassette']['length_bp']} bp")

    try:
        genome = fetch_genome(cfg, cache_dir=args.cache_dir)
    except Exception as exc:
        print(f"error: could not get the reference genome: {exc}", file=sys.stderr)
        return 2
    print(f"reference {cfg['reference']['accession']}: {len(genome):,} bp\n")

    jobs = args.jobs if args.jobs and args.jobs > 0 else (os.cpu_count() or 1)
    jobs = max(1, min(jobs, len(entries)))
    if jobs > 1:
        print(f"designing on {jobs} cores\n")

    def on_locate(gene: str, ok: bool, message: str) -> None:
        if not ok:
            print(f"  FAILED {gene}: {message}")

    def on_result(res) -> None:
        mark = "ok   " if res.ok else "check"
        shift = f" ({res.outer_shift:+d} bp)" if res.outer_shift is not None else ""
        ctrl = f", control {res.control.product_bp} bp" if res.control and res.control.ok else ""
        kept = " [kept]" if res.source == "existing" else ""
        print(f"  {mark} {res.gene:<8} {res.strand}  ORF {res.orf_len:>5} bp  "
              f"inner {res.inner_product_wt or '-'}  "
              f"outer {res.outer_product_wt or '-'}->{res.outer_product_ko or '-'}"
              f"{shift}{ctrl}{kept}")
        for w in res.warnings:
            print(f"          - {w}")

    outcome = run_batch(entries, cfg, genome, jobs=jobs,
                        on_locate=on_locate, on_result=on_result)
    results, failures = outcome.results, outcome.failures

    if results:
        write_results_csv(results, args.out)
        print(f"\nwrote {args.out}  ({len(results)} gene(s), {len(results) * 4} primers)")
    if failures:
        fail_path = args.out.with_name(args.out.stem + "_failed.csv")
        write_failures_csv(failures, fail_path)
        print(f"wrote {fail_path}  ({len(failures)} gene(s) could not be designed)")
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps({
            "criteria": {k: v for k, v in cfg.items() if not k.startswith("_")},
            "n_genes": len(entries),
            "n_ok": sum(1 for r in results if r.ok),
            "results": [r.as_dict() for r in results],
            "failures": [{"gene": g, "reason": why} for g, why in failures],
        }, indent=1, default=str), encoding="utf-8")
        print(f"wrote {args.report}")

    n_ok = sum(1 for r in results if r.ok)
    n_check = len(results) - n_ok
    print(f"\n{n_ok} ready, {n_check} to check, {len(failures)} failed")
    if n_check:
        print("'to check' means a warning was raised, not that the primers are unusable "
              "— read the Notes column.")
    return 1 if failures and not results else 0


if __name__ == "__main__":
    sys.exit(main())
