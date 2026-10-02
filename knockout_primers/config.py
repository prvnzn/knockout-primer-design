"""Configuration loading.

Every threshold the designer uses comes from a YAML file. The defaults ship with
the package (``knockout_primers/data/criteria.yaml``); to change them, either pass
``--set section.key=value`` on the command line or take your own copy:

    knockout-primers --write-config my_criteria.yaml
    knockout-primers genes.txt --config my_criteria.yaml

Code never hard-codes a cut-off, so every design decision is visible in one file.
"""
from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any

import yaml

PACKAGE_ROOT = Path(__file__).resolve().parent
DATA_DIR = PACKAGE_ROOT / "data"
DEFAULT_CONFIG = DATA_DIR / "criteria.yaml"


def load_config(path: str | Path | None = None) -> dict:
    """Load a criteria YAML, defaulting to the one shipped with the package."""
    p = Path(path) if path else DEFAULT_CONFIG
    if not p.exists():
        raise FileNotFoundError(f"config file not found: {p}")
    with open(p, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh) or {}
    cfg["_config_path"] = str(p)
    cfg["_config_dir"] = str(p.resolve().parent)
    return cfg


def resolve_path(rel: str | Path, cfg: dict | None = None) -> Path:
    """Resolve a path from the config.

    Absolute paths are used as given; a relative one is read next to the config
    file itself, so a copied-and-edited config keeps working wherever it lives.
    """
    p = Path(rel)
    if p.is_absolute():
        return p
    base = Path(cfg["_config_dir"]) if cfg and cfg.get("_config_dir") else DATA_DIR
    candidate = base / p
    return candidate if candidate.exists() else DATA_DIR / p


def write_default_config(dest: str | Path) -> Path:
    """Copy the default criteria file so it can be edited. Never overwrites."""
    out = Path(dest)
    if out.exists():
        raise FileExistsError(f"{out} already exists")
    out.parent.mkdir(parents=True, exist_ok=True)
    text = DEFAULT_CONFIG.read_text(encoding="utf-8")
    # A relative cassette path is read next to the config file, so point the copy
    # at the packaged cassette — otherwise it would break as soon as it is moved.
    cassette = (DATA_DIR / "FRT_aphA1_FRT_pKD13.fasta").as_posix()
    text = text.replace(
        'sequence_fasta:  "FRT_aphA1_FRT_pKD13.fasta"   # relative to the config file',
        f'sequence_fasta:  "{cassette}"')
    out.write_text(text, encoding="utf-8")
    return out


def _coerce(text: str) -> Any:
    """Turn a command-line string into a bool / number / string."""
    low = text.strip().lower()
    if low in ("true", "yes"):
        return True
    if low in ("false", "no"):
        return False
    if low in ("none", "null", ""):
        return None
    try:
        return int(text)
    except ValueError:
        pass
    try:
        return float(text)
    except ValueError:
        return text


def apply_overrides(cfg: dict, overrides: list[str] | None) -> dict:
    """Apply ``dotted.key=value`` overrides to a copy of the config.

    >>> cfg = apply_overrides(cfg, ["primers.tm.opt=62", "outer.min_shift_bp=100"])

    A path that is not already in the config raises ``KeyError``, so a typo is
    reported rather than silently creating a setting that nothing reads.
    """
    if not overrides:
        return cfg
    out = deepcopy(cfg)
    for item in overrides:
        if "=" not in item:
            raise ValueError(f"--set expects key=value, got {item!r}")
        dotted, raw = item.split("=", 1)
        keys = dotted.strip().split(".")
        node = out
        for depth, k in enumerate(keys[:-1]):
            if not isinstance(node.get(k), dict):
                raise KeyError(f"unknown config section {'.'.join(keys[:depth + 1])!r}")
            node = node[k]
        if keys[-1] not in node:
            raise KeyError(f"unknown config key {dotted!r}")
        node[keys[-1]] = _coerce(raw)
    return out
