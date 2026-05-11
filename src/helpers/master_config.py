#!/usr/bin/env python3
"""Resolve master configuration paths for shell entrypoints."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import yaml


FALLBACK_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MASTER_CONFIG = FALLBACK_ROOT / "config" / "master.yaml"
STAGE_DEFAULTS = {
    "file_manifest": "file_manifest.json",
    "preprocess_params": "config/ribo_preprocess.yaml",
    "metagenome_preprocess_params": "config/metagenome_preprocess.yaml",
    "contig_annotation_params": "config/contig_annotation.yaml",
}


class MasterConfigError(RuntimeError):
    """Raised when master configuration cannot be resolved."""


def load_yaml(path: Path) -> dict:
    with path.open() as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise MasterConfigError(f"Expected mapping in YAML file: {path}")
    return data


def resolve_from_root(root: Path, value: object, default: str) -> Path:
    path = Path(str(value if value is not None else default)).expanduser()
    return path if path.is_absolute() else root / path


def require_key(config: dict, key: str, source: Path) -> object:
    if key not in config:
        raise MasterConfigError(f"{source} is missing required key: {key}")
    return config[key]


def load_master(master_config: Path) -> dict[str, Path]:
    cfg = load_yaml(master_config)
    if "root" not in cfg:
        raise MasterConfigError(f"{master_config} is missing required key: root")
    if "scratch" not in cfg:
        raise MasterConfigError(f"{master_config} is missing required key: scratch")

    root = Path(str(cfg["root"])).expanduser()
    scratch = Path(str(cfg["scratch"])).expanduser()
    values = {"root": root, "scratch": scratch}
    for key, default in STAGE_DEFAULTS.items():
        values[key] = resolve_from_root(root, cfg.get(key), default)
    return values


def metagenome_values(master_config: Path) -> list[Path]:
    master = load_master(master_config)
    stage_cfg = load_yaml(master["metagenome_preprocess_params"])
    log_path = resolve_from_root(
        master["root"],
        require_key(stage_cfg, "log", master["metagenome_preprocess_params"]),
        "",
    )
    return [
        master["root"],
        master["file_manifest"],
        master["preprocess_params"],
        master["metagenome_preprocess_params"],
        master["contig_annotation_params"],
        log_path,
    ]


def annotation_values(master_config: Path) -> list[Path]:
    master = load_master(master_config)
    metagenome_cfg = load_yaml(master["metagenome_preprocess_params"])
    annotation_cfg = load_yaml(master["contig_annotation_params"])
    samplesheet = resolve_from_root(
        master["root"],
        metagenome_cfg.get("samplesheet"),
        "config/samples.csv",
    )
    log_path = resolve_from_root(
        master["root"],
        require_key(annotation_cfg, "log", master["contig_annotation_params"]),
        "",
    )
    return [
        master["root"],
        master["file_manifest"],
        master["preprocess_params"],
        master["metagenome_preprocess_params"],
        master["contig_annotation_params"],
        samplesheet,
        log_path,
    ]


def ribo_values(master_config: Path) -> list[Path]:
    master = load_master(master_config)
    stage_cfg = load_yaml(master["preprocess_params"])
    log_path = resolve_from_root(
        master["root"],
        require_key(stage_cfg, "log", master["preprocess_params"]),
        "",
    )
    return [
        master["root"],
        master["file_manifest"],
        master["preprocess_params"],
        master["metagenome_preprocess_params"],
        master["contig_annotation_params"],
        log_path,
    ]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="resolve MetaRibo master config")
    parser.add_argument("stage", choices=("metagenome", "annotation", "ribo"))
    parser.add_argument("--master-config", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        if args.stage == "metagenome":
            values = metagenome_values(args.master_config)
        elif args.stage == "annotation":
            values = annotation_values(args.master_config)
        else:
            values = ribo_values(args.master_config)
    except (KeyError, OSError, yaml.YAMLError, MasterConfigError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    for value in values:
        print(value)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
