"""Resolve master configuration paths for shell entrypoints."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import yaml


FALLBACK_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MASTER_CONFIG = FALLBACK_ROOT / "config" / "master.yaml"
MASTER_PATH_KEYS = (
    "results_dir",
    "logs_dir",
    "file_manifest",
    "ribo_preprocess_params",
    "metagenome_preprocess_params",
    "contig_annotation_params",
)


class MasterConfigError(RuntimeError):
    """Raised when master configuration cannot be resolved."""


def load_yaml(path: Path) -> dict:
    with path.open() as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise MasterConfigError(f"Expected mapping in YAML file: {path}")
    return data


def resolve_from_root(root: Path, value: object) -> Path:
    path = Path(str(value)).expanduser()
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
    missing = [key for key in MASTER_PATH_KEYS if key not in cfg]
    if missing:
        raise MasterConfigError(f"{master_config} is missing required keys: {', '.join(missing)}")

    root = Path(str(cfg["root"])).expanduser()
    scratch = Path(str(cfg["scratch"])).expanduser()
    results_dir = resolve_from_root(root, cfg["results_dir"])
    logs_dir = resolve_from_root(root, cfg["logs_dir"])
    values = {
        "root": root,
        "scratch": scratch,
        "results_dir": results_dir,
        "logs_dir": logs_dir,
        "file_manifest": resolve_from_root(root, cfg["file_manifest"]),
        "ribo_preprocess_params": resolve_from_root(root, cfg["ribo_preprocess_params"]),
        "metagenome_preprocess_params": resolve_from_root(root, cfg["metagenome_preprocess_params"]),
        "contig_annotation_params": resolve_from_root(root, cfg["contig_annotation_params"]),
        "metagenome_scratch": scratch / "metagenome_preprocess",
        "annotation_scratch": scratch / "metagenome_annotation",
        "ribo_scratch": scratch / "ribo_preprocess",
        "metagenome_work": scratch / "work" / "metagenome_preprocess",
        "annotation_work": scratch / "work" / "metagenome_annotation",
        "ribo_work": scratch / "work" / "ribo_preprocess",
        "metagenome_results": results_dir / "metagenome",
        "ribo_results": results_dir / "ribo",
        "coverage_results": results_dir / "coverage",
        "metagenome_log": logs_dir / "metagenome_preprocess.log",
        "annotation_log": logs_dir / "contig_annotation.log",
        "ribo_log": logs_dir / "ribo_preprocess.log",
    }
    return values


def metagenome_values(master_config: Path) -> list[Path]:
    master = load_master(master_config)
    return [
        master["root"],
        master["file_manifest"],
        master["metagenome_preprocess_params"],
        master["metagenome_log"],
    ]


def annotation_values(master_config: Path) -> list[Path]:
    master = load_master(master_config)
    metagenome_cfg = load_yaml(master["metagenome_preprocess_params"])
    samplesheet = resolve_from_root(
        master["root"],
        require_key(metagenome_cfg, "samplesheet", master["metagenome_preprocess_params"]),
    )
    return [
        master["root"],
        master["file_manifest"],
        master["contig_annotation_params"],
        samplesheet,
        master["annotation_log"],
    ]


def ribo_values(master_config: Path) -> list[Path]:
    master = load_master(master_config)
    return [
        master["root"],
        master["file_manifest"],
        master["ribo_preprocess_params"],
        master["ribo_log"],
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
