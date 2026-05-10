#!/usr/bin/env python3
"""Command-line helper for Ribo-ITP preprocessing workflow stages."""

from __future__ import annotations

import argparse
import csv
import gzip
import shlex
import shutil
import sys
from pathlib import Path

from utils import (
    FILE_MANIFEST,
    PREPROCESS_PARAMS,
    ROOT,
    WorkflowError,
    copy_bam_index,
    copy_if_needed,
    ensure_dir,
    fail,
    load_manifest,
    load_yaml,
    require_path,
    resolve_repo_path,
    resolve_scratch_path,
    save_manifest,
    shell_assignments,
)


RIBO_SAMPLE_COLUMNS = ["id", "meta_id", "path"]
FASTQ_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz")


def as_bool(value: object, key: str) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        lowered = value.strip().lower()
        if lowered in {"true", "yes", "1"}:
            return True
        if lowered in {"false", "no", "0"}:
            return False
    fail(f"{key} must be true or false")


def read_ribo_samples(samplesheet_path: Path, check_files: bool = True) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    require_path(samplesheet_path, "samples_ribo.csv", "file")
    with samplesheet_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != RIBO_SAMPLE_COLUMNS:
            fail(f"{samplesheet_path} columns must be exactly: {','.join(RIBO_SAMPLE_COLUMNS)}")
        for row in reader:
            rows.append({col: (row.get(col) or "") for col in RIBO_SAMPLE_COLUMNS})
    if not rows:
        fail(f"{samplesheet_path} must contain at least one ribo sample row")
    seen: set[str] = set()
    for row in rows:
        sample = row["id"]
        for col in RIBO_SAMPLE_COLUMNS:
            if not row[col]:
                fail(f"Ribo sample {sample or '<blank>'} is missing required field '{col}'")
        if sample in seen:
            fail(f"Duplicate ribo id in {samplesheet_path}: {sample}")
        seen.add(sample)
        fastq = Path(row["path"]).expanduser()
        if not str(fastq).endswith(FASTQ_SUFFIXES):
            fail(f"Ribo sample {sample} path must be FASTQ/FASTQ.GZ: {fastq}")
        if check_files:
            require_path(fastq, f"Ribo sample {sample} FASTQ", "file")
        row["path"] = str(fastq)
    return rows


def ribo_config(path: Path) -> dict:
    cfg = load_yaml(path)
    required = [
        "samples",
        "cutadapt_args",
        "contains_umi",
        "umitools_extract_args",
        "bt2_args",
        "umitools_dedup_args",
        "remove_ncrna",
        "remove_multimappers",
        "outdir",
        "local_outdir",
        "log",
        "seed",
        "threads",
    ]
    missing = [key for key in required if key not in cfg]
    if missing:
        fail(f"{path} is missing required keys: {', '.join(missing)}")
    cfg["contains_umi"] = as_bool(cfg["contains_umi"], "contains_umi")
    cfg["remove_ncrna"] = as_bool(cfg["remove_ncrna"], "remove_ncrna")
    cfg["remove_multimappers"] = as_bool(cfg["remove_multimappers"], "remove_multimappers")
    return cfg


def samplesheet_from_config(cfg: dict) -> Path:
    return resolve_repo_path(cfg["samples"])


RIBO_COMMANDS = {
    "cutadapt": "cutadapt",
    "umi_tools": "umi_tools",
    "bowtie2": "bowtie2",
    "bowtie2_build": "bowtie2-build",
    "samtools": "samtools",
    "bedtools": "bedtools",
}


def command_name(config: dict, key: str) -> str:
    commands = config.get("commands", {})
    return str(commands.get(key, RIBO_COMMANDS[key]))


def validate_manifest_inputs(manifest_path: Path, samples: list[dict[str, str]], remove_ncrna: bool, check_paths: bool = True) -> None:
    manifest = load_manifest(manifest_path)
    metagenome = manifest.get("metagenome", {})
    if "contigs" not in metagenome:
        fail("Manifest metagenome.contigs is missing; run metagenome preprocessing collection first")
    if remove_ncrna and "bakta_ncrna_bed" not in metagenome:
        fail("Manifest metagenome.bakta_ncrna_bed is missing; run metagenome annotation before ncRNA filtering")
    for row in samples:
        ribo_id = row["id"]
        meta_id = row["meta_id"]
        if meta_id not in metagenome["contigs"]:
            fail(f"Ribo sample {ribo_id} meta_id {meta_id} is missing from manifest metagenome.contigs")
        if check_paths:
            require_path(metagenome["contigs"][meta_id], f"{ribo_id} metagenome contigs for {meta_id}", "file")
        if remove_ncrna:
            if meta_id not in metagenome["bakta_ncrna_bed"]:
                fail(f"Ribo sample {ribo_id} meta_id {meta_id} is missing from manifest metagenome.bakta_ncrna_bed")
            if check_paths:
                require_path(metagenome["bakta_ncrna_bed"][meta_id], f"{ribo_id} ncRNA BED for {meta_id}", "file")


def check_paths(manifest_path: Path, config_path: Path, version_config: Path | None = None) -> None:
    cfg = ribo_config(config_path)
    samples = read_ribo_samples(samplesheet_from_config(cfg), check_files=True)
    validate_manifest_inputs(manifest_path, samples, cfg["remove_ncrna"], check_paths=True)
    ensure_dir(resolve_scratch_path(cfg["outdir"]))
    ensure_dir(resolve_repo_path(cfg["local_outdir"]))
    print(f"paths are all valid for {len(samples)} ribo preprocessing samples")


def config_shell(config_path: Path, version_config: Path | None = None) -> str:
    cfg = ribo_config(config_path)
    outdir = resolve_scratch_path(cfg["outdir"])
    local_outdir = resolve_repo_path(cfg["local_outdir"])
    values = {
        "CUTADAPT": command_name(cfg, "cutadapt"),
        "UMI_TOOLS": command_name(cfg, "umi_tools"),
        "BOWTIE2": command_name(cfg, "bowtie2"),
        "BOWTIE2_BUILD": command_name(cfg, "bowtie2_build"),
        "SAMTOOLS": command_name(cfg, "samtools"),
        "BEDTOOLS": command_name(cfg, "bedtools"),
        "CUTADAPT_ARGS": cfg.get("cutadapt_args", ""),
        "UMITOOLS_EXTRACT_ARGS": cfg.get("umitools_extract_args", ""),
        "BT2_ARGS": cfg.get("bt2_args", ""),
        "UMITOOLS_DEDUP_ARGS": cfg.get("umitools_dedup_args", ""),
        "CONTAINS_UMI": int(cfg["contains_umi"]),
        "REMOVE_NCRNA": int(cfg["remove_ncrna"]),
        "REMOVE_MULTIMAPPERS": int(cfg["remove_multimappers"]),
        "RIBO_OUTDIR": outdir,
        "RIBO_DIR": outdir / "ribo",
        "BT2_DIR": outdir / "bt2",
        "REFERENCE_DIR": outdir / "references",
        "LOCAL_OUTDIR": local_outdir,
        "LOCAL_ALIGNMENT_DIR": local_outdir / "alignment",
        "LOCAL_STATS_DIR": local_outdir / "stats",
        "SEED": cfg["seed"],
        "THREADS": cfg["threads"],
    }
    return shell_assignments(values)


def list_inputs(manifest_path: Path, config_path: Path) -> None:
    cfg = ribo_config(config_path)
    samples = read_ribo_samples(samplesheet_from_config(cfg), check_files=False)
    validate_manifest_inputs(manifest_path, samples, cfg["remove_ncrna"], check_paths=False)
    metagenome = load_manifest(manifest_path)["metagenome"]
    print("ribo_id\tmeta_id\tribo_fq\tcontigs\tncrna_bed")
    for row in samples:
        meta_id = row["meta_id"]
        ncrna_bed = metagenome.get("bakta_ncrna_bed", {}).get(meta_id, "")
        print(f"{row['id']}\t{meta_id}\t{row['path']}\t{metagenome['contigs'][meta_id]}\t{ncrna_bed}")


def strip_gzip_suffix(path: Path) -> str:
    name = path.name
    return name[:-3] if name.endswith(".gz") else name


def reference_path(contigs: Path, ribo_id: str, reference_dir: Path) -> Path:
    contigs = Path(contigs)
    if contigs.name.endswith(".gz"):
        return reference_dir / ribo_id / strip_gzip_suffix(contigs)
    return contigs


def prepare_reference(contigs: Path, ribo_id: str, reference_dir: Path) -> Path:
    source = require_path(contigs, f"{ribo_id} metagenome contigs", "file")
    target = reference_path(source, ribo_id, reference_dir)
    if target == source:
        return source
    ensure_dir(target.parent)
    if target.exists() and target.stat().st_mtime >= source.stat().st_mtime:
        return target
    with gzip.open(source, "rt") as src, target.open("w") as dst:
        shutil.copyfileobj(src, dst)
    return target


def count_fastq(path: Path) -> int:
    require_path(path, "FASTQ", "file")
    opener = gzip.open if str(path).endswith(".gz") else open
    lines = 0
    with opener(path, "rt") as handle:
        for lines, _line in enumerate(handle, start=1):
            pass
    if lines % 4 != 0:
        fail(f"FASTQ line count is not divisible by four: {path}")
    return lines // 4


def init_stats(path: Path, sample: str) -> None:
    ensure_dir(path.parent)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "stage", "before_count", "after_count", "removed_count", "removed_proportion"])


def append_stats(path: Path, sample: str, stage: str, before: int, after: int) -> None:
    # Record how many reads each ribo preprocessing step removed.
    removed = before - after
    proportion = 0.0 if before == 0 else removed / before
    with path.open("a", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([sample, stage, before, after, removed, f"{proportion:.6f}"])


def publish_outputs(alignment_bam: Path, stats: Path, local_outdir: Path) -> tuple[Path, Path]:
    # Keep the final ribo BAM and stats in the local output tree for analysis.
    local_outdir = resolve_repo_path(local_outdir)
    target_bam = copy_if_needed(alignment_bam, local_outdir / "alignment" / alignment_bam.name)
    copy_bam_index(alignment_bam, target_bam)
    target_stats = copy_if_needed(stats, local_outdir / "stats" / stats.name)
    return target_bam, target_stats


def update_manifest(manifest_path: Path, sample: str, alignment_bam: Path, stats: Path, local_outdir: Path | None = None) -> None:
    # Add the final ribo alignment and stats file for downstream notebooks.
    if local_outdir is not None:
        alignment_bam, stats = publish_outputs(alignment_bam, stats, local_outdir)
    manifest = load_manifest(manifest_path)
    ribo = manifest.setdefault("ribo", {})
    ribo.setdefault("alignment", {})
    ribo.setdefault("stats", {})
    ribo["alignment"][sample] = str(alignment_bam)
    ribo["stats"][sample] = str(stats)
    save_manifest(manifest_path, manifest)


def split_args(args: str) -> None:
    for value in shlex.split(args or ""):
        print(value)


def main() -> int:
    parser = argparse.ArgumentParser(description="ribo preprocessing helper")
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("check-paths")
    p.add_argument("--manifest", type=Path, default=FILE_MANIFEST)
    p.add_argument("--ribo-config", type=Path, default=PREPROCESS_PARAMS)
    p.add_argument("--versions", type=Path, default=None)
    p = sub.add_parser("config-shell")
    p.add_argument("--ribo-config", type=Path, default=PREPROCESS_PARAMS)
    p.add_argument("--versions", type=Path, default=None)
    p = sub.add_parser("list-inputs")
    p.add_argument("--manifest", type=Path, default=FILE_MANIFEST)
    p.add_argument("--ribo-config", type=Path, default=PREPROCESS_PARAMS)
    p = sub.add_parser("reference-path")
    p.add_argument("--contigs", type=Path, required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--reference-dir", type=Path, required=True)
    p = sub.add_parser("prepare-reference")
    p.add_argument("--contigs", type=Path, required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--reference-dir", type=Path, required=True)
    p = sub.add_parser("count-fastq")
    p.add_argument("--fastq", type=Path, required=True)
    p = sub.add_parser("init-stats")
    p.add_argument("--stats", type=Path, required=True)
    p.add_argument("--sample", required=True)
    p = sub.add_parser("append-stats")
    p.add_argument("--stats", type=Path, required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--stage", required=True)
    p.add_argument("--before", type=int, required=True)
    p.add_argument("--after", type=int, required=True)
    p = sub.add_parser("update-manifest")
    p.add_argument("--manifest", type=Path, default=FILE_MANIFEST)
    p.add_argument("--sample", required=True)
    p.add_argument("--alignment-bam", type=Path, required=True)
    p.add_argument("--stats", type=Path, required=True)
    p.add_argument("--local-outdir", type=Path)
    p = sub.add_parser("split-args")
    p.add_argument("--args", default="")
    args = parser.parse_args()
    try:
        if args.command == "check-paths":
            check_paths(args.manifest, args.ribo_config, args.versions)
        elif args.command == "config-shell":
            print(config_shell(args.ribo_config, args.versions))
        elif args.command == "list-inputs":
            list_inputs(args.manifest, args.ribo_config)
        elif args.command == "reference-path":
            print(reference_path(args.contigs, args.sample, args.reference_dir))
        elif args.command == "prepare-reference":
            print(prepare_reference(args.contigs, args.sample, args.reference_dir))
        elif args.command == "count-fastq":
            print(count_fastq(args.fastq))
        elif args.command == "init-stats":
            init_stats(args.stats, args.sample)
        elif args.command == "append-stats":
            append_stats(args.stats, args.sample, args.stage, args.before, args.after)
        elif args.command == "update-manifest":
            update_manifest(args.manifest, args.sample, args.alignment_bam, args.stats, args.local_outdir)
        elif args.command == "split-args":
            split_args(args.args)
    except WorkflowError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
