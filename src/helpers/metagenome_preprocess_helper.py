#!/usr/bin/env python3
"""Command-line helper for metagenome preprocessing workflow stages."""

from __future__ import annotations

import argparse
import csv
import gzip
import shutil
import sys
from pathlib import Path

from utils import (
    FILE_MANIFEST,
    METAGENOME_PREPROCESS_PARAMS,
    ROOT,
    WorkflowError,
    copy_bam_index,
    copy_if_needed,
    ensure_dir,
    fail,
    load_manifest,
    load_yaml,
    read_samplesheet,
    require_path,
    resolve_command,
    resolve_repo_path,
    resolve_scratch_path,
    resolve_work_path,
    save_manifest,
    shell_assignments,
    write_contig_to_genus,
)


NEXTFLOW_COMMAND = "nextflow"


def check_paths(metagenome_preprocess_config: Path, version_config: Path | None = None, check_fastq_integrity: bool = False) -> None:
    cfg = load_yaml(metagenome_preprocess_config)
    if cfg.get("gtdb_db"):
        require_path(cfg["gtdb_db"], "GTDB-Tk database", "dir")
    ensure_dir(resolve_scratch_path(cfg["apptainer_cache"]))
    ensure_dir(resolve_work_path(cfg.get("work_dir", "metagenome_preprocess")))
    rows = read_samplesheet(resolve_repo_path(cfg["samplesheet"]), check_files=True, check_integrity=check_fastq_integrity)
    print(f"paths are all valid for {len(rows)} nf-core/mag samples")


def command_shell(metagenome_preprocess_config: Path) -> str:
    cfg = load_yaml(metagenome_preprocess_config)
    local_outdir = resolve_repo_path(cfg.get("local_outdir", ROOT / "output" / "metagenome"))
    values = {
        "NEXTFLOW": resolve_command(NEXTFLOW_COMMAND),
        "PARAMS_FILE": resolve_repo_path(cfg["params_file"]),
        "NEXTFLOW_CONFIG": resolve_repo_path(cfg["nextflow_config"]),
        "SAMPLESHEET": resolve_repo_path(cfg["samplesheet"]),
        "OUTDIR": resolve_scratch_path(cfg["outdir"]),
        "LOCAL_OUTDIR": local_outdir,
        "GTDB_DB": cfg.get("gtdb_db", ""),
        "APPTAINER_CACHEDIR": resolve_scratch_path(cfg["apptainer_cache"]),
        "WORK_DIR": resolve_work_path(cfg.get("work_dir", "metagenome_preprocess")),
    }
    return shell_assignments(values)


def strip_fasta_suffix(path_or_name: str) -> str:
    name = Path(path_or_name).name
    for suffix in (".fa.gz", ".fna.gz", ".fasta.gz", ".fa", ".fna", ".fasta"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def parse_gtdb_summary(path: Path, allowed_bins: set[str]) -> dict[str, str]:
    bin_to_genus: dict[str, str] = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            genome = strip_fasta_suffix(row.get("user_genome", ""))
            if genome not in allowed_bins:
                continue
            for part in row.get("classification", "").split(";"):
                if part.startswith("g__") and len(part) > 3:
                    bin_to_genus[genome] = part[3:]
                    break
    return bin_to_genus


def resolve_bin_genus(bin_id: str, bin_to_genus: dict[str, str]) -> str | None:
    candidates = [strip_fasta_suffix(bin_id), Path(bin_id).name, strip_fasta_suffix(Path(bin_id).name)]
    for candidate in candidates:
        if candidate in bin_to_genus:
            return bin_to_genus[candidate]
    stripped = strip_fasta_suffix(bin_id)
    for known, genus in bin_to_genus.items():
        if known.endswith("." + stripped) or known.endswith("-" + stripped):
            return genus
    return None


def gunzip_to_fasta(gz_path: Path, fasta_path: Path | None = None) -> Path:
    if fasta_path is None:
        fasta_path = gz_path.with_suffix("")
    ensure_dir(fasta_path.parent)
    if fasta_path.exists() and fasta_path.stat().st_mtime >= gz_path.stat().st_mtime:
        return fasta_path
    with gzip.open(gz_path, "rt") as src, fasta_path.open("w") as dst:
        shutil.copyfileobj(src, dst)
    return fasta_path


def find_one(paths: list[Path], label: str) -> Path:
    if not paths:
        fail(f"No {label} found")
    return sorted(paths)[0]


def collect_outputs(manifest_path: Path, samplesheet: Path, outdir: Path, output_dir: Path) -> None:
    # After nf-core/mag finishes, keep the contigs, BAMs, bins, and genus calls the rest of the workflow needs.
    output_dir = resolve_repo_path(output_dir)
    contigs_dir = output_dir / "contigs"
    bins_dir = output_dir / "bins"
    bowtie2_dir = output_dir / "bowtie2"
    contig_to_genus_dir = output_dir / "contig_to_genus"
    for directory in (contigs_dir, bins_dir, bowtie2_dir):
        ensure_dir(directory)

    manifest = load_manifest(manifest_path)
    samples = read_samplesheet(samplesheet, check_files=False)
    metagenome = manifest.setdefault("metagenome", {})
    for key in ("contigs", "bins", "bowtie2"):
        metagenome.setdefault(key, {})

    gtdb_dir = outdir / "Taxonomy" / "GTDB-Tk"
    gtdb_summaries = [p for p in [gtdb_dir / "gtdbtk_summary.tsv", *sorted(gtdb_dir.glob("*summary.tsv"))] if p.exists()]
    gtdb_summary = sorted(gtdb_summaries)[0] if gtdb_summaries else None
    if gtdb_summary is not None:
        metagenome.setdefault("contig_to_genus", {})
        ensure_dir(contig_to_genus_dir)

    for row in samples:
        sample = row["sample"]
        contig_gz = require_path(outdir / "Assembly" / "MEGAHIT" / f"MEGAHIT-{sample}.contigs.fa.gz", f"{sample} MEGAHIT contigs", "file")
        contigs_fa = gunzip_to_fasta(contig_gz, contigs_dir / contig_gz.with_suffix("").name)
        qc_dir = outdir / "Assembly" / "MEGAHIT" / "QC" / sample
        source_bam = find_one(sorted(qc_dir.glob(f"MEGAHIT-{sample}-*.bam")), f"{sample} Bowtie2 BAM")
        bam = copy_if_needed(source_bam, bowtie2_dir / source_bam.name)
        copy_bam_index(source_bam, bam)
        dastool_dir = outdir / "GenomeBinning" / "DASTool"
        source_bins = sorted((dastool_dir / "bins").glob(f"MEGAHIT-*Refined-{sample}.*.fa")) or sorted((dastool_dir / "bins").glob(f"*{sample}*.fa"))
        bins = [copy_if_needed(path, bins_dir / path.name) for path in source_bins]
        if gtdb_summary is not None:
            if not source_bins:
                fail(f"No DASTool bins found for {sample}")
            bin_to_genus = parse_gtdb_summary(gtdb_summary, {strip_fasta_suffix(p.name) for p in source_bins})
            contig2bin = find_one(sorted(dastool_dir.glob(f"MEGAHIT-*{sample}*_DASTool_contig2bin.tsv")), f"{sample} DASTool contig2bin TSV")
            contig_to_genus: dict[str, str] = {}
            with contig2bin.open() as handle:
                for line in handle:
                    if not line.strip():
                        continue
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 2:
                        fail(f"Malformed DASTool contig2bin row: {line.rstrip()}")
                    genus = resolve_bin_genus(fields[1], bin_to_genus)
                    if genus:
                        contig_to_genus[fields[0]] = genus
            contig_to_genus_path = contig_to_genus_dir / f"contig_to_genus_{sample}.tsv"
            write_contig_to_genus(contig_to_genus_path, contig_to_genus)
            metagenome["contig_to_genus"][sample] = str(contig_to_genus_path)
        metagenome["contigs"][sample] = str(contigs_fa)
        metagenome["bins"][sample] = [str(p) for p in bins]
        metagenome["bowtie2"][sample] = str(bam)
        print(f"{sample} finished metagenome-preprocess.")
    save_manifest(manifest_path, manifest)


def main() -> int:
    parser = argparse.ArgumentParser(description="metagenome preprocessing helper")
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("validate-samples")
    p.add_argument("--samplesheet", type=Path, default=ROOT / "config" / "samples.csv")
    p.add_argument("--check-fastq-integrity", action="store_true")
    p = sub.add_parser("check-paths")
    p.add_argument("--metagenome-preprocess-config", type=Path, default=METAGENOME_PREPROCESS_PARAMS)
    p.add_argument("--versions", type=Path, default=None)
    p.add_argument("--check-fastq-integrity", action="store_true")
    p = sub.add_parser("command-shell")
    p.add_argument("--metagenome-preprocess-config", type=Path, default=METAGENOME_PREPROCESS_PARAMS)
    p = sub.add_parser("collect-outputs")
    p.add_argument("--manifest", type=Path, default=FILE_MANIFEST)
    p.add_argument("--samplesheet", type=Path, default=ROOT / "config" / "samples.csv")
    p.add_argument("--outdir", type=Path, default=None)
    p.add_argument("--output-dir", type=Path, default=None)
    p.add_argument("--metagenome-preprocess-config", type=Path, default=METAGENOME_PREPROCESS_PARAMS)
    args = parser.parse_args()
    try:
        if args.command == "validate-samples":
            rows = read_samplesheet(args.samplesheet, check_files=True, check_integrity=args.check_fastq_integrity)
            print(f"Validated {args.samplesheet} with {len(rows)} samples")
        elif args.command == "check-paths":
            check_paths(args.metagenome_preprocess_config, args.versions, check_fastq_integrity=args.check_fastq_integrity)
        elif args.command == "command-shell":
            print(command_shell(args.metagenome_preprocess_config))
        elif args.command == "collect-outputs":
            cfg = load_yaml(args.metagenome_preprocess_config)
            outdir = args.outdir if args.outdir is not None else resolve_scratch_path(cfg["outdir"])
            output_dir = args.output_dir if args.output_dir is not None else resolve_repo_path(cfg.get("local_outdir", ROOT / "output" / "metagenome"))
            collect_outputs(args.manifest, args.samplesheet, outdir, output_dir)
    except WorkflowError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
