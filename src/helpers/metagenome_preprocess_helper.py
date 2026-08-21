"""Command-line helper for metagenome preprocessing workflow stages."""

from __future__ import annotations

import argparse
import csv
import gzip
import shutil
import sys
from pathlib import Path

from utils import (
    WorkflowError,
    copy_bam_index,
    copy_if_needed,
    ensure_dir,
    fail,
    load_manifest,
    load_yaml,
    master_config,
    read_samplesheet,
    require_keys,
    require_path,
    resolve_command,
    resolve_repo_path,
    resolve_scratch_path,
    save_manifest,
    shell_assignments,
    write_contig_to_genus,
)


NEXTFLOW_COMMAND = "nextflow"
METAGENOME_CONFIG_KEYS = (
    "params_file",
    "nextflow_config",
    "samplesheet",
    "gtdb_db",
    "apptainer_cache",
)


def load_metagenome_config(path: Path) -> dict:
    cfg = load_yaml(path)
    require_keys(cfg, METAGENOME_CONFIG_KEYS, path)
    return cfg


def gtdb_db_path(cfg: dict) -> str:
    raw_gtdb_db = cfg.get("gtdb_db")
    if raw_gtdb_db is None:
        return ""
    value = str(raw_gtdb_db).strip()
    return str(Path(value).expanduser()) if value else ""


def check_paths(metagenome_preprocess_config: Path, master: dict) -> None:
    cfg = load_metagenome_config(metagenome_preprocess_config)
    gtdb_db = gtdb_db_path(cfg)
    if gtdb_db:
        require_path(gtdb_db, "GTDB-Tk database", "dir")
    require_path(
        resolve_repo_path(cfg["params_file"], master),
        "nf-core/mag parameters",
        "file",
    )
    require_path(
        resolve_repo_path(cfg["nextflow_config"], master),
        "Nextflow configuration",
        "file",
    )
    ensure_dir(resolve_scratch_path(cfg["apptainer_cache"], master))
    ensure_dir(master["metagenome_work"])
    rows = read_samplesheet(resolve_repo_path(cfg["samplesheet"], master), check_files=True)
    print(f"paths are all valid for {len(rows)} nf-core/mag samples")


def command_shell(metagenome_preprocess_config: Path, master: dict) -> str:
    cfg = load_metagenome_config(metagenome_preprocess_config)
    values = {
        "NEXTFLOW": resolve_command(NEXTFLOW_COMMAND),
        "PARAMS_FILE": resolve_repo_path(cfg["params_file"], master),
        "NEXTFLOW_CONFIG": resolve_repo_path(cfg["nextflow_config"], master),
        "SAMPLESHEET": resolve_repo_path(cfg["samplesheet"], master),
        "OUTDIR": master["metagenome_scratch"],
        "GTDB_DB": gtdb_db_path(cfg),
        "APPTAINER_CACHEDIR": resolve_scratch_path(cfg["apptainer_cache"], master),
        "WORK_DIR": master["metagenome_work"],
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
        required_columns = {"user_genome", "classification"}
        if not required_columns.issubset(reader.fieldnames or []):
            fail(f"GTDB-Tk summary lacks required columns: {path}")
        for row in reader:
            genome = strip_fasta_suffix(row.get("user_genome") or "")
            if genome not in allowed_bins:
                continue
            for raw_part in (row.get("classification") or "").split(";"):
                part = raw_part.strip()
                if part.startswith("g__") and len(part) > 3:
                    bin_to_genus[genome] = part[3:]
                    break
    return bin_to_genus


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


def collect_outputs(
    manifest_path: Path,
    samplesheet: Path,
    outdir: Path,
    output_dir: Path,
    master: dict,
    metagenome_preprocess_config: Path,
) -> None:
    # After nf-core/mag finishes, keep the contigs, BAMs, bins, and genus calls the rest of the workflow needs.
    output_dir = resolve_repo_path(output_dir, master)
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
        if not isinstance(metagenome[key], dict):
            fail(f"Manifest metagenome.{key} must be a JSON object")

    gtdb_enabled = bool(
        gtdb_db_path(load_metagenome_config(metagenome_preprocess_config))
    )
    gtdb_summary = outdir / "Taxonomy" / "GTDB-Tk" / "gtdbtk_summary.tsv"
    if gtdb_enabled:
        require_path(gtdb_summary, "GTDB-Tk classification summary", "file")
        metagenome.setdefault("contig_to_genus", {})
        if not isinstance(metagenome["contig_to_genus"], dict):
            fail("Manifest metagenome.contig_to_genus must be a JSON object")
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
        source_bins = sorted(
            (dastool_dir / "bins").glob(f"MEGAHIT-*Refined-{sample}.*.fa")
        )
        bins = [copy_if_needed(path, bins_dir / path.name) for path in source_bins]
        if gtdb_enabled:
            if not source_bins:
                fail(f"No DASTool bins found for {sample}")
            bin_to_genus = parse_gtdb_summary(
                gtdb_summary,
                {strip_fasta_suffix(path.name) for path in source_bins},
            )
            contig2bin = find_one(sorted(dastool_dir.glob(f"MEGAHIT-*{sample}*_DASTool_contig2bin.tsv")), f"{sample} DASTool contig2bin TSV")
            contig_to_genus: dict[str, str] = {}
            with contig2bin.open() as handle:
                for line in handle:
                    if not line.strip():
                        continue
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 2:
                        fail(f"Malformed DASTool contig2bin row: {line.rstrip()}")
                    genus = bin_to_genus.get(strip_fasta_suffix(fields[1]))
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
    p.add_argument("--samplesheet", type=Path, default=None)
    p = sub.add_parser("check-paths")
    p.add_argument("--metagenome-preprocess-config", type=Path, default=None)
    p = sub.add_parser("command-shell")
    p.add_argument("--metagenome-preprocess-config", type=Path, default=None)
    p = sub.add_parser("collect-outputs")
    p.add_argument("--manifest", type=Path, default=None)
    p.add_argument("--samplesheet", type=Path, required=True)
    p.add_argument("--metagenome-preprocess-config", type=Path, default=None)
    args = parser.parse_args()
    try:
        master = master_config()
        metagenome_config = getattr(args, "metagenome_preprocess_config", None) or master["metagenome_preprocess_params"]
        if args.command == "validate-samples":
            samplesheet = args.samplesheet or master["root"] / "config" / "samples.csv"
            rows = read_samplesheet(samplesheet, check_files=True)
            print(f"Validated {samplesheet} with {len(rows)} samples")
        elif args.command == "check-paths":
            check_paths(metagenome_config, master)
        elif args.command == "command-shell":
            print(command_shell(metagenome_config, master))
        elif args.command == "collect-outputs":
            collect_outputs(
                args.manifest or master["file_manifest"],
                args.samplesheet,
                master["metagenome_scratch"],
                master["metagenome_results"],
                master,
                metagenome_config,
            )
    except WorkflowError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
