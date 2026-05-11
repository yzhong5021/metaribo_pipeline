#!/usr/bin/env python3
"""Command-line helper for metagenome annotation workflow stages."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

from utils import (
    CONTIG_ANNOTATION_PARAMS,
    FILE_MANIFEST,
    ROOT,
    WorkflowError,
    ensure_dir,
    fail,
    load_manifest,
    load_yaml,
    read_contig_to_genus,
    read_samplesheet,
    require_path,
    resolve_command,
    resolve_repo_path,
    resolve_scratch_path,
    resolve_work_path,
    save_manifest,
    shell_assignments,
    write_contig_to_genus,
    write_text,
)


GENUS_RE = re.compile(r"([^\t;]+?)\s+\(genus\):\s*([0-9]+(?:\.[0-9]+)?)")
ANNOTATION_COMMANDS = {
    "cat_pack": "CAT_pack",
    "bakta": "bakta",
    "bedtools": "bedtools",
    "exec_annotation": "exec_annotation",
}


def annotation_local_outdir(cfg: dict) -> Path:
    return resolve_repo_path(cfg.get("local_outdir", ROOT / "output" / "metagenome"))


def command_name(config: dict, key: str) -> str:
    commands = config.get("commands", {})
    return str(commands.get(key, ANNOTATION_COMMANDS[key]))


def check_paths(manifest_path: Path, samplesheet: Path, annotation_config: Path, version_config: Path | None = None) -> None:
    cfg = load_yaml(annotation_config)
    samples = read_samplesheet(samplesheet, check_files=False)
    require_path(cfg["cat_db"], "CAT database folder", "dir")
    require_path(cfg["cat_tax"], "CAT taxonomy folder", "dir")
    for command_key in ANNOTATION_COMMANDS:
        resolve_command(command_name(cfg, command_key))
    require_path(cfg["bakta_db"], "Bakta database", "dir")
    require_path(cfg["kofam_profiles"], "Kofam profiles", "dir")
    require_path(cfg["kofam_list"], "Kofam ko_list", "file")
    metagenome = load_manifest(manifest_path).get("metagenome", {})
    required_keys = ["contigs", "bowtie2"]
    for key in required_keys:
        if key not in metagenome:
            fail(f"Manifest metagenome.{key} is missing; run preprocessing collection first")
        missing = [row["sample"] for row in samples if row["sample"] not in metagenome[key]]
        if missing:
            fail(f"Manifest metagenome.{key} is missing samples from {samplesheet}: {', '.join(missing)}")
    print("paths are all valid for metagenome annotation")


def config_shell(annotation_config: Path, version_config: Path | None = None) -> str:
    cfg = load_yaml(annotation_config)
    outdir = resolve_scratch_path(cfg["outdir"])
    workdir = resolve_work_path("metagenome_annotation")
    local_outdir = annotation_local_outdir(cfg)
    values = {
        "CAT_PATH": resolve_command(command_name(cfg, "cat_pack")),
        "CAT_PARAMS": cfg.get("cat_params", ""),
        "CAT_DB": cfg.get("cat_db", ""),
        "CAT_TAX": cfg.get("cat_tax", ""),
        "BAKTA": resolve_command(command_name(cfg, "bakta")),
        "BAKTA_DB": cfg["bakta_db"],
        "BAKTA_PARAMS": cfg.get("bakta_params", ""),
        "BEDTOOLS": resolve_command(command_name(cfg, "bedtools")),
        "EXEC_ANNOTATION": resolve_command(command_name(cfg, "exec_annotation")),
        "KOFAM_PROFILES": cfg["kofam_profiles"],
        "KOFAM_LIST": cfg["kofam_list"],
        "THREADS": cfg["threads"],
        "ANNOT_OUTDIR": outdir,
        "CAT_DIR": outdir / "CAT",
        "BAKTA_DIR": outdir / "bakta",
        "ANNOT_WORKDIR": workdir,
        "LOCAL_OUTDIR": local_outdir,
        "CONTIG_TO_GENUS_DIR": local_outdir / "contig_to_genus",
        "BAKTA_CDS_GFF_DIR": local_outdir / "bakta_cds_annotations_gff",
        "BAKTA_NCRNA_BED_DIR": local_outdir / "bakta_ncrna_bed",
        "BAKTA_CDS_BED_DIR": local_outdir / "bakta_cds_bed",
        "BAKTA_CDS_FAA_DIR": local_outdir / "bakta_cds_faa",
        "BAKTA_COVERAGE_DIR": local_outdir / "bakta_gene_mean_coverage",
        "KOFAM_DIR": local_outdir / "kofam_annotations",
        "MASTER_SHEET_DIR": local_outdir / "master_sheets",
    }
    return shell_assignments(values)


def list_inputs(manifest_path: Path, samplesheet: Path, annotation_config: Path) -> None:
    manifest = load_manifest(manifest_path)
    samples = read_samplesheet(samplesheet, check_files=False)
    metagenome = manifest["metagenome"]
    print("sample_id\tcontigs\tcontig_to_genus\tbam")
    for row in samples:
        sample = row["sample"]
        print(f"{sample}\t{metagenome['contigs'][sample]}\tNA\t{metagenome['bowtie2'][sample]}")


def parse_cat_translated(path: Path) -> dict[str, tuple[str, float]]:
    assignments: dict[str, tuple[str, float]] = {}
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#") or "(genus)" not in line:
                continue
            contig = line.split("\t", 1)[0]
            matches = [(m.group(1).strip(), float(m.group(2))) for m in GENUS_RE.finditer(line)]
            if matches:
                assignments[contig] = max(matches, key=lambda item: item[1])
    return assignments


def merge_cat_assignments(contig_to_genus_path: Path, cat_translated_path: Path) -> int:
    # Keep MAG-propagated genus calls, and let CAT fill contigs still missing a genus.
    existing = read_contig_to_genus(contig_to_genus_path)
    updated = 0
    for contig, (genus, _score) in parse_cat_translated(cat_translated_path).items():
        if contig not in existing:
            existing[contig] = genus
            updated += 1
    write_contig_to_genus(contig_to_genus_path, existing)
    return updated


def parse_gff_attributes(attributes: str) -> dict[str, str]:
    parsed = {}
    for item in attributes.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            parsed[key] = value
    return parsed


def cds_records(cds_gff: Path) -> list[dict]:
    records = []
    with cds_gff.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            attrs = parse_gff_attributes(fields[8])
            start = int(fields[3])
            end = int(fields[4])
            records.append(
                {
                    "CDS_id": attrs.get("ID") or attrs.get("Name") or f"CDS_{len(records) + 1}",
                    "contig_id": fields[0],
                    "start": start,
                    "end": end,
                    "strand": fields[6],
                    "gene_length": end - start + 1,
                }
            )
    return records


def process_bakta_outputs(sample: str, bakta_run_dir: Path, bakta_outdir: Path, local_outdir: Path | None = None) -> tuple[Path, Path, Path, Path, int, int]:
    # Pull out the Bakta CDS, ncRNA, BED, and protein files used by coverage and Kofam.
    source_gff = require_path(bakta_run_dir / f"{sample}.gff3", f"{sample} Bakta GFF3", "file")
    source_faa = require_path(bakta_run_dir / f"{sample}.faa", f"{sample} Bakta FAA", "file")
    if local_outdir is None:
        ensure_dir(bakta_outdir)
        cds_gff = bakta_outdir / f"{sample}_cds.gff3"
        ncrna_bed = bakta_outdir / f"{sample}_ncrna.bed"
        cds_bed = bakta_outdir / f"{sample}_cds.bed"
        cds_faa = bakta_outdir / f"{sample}_cds.faa"
    else:
        cds_gff = local_outdir / "bakta_cds_annotations_gff" / f"{sample}_cds.gff3"
        ncrna_bed = local_outdir / "bakta_ncrna_bed" / f"{sample}_ncrna.bed"
        cds_bed = local_outdir / "bakta_cds_bed" / f"{sample}_cds.bed"
        cds_faa = local_outdir / "bakta_cds_faa" / f"{sample}_cds.faa"
    cds_lines = ["##gff-version 3\n"]
    ncrna_rows = []
    with source_gff.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            start = int(fields[3])
            end = int(fields[4])
            attrs = parse_gff_attributes(fields[8])
            feature_id = attrs.get("ID") or attrs.get("Name") or fields[2]
            if fields[2] == "CDS":
                cds_lines.append(line)
            elif fields[2] in {"rRNA", "tRNA", "tmRNA"}:
                ncrna_rows.append((fields[0], start - 1, end, feature_id, "0", fields[6]))
    write_text(cds_gff, "".join(cds_lines))
    cds = cds_records(cds_gff)
    write_text(cds_bed, "".join("\t".join(map(str, row)) + "\n" for row in sorted((r["contig_id"], r["start"] - 1, r["end"], r["CDS_id"], "0", r["strand"]) for r in cds)))
    write_text(ncrna_bed, "".join("\t".join(map(str, row)) + "\n" for row in sorted(ncrna_rows)))
    write_text(cds_faa, source_faa.read_text())
    return cds_gff, ncrna_bed, cds_bed, cds_faa, len(cds), len(ncrna_rows)


def fasta_lengths(fasta: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    current: str | None = None
    with fasta.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                lengths[current] = 0
            elif current is not None:
                lengths[current] += len(line.strip())
    return lengths


def sort_bed_for_coverage(fasta: Path, bed: Path, sorted_bed: Path, genome_file: Path) -> None:
    lengths = fasta_lengths(fasta)
    order = {contig: idx for idx, contig in enumerate(lengths)}
    rows = []
    with bed.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[0] not in order:
                fail(f"BED contig {fields[0]} is absent from FASTA {fasta}")
            rows.append((order[fields[0]], int(fields[1]), int(fields[2]), fields))
    write_text(genome_file, "".join(f"{contig}\t{length}\n" for contig, length in lengths.items()))
    write_text(sorted_bed, "".join("\t".join(row[3]) + "\n" for row in sorted(rows)))


def parse_coverage(path: Path) -> dict[str, str]:
    coverage = {}
    if not path.exists():
        return coverage
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 7:
                coverage[fields[3]] = fields[-1]
    return coverage


def parse_kofam(path: Path) -> dict[str, str]:
    annotations: dict[str, list[str]] = {}
    ko_re = re.compile(r"^K\d{5}$")
    if not path.exists():
        return {}
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            cds_id = fields[1] if fields[0] == "*" and len(fields) > 2 else fields[0]
            search_fields = fields[2:] if fields[0] == "*" and len(fields) > 2 else fields[1:]
            ko = next((field for field in search_fields if ko_re.match(field)), "")
            if ko:
                annotations.setdefault(cds_id, []).append(ko)
    return {cds_id: ",".join(sorted(set(kos))) for cds_id, kos in annotations.items()}


def make_master_sheet(cds_gff: Path, contig_to_genus_path: Path | None, coverage_path: Path, kofam_path: Path, output_path: Path) -> None:
    # One row per Bakta CDS: coordinates, KO, genus, and metagenome coverage for notebook input.
    genus = read_contig_to_genus(contig_to_genus_path) if contig_to_genus_path else {}
    coverage = parse_coverage(coverage_path)
    kofam = parse_kofam(kofam_path)
    ensure_dir(output_path.parent)
    with output_path.open("w", newline="") as handle:
        fieldnames = ["CDS_id", "contig_id", "start", "end", "strand", "gene_length", "kofam_ko", "genus", "coverage"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for record in cds_records(cds_gff):
            writer.writerow(
                {
                    "CDS_id": record["CDS_id"],
                    "contig_id": record["contig_id"],
                    "start": record["start"],
                    "end": record["end"],
                    "strand": record["strand"],
                    "gene_length": record["gene_length"],
                    "kofam_ko": kofam.get(record["CDS_id"], ""),
                    "genus": genus.get(record["contig_id"], ""),
                    "coverage": coverage.get(record["CDS_id"], ""),
                }
            )


make_annotation_summary = make_master_sheet


def update_manifest(manifest_path: Path, sample: str, cds_gff: Path, ncrna_bed: Path, kofam_tsv: Path, master_sheet_tsv: Path) -> None:
    # Add the annotation files the ribo step and analysis notebooks read later.
    manifest = load_manifest(manifest_path)
    metagenome = manifest.setdefault("metagenome", {})
    for key in ("bakta_cds_annotations_gff", "bakta_ncrna_bed", "kofam_annotations", "master_sheets"):
        metagenome.setdefault(key, {})
    metagenome["bakta_cds_annotations_gff"][sample] = str(cds_gff)
    metagenome["bakta_ncrna_bed"][sample] = str(ncrna_bed)
    metagenome["kofam_annotations"][sample] = str(kofam_tsv)
    metagenome["master_sheets"][sample] = str(master_sheet_tsv)
    save_manifest(manifest_path, manifest)


def main() -> int:
    parser = argparse.ArgumentParser(description="metagenome annotation helper")
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("check-paths")
    p.add_argument("--manifest", type=Path, default=FILE_MANIFEST)
    p.add_argument("--samplesheet", type=Path, default=ROOT / "config" / "samples.csv")
    p.add_argument("--annotation-config", type=Path, default=CONTIG_ANNOTATION_PARAMS)
    p.add_argument("--versions", type=Path, default=None)
    p = sub.add_parser("config-shell")
    p.add_argument("--annotation-config", type=Path, default=CONTIG_ANNOTATION_PARAMS)
    p.add_argument("--versions", type=Path, default=None)
    p = sub.add_parser("list-inputs")
    p.add_argument("--manifest", type=Path, default=FILE_MANIFEST)
    p.add_argument("--samplesheet", type=Path, default=ROOT / "config" / "samples.csv")
    p.add_argument("--annotation-config", type=Path, default=CONTIG_ANNOTATION_PARAMS)
    p = sub.add_parser("merge-cat")
    p.add_argument("--contig-to-genus", type=Path, required=True)
    p.add_argument("--cat-translated", type=Path, required=True)
    p = sub.add_parser("process-bakta")
    p.add_argument("--sample", required=True)
    p.add_argument("--bakta-run-dir", type=Path, required=True)
    p.add_argument("--bakta-outdir", type=Path, required=True)
    p.add_argument("--local-outdir", type=Path)
    p.add_argument("--shell", action="store_true")
    for command_name in ("make-master-sheet", "make-summary"):
        p = sub.add_parser(command_name)
        p.add_argument("--cds-gff", type=Path, required=True)
        p.add_argument("--contig-to-genus", type=Path)
        p.add_argument("--coverage", type=Path, required=True)
        p.add_argument("--kofam", type=Path, required=True)
        p.add_argument("--output", type=Path, required=True)
    p = sub.add_parser("update-manifest")
    p.add_argument("--manifest", type=Path, default=FILE_MANIFEST)
    p.add_argument("--sample", required=True)
    p.add_argument("--cds-gff", type=Path, required=True)
    p.add_argument("--ncrna-bed", type=Path, required=True)
    p.add_argument("--kofam", type=Path, required=True)
    p.add_argument("--master-sheet", type=Path, required=True)
    p = sub.add_parser("sort-bed-for-coverage")
    p.add_argument("--fasta", type=Path, required=True)
    p.add_argument("--bed", type=Path, required=True)
    p.add_argument("--sorted-bed", type=Path, required=True)
    p.add_argument("--genome", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.command == "check-paths":
            check_paths(args.manifest, args.samplesheet, args.annotation_config, args.versions)
        elif args.command == "config-shell":
            print(config_shell(args.annotation_config, args.versions))
        elif args.command == "list-inputs":
            list_inputs(args.manifest, args.samplesheet, args.annotation_config)
        elif args.command == "merge-cat":
            print(merge_cat_assignments(args.contig_to_genus, args.cat_translated))
        elif args.command == "process-bakta":
            cds_gff, ncrna_bed, cds_bed, cds_faa, num_cds, num_ncrna = process_bakta_outputs(args.sample, args.bakta_run_dir, args.bakta_outdir, args.local_outdir)
            if args.shell:
                print(shell_assignments({"CDS_GFF": cds_gff, "NCRNA_BED": ncrna_bed, "CDS_BED": cds_bed, "CDS_FAA": cds_faa, "NUM_CDS": num_cds, "NUM_NCRNA": num_ncrna}))
            else:
                print(f"CDS_COUNT={num_cds}\nNCRNA_COUNT={num_ncrna}")
        elif args.command in {"make-master-sheet", "make-summary"}:
            make_master_sheet(args.cds_gff, args.contig_to_genus, args.coverage, args.kofam, args.output)
        elif args.command == "update-manifest":
            update_manifest(args.manifest, args.sample, args.cds_gff, args.ncrna_bed, args.kofam, args.master_sheet)
        elif args.command == "sort-bed-for-coverage":
            sort_bed_for_coverage(args.fasta, args.bed, args.sorted_bed, args.genome)
    except WorkflowError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
