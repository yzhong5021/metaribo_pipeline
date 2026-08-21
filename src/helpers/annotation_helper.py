"""CLI helper functions for metagenome annotation"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

from utils import (
    WorkflowError,
    ensure_dir,
    fail,
    load_manifest,
    load_yaml,
    master_config,
    parse_shell_args,
    read_contig_to_genus,
    read_samplesheet,
    require_keys,
    require_path,
    resolve_command,
    save_manifest,
    shell_assignments,
    write_contig_to_genus,
    write_text,
)

# for cat output data matching
GENUS_RE = re.compile(r"([^\t;]+?)\s+\(genus\):\s*([0-9]+(?:\.[0-9]+)?)")
ANNOTATION_CONFIG_KEYS = (
    "cat_params",
    "cat_db",
    "cat_tax",
    "bakta_db",
    "bakta_params",
    "kofam_profiles",
    "kofam_list",
    "threads",
)


def load_annotation_config(path: Path) -> dict:
    cfg = load_yaml(path)
    require_keys(cfg, ANNOTATION_CONFIG_KEYS, path)
    for key in ("cat_db", "cat_tax", "bakta_db", "kofam_profiles", "kofam_list"):
        raw_value = cfg[key]
        value = "" if raw_value is None else str(raw_value).strip()
        if not value:
            if key in {"cat_db", "cat_tax"}:
                cfg[key] = None
                continue
            fail(f"{path} parameter {key} must be a path")
        cfg[key] = str(Path(value).expanduser())
    if (
        isinstance(cfg["threads"], bool)
        or not isinstance(cfg["threads"], int)
        or cfg["threads"] < 1
    ):
        fail(f"{path} parameter threads must be a positive integer")
    return cfg


def cat_is_enabled(cfg: dict) -> bool:
    if (cfg["cat_db"] is None) != (cfg["cat_tax"] is None):
        fail("cat_db and cat_tax must both be configured or both be null")
    return cfg["cat_db"] is not None


def validate_annotation_inputs(
    manifest_path: Path,
    samples: list[dict],
    check_files: bool,
) -> dict:
    metagenome = load_manifest(manifest_path).get("metagenome", {})
    for key in ("contigs", "bowtie2"):
        values = metagenome.get(key)
        if not isinstance(values, dict):
            fail(f"Manifest metagenome.{key} is missing; run preprocessing collection first")
        missing = [row["sample"] for row in samples if row["sample"] not in values]
        if missing:
            fail(
                f"Manifest metagenome.{key} is missing samples: "
                f"{', '.join(missing)}"
            )
        if check_files:
            for row in samples:
                sample = row["sample"]
                require_path(values[sample], f"{sample} {key}", "file")
    genus_paths = metagenome.get("contig_to_genus", {})
    if not isinstance(genus_paths, dict):
        fail("Manifest metagenome.contig_to_genus must be a JSON object")
    if check_files:
        for row in samples:
            sample = row["sample"]
            if sample in genus_paths:
                require_path(
                    genus_paths[sample],
                    f"{sample} contig-to-genus table",
                    "file",
                )
    return metagenome


def check_paths(manifest_path: Path, samplesheet: Path, annotation_config: Path) -> None:
    cfg = load_annotation_config(annotation_config)
    samples = read_samplesheet(samplesheet, check_files=False)
    if cat_is_enabled(cfg):
        resolve_command("CAT_pack")
        require_path(cfg["cat_db"], "CAT database", "dir")
        require_path(cfg["cat_tax"], "CAT taxonomy", "dir")
    else:
        print("CAT classification disabled (cat_db and cat_tax are null)")
    for command in ("bakta", "bedtools", "exec_annotation"):
        resolve_command(command)
    require_path(cfg["bakta_db"], "Bakta database", "dir")
    require_path(cfg["kofam_profiles"], "Kofam profiles", "dir")
    require_path(cfg["kofam_list"], "Kofam ko_list", "file")
    validate_annotation_inputs(manifest_path, samples, check_files=True)
    print("paths are all valid for metagenome annotation")


def config_shell(annotation_config: Path, master: dict) -> str:
    cfg = load_annotation_config(annotation_config)
    run_cat = cat_is_enabled(cfg)
    values = {
        "RUN_CAT": int(run_cat),
        "CAT_PATH": resolve_command("CAT_pack") if run_cat else "",
        "CAT_PARAM_ARGS": parse_shell_args(
            cfg["cat_params"], "cat_params", annotation_config
        ),
        "CAT_DB": cfg["cat_db"] or "",
        "CAT_TAX": cfg["cat_tax"] or "",
        "BAKTA": resolve_command("bakta"),
        "BAKTA_DB": cfg["bakta_db"],
        "BAKTA_PARAM_ARGS": parse_shell_args(
            cfg["bakta_params"], "bakta_params", annotation_config
        ),
        "BEDTOOLS": resolve_command("bedtools"),
        "EXEC_ANNOTATION": resolve_command("exec_annotation"),
        "KOFAM_PROFILES": cfg["kofam_profiles"],
        "KOFAM_LIST": cfg["kofam_list"],
        "THREADS": cfg["threads"],
        "CAT_DIR": master["annotation_scratch"] / "CAT",
        "BAKTA_DIR": master["annotation_scratch"] / "bakta",
        "ANNOT_WORKDIR": master["annotation_work"],
        "LOCAL_OUTDIR": master["metagenome_results"],
        "CONTIG_TO_GENUS_DIR": master["metagenome_results"] / "contig_to_genus",
        "BAKTA_CDS_GFF_DIR": master["metagenome_results"] / "bakta_cds_annotations_gff",
        "BAKTA_NCRNA_BED_DIR": master["metagenome_results"] / "bakta_ncrna_bed",
        "BAKTA_CDS_BED_DIR": master["metagenome_results"] / "bakta_cds_bed",
        "BAKTA_CDS_FAA_DIR": master["metagenome_results"] / "bakta_cds_faa",
        "BAKTA_COVERAGE_DIR": master["metagenome_results"] / "bakta_gene_mean_coverage",
        "KOFAM_DIR": master["metagenome_results"] / "kofam_annotations",
        "MASTER_SHEET_DIR": master["metagenome_results"] / "master_sheets",
    }
    return shell_assignments(values)


def list_inputs(manifest_path: Path, samplesheet: Path) -> None:
    samples = read_samplesheet(samplesheet, check_files=False)
    metagenome = validate_annotation_inputs(
        manifest_path, samples, check_files=False
    )
    print("sample_id\tcontigs\tcontig_to_genus\tbam")
    for row in samples:
        sample = row["sample"]
        genus_path = metagenome.get("contig_to_genus", {}).get(sample, "NA")
        print(
            f"{sample}\t{metagenome['contigs'][sample]}\t"
            f"{genus_path}\t{metagenome['bowtie2'][sample]}"
        )


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


def merge_cat_assignments(
    contig_to_genus_path: Path,
    cat_translated_path: Path,
    reset_existing: bool = False,
) -> int:
    # Keep MAG-propagated genus calls, and let CAT fill contigs still missing a genus.
    existing = (
        {} if reset_existing else read_contig_to_genus(contig_to_genus_path)
    )
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
            if len(fields) < 3 or fields[0] != "*":
                continue
            cds_id = fields[1].strip()
            if not cds_id:
                continue
            ko = next(
                (
                    field.strip()
                    for field in fields[2:]
                    if ko_re.fullmatch(field.strip())
                ),
                "",
            )
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


def update_manifest(
    manifest_path: Path,
    sample: str,
    cds_gff: Path,
    ncrna_bed: Path,
    kofam_tsv: Path,
    master_sheet_tsv: Path,
    contig_to_genus_path: Path | None = None,
) -> None:
    # Add the annotation files the ribo step and analysis notebooks read later.
    manifest = load_manifest(manifest_path)
    metagenome = manifest.setdefault("metagenome", {})
    outputs = {
        "bakta_cds_annotations_gff": cds_gff,
        "bakta_ncrna_bed": ncrna_bed,
        "kofam_annotations": kofam_tsv,
        "master_sheets": master_sheet_tsv,
    }
    if contig_to_genus_path is not None:
        outputs["contig_to_genus"] = contig_to_genus_path
    for key, path in outputs.items():
        entries = metagenome.setdefault(key, {})
        if not isinstance(entries, dict):
            fail(f"Manifest metagenome.{key} must be a JSON object")
        entries[sample] = str(path)
    save_manifest(manifest_path, manifest)


def main() -> int:
    parser = argparse.ArgumentParser(description="metagenome annotation helper")
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("check-paths")
    p.add_argument("--manifest", type=Path, default=None)
    p.add_argument("--samplesheet", type=Path, default=None)
    p.add_argument("--annotation-config", type=Path, default=None)
    p = sub.add_parser("config-shell")
    p.add_argument("--annotation-config", type=Path, default=None)
    p = sub.add_parser("list-inputs")
    p.add_argument("--manifest", type=Path, default=None)
    p.add_argument("--samplesheet", type=Path, default=None)
    p.add_argument("--annotation-config", type=Path, default=None)
    p = sub.add_parser("merge-cat")
    p.add_argument("--contig-to-genus", type=Path, required=True)
    p.add_argument("--cat-translated", type=Path, required=True)
    p.add_argument("--reset-existing", action="store_true")
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
    p.add_argument("--manifest", type=Path, default=None)
    p.add_argument("--sample", required=True)
    p.add_argument("--cds-gff", type=Path, required=True)
    p.add_argument("--ncrna-bed", type=Path, required=True)
    p.add_argument("--kofam", type=Path, required=True)
    p.add_argument("--master-sheet", type=Path, required=True)
    p.add_argument("--contig-to-genus", type=Path)
    p = sub.add_parser("sort-bed-for-coverage")
    p.add_argument("--fasta", type=Path, required=True)
    p.add_argument("--bed", type=Path, required=True)
    p.add_argument("--sorted-bed", type=Path, required=True)
    p.add_argument("--genome", type=Path, required=True)
    args = parser.parse_args()
    try:
        if args.command in {"check-paths", "config-shell", "list-inputs"}:
            master = master_config()
            annotation_config = (
                args.annotation_config or master["contig_annotation_params"]
            )
            manifest = getattr(args, "manifest", None) or master["file_manifest"]
            samplesheet = (
                getattr(args, "samplesheet", None)
                or master["root"] / "config" / "samples.csv"
            )
        if args.command == "check-paths":
            check_paths(manifest, samplesheet, annotation_config)
        elif args.command == "config-shell":
            print(config_shell(annotation_config, master))
        elif args.command == "list-inputs":
            list_inputs(manifest, samplesheet)
        elif args.command == "merge-cat":
            print(
                merge_cat_assignments(
                    args.contig_to_genus,
                    args.cat_translated,
                    args.reset_existing,
                )
            )
        elif args.command == "process-bakta":
            cds_gff, ncrna_bed, cds_bed, cds_faa, num_cds, num_ncrna = process_bakta_outputs(args.sample, args.bakta_run_dir, args.bakta_outdir, args.local_outdir)
            if args.shell:
                print(shell_assignments({"CDS_GFF": cds_gff, "NCRNA_BED": ncrna_bed, "CDS_BED": cds_bed, "CDS_FAA": cds_faa, "NUM_CDS": num_cds, "NUM_NCRNA": num_ncrna}))
            else:
                print(f"CDS_COUNT={num_cds}\nNCRNA_COUNT={num_ncrna}")
        elif args.command in {"make-master-sheet", "make-summary"}:
            make_master_sheet(args.cds_gff, args.contig_to_genus, args.coverage, args.kofam, args.output)
        elif args.command == "update-manifest":
            manifest = args.manifest or master_config()["file_manifest"]
            update_manifest(
                manifest,
                args.sample,
                args.cds_gff,
                args.ncrna_bed,
                args.kofam,
                args.master_sheet,
                args.contig_to_genus,
            )
        elif args.command == "sort-bed-for-coverage":
            sort_bed_for_coverage(args.fasta, args.bed, args.sorted_bed, args.genome)
    except WorkflowError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
