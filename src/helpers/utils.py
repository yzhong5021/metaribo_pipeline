"""Shared helpers for MetaRibo pipeline command-line modules."""

from __future__ import annotations

import csv
import json
import os
import shlex
import shutil
from pathlib import Path
from typing import Sequence

from master_config import DEFAULT_MASTER_CONFIG, MasterConfigError, load_master

try:
    import yaml
except ImportError:
    yaml = None


SAMPLE_COLUMNS = [
    "sample",
    "group",
    "short_reads_1",
    "short_reads_2",
    "long_reads",
    "short_reads_platform",
    "long_reads_platform",
]


class WorkflowError(RuntimeError):
    """Raised for expected workflow validation failures."""


def fail(message: str) -> None:
    raise WorkflowError(message)


def load_yaml(path: Path) -> dict:
    if yaml is None:
        fail("PyYAML is required to read workflow YAML config files")
    try:
        with path.open() as handle:
            data = yaml.safe_load(handle) or {}
    except OSError as exc:
        fail(f"Cannot read YAML file {path}: {exc}")
    except yaml.YAMLError as exc:
        fail(f"Malformed YAML in {path}: {exc}")
    if not isinstance(data, dict):
        fail(f"Expected mapping in YAML file: {path}")
    return data


def default_master_config_path() -> Path:
    return Path(os.environ.get("MASTER_CONFIG", DEFAULT_MASTER_CONFIG)).expanduser()


def master_config(path: Path | None = None) -> dict:
    try:
        return load_master(path or default_master_config_path())
    except (MasterConfigError, OSError, yaml.YAMLError) as exc:
        fail(str(exc))


def require_keys(config: dict, keys: Sequence[str], source: Path) -> None:
    missing = [key for key in keys if key not in config]
    if missing:
        fail(f"{source} is missing required keys: {', '.join(missing)}")


def resolve_repo_path(value: str | Path, master: dict) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else master["root"] / path


def resolve_scratch_path(value: str | Path, master: dict) -> Path:
    path = Path(value).expanduser()
    if path.is_absolute():
        return path
    return master["scratch"] / path


def resolve_command(command: str) -> str:
    if "/" in command:
        return str(Path(command).expanduser())
    resolved = shutil.which(command)
    if resolved is None:
        fail(f"Command not found on PATH: {command}")
    return resolved


def parse_shell_args(value: object, key: str, source: Path) -> list[str]:
    if not isinstance(value, str):
        fail(f"{source} parameter {key} must be a string")
    try:
        return shlex.split(value)
    except ValueError as exc:
        fail(f"Malformed shell arguments for {key} in {source}: {exc}")


def ensure_dir(path: Path) -> None:
    if path.exists():
        if not path.is_dir():
            fail(f"Expected directory but found non-directory path: {path}")
        return
    path.mkdir(parents=True, exist_ok=True)


def copy_if_needed(source: Path, target: Path) -> Path:
    source = require_path(source, "source file", "file")
    ensure_dir(target.parent)
    if (
        target.exists()
        and target.stat().st_size == source.stat().st_size
        and target.stat().st_mtime >= source.stat().st_mtime
    ):
        return target
    shutil.copy2(source, target)
    return target


def copy_bam_index(source_bam: Path, target_bam: Path) -> None:
    candidates = [
        Path(str(source_bam) + ".bai"),
        source_bam.with_suffix(".bai"),
    ]
    source_index = next((path for path in candidates if path.exists()), None)
    if source_index is None:
        fail(f"No BAM index found for {source_bam}")
    copy_if_needed(source_index, Path(str(target_bam) + ".bai"))


def write_text(path: Path, text: str) -> None:
    ensure_dir(path.parent)
    path.write_text(text)


def load_manifest(path: Path) -> dict:
    # The manifest is only stage outputs; sample inputs stay in the config CSVs.
    if not path.exists():
        return {"metagenome": {}}
    try:
        text = path.read_text()
    except OSError as exc:
        fail(f"Cannot read manifest {path}: {exc}")
    if not text.strip():
        return {"metagenome": {}}
    try:
        data = json.loads(text)
    except json.JSONDecodeError as exc:
        fail(f"Malformed JSON in manifest {path}: {exc}")
    if not isinstance(data, dict):
        fail(f"Manifest must be a JSON object: {path}")
    if "samples" in data:
        fail(f"{path} must not contain sample inputs; config/samples.csv is the only sample entry point")
    data.setdefault("metagenome", {})
    if not isinstance(data["metagenome"], dict):
        fail(f"Manifest metagenome entry must be a JSON object: {path}")
    if "ribo" in data and not isinstance(data["ribo"], dict):
        fail(f"Manifest ribo entry must be a JSON object: {path}")
    return data


def save_manifest(path: Path, manifest: dict) -> None:
    # Write the handoff file used by the next workflow stage.
    ensure_dir(path.parent)
    path.write_text(json.dumps(manifest, indent=2, sort_keys=False) + "\n")


def require_path(path: str | Path, label: str, kind: str | None = None) -> Path:
    p = Path(path).expanduser()
    if not p.exists():
        fail(f"{label} does not exist: {p}")
    if kind == "file" and not p.is_file():
        fail(f"{label} is not a file: {p}")
    if kind == "dir" and not p.is_dir():
        fail(f"{label} is not a directory: {p}")
    return p


def read_samplesheet(samplesheet_path: Path, check_files: bool = True) -> list[dict]:
    rows: list[dict] = []
    require_path(samplesheet_path, "samples.csv", "file")
    try:
        with samplesheet_path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames != SAMPLE_COLUMNS:
                fail(f"{samplesheet_path} columns must be exactly: {','.join(SAMPLE_COLUMNS)}")
            for line_no, row in enumerate(reader, start=2):
                if None in row:
                    fail(f"Malformed row {line_no} in {samplesheet_path}: too many fields")
                rows.append({col: (row.get(col) or "").strip() for col in SAMPLE_COLUMNS})
    except csv.Error as exc:
        fail(f"Malformed CSV in {samplesheet_path}: {exc}")
    except OSError as exc:
        fail(f"Cannot read sample sheet {samplesheet_path}: {exc}")
    if not rows:
        fail(f"{samplesheet_path} must contain at least one sample row")
    seen: set[str] = set()
    for row in rows:
        sample = row["sample"]
        for col in ("sample", "group", "short_reads_1", "short_reads_2", "short_reads_platform"):
            if not row[col]:
                fail(f"Sample {sample or '<blank>'} is missing required field '{col}'")
        if sample in seen:
            fail(f"Duplicate sample id in {samplesheet_path}: {sample}")
        seen.add(sample)
        if row["short_reads_platform"] != "ILLUMINA":
            fail(f"Sample {sample} short_reads_platform must be ILLUMINA")
        if row["group"] != sample:
            fail(f"Sample {sample} must use a separate group matching the sample id")
    if check_files:
        for row in rows:
            for col in ("short_reads_1", "short_reads_2"):
                path = Path(row[col]).expanduser()
                require_path(
                    path,
                    f"Sample {row['sample']} {col}",
                    "file",
                )
    return rows


def read_contig_to_genus(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    if not path.exists():
        return mapping
    with path.open() as handle:
        for line_no, line in enumerate(handle, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if line_no == 1 and fields[:2] == ["contig_id", "genus"]:
                continue
            if len(fields) < 2:
                fail(f"Malformed contig-to-genus row in {path}: {line}")
            mapping[fields[0]] = fields[1]
    return mapping


def write_contig_to_genus(path: Path, mapping: dict[str, str]) -> None:
    ensure_dir(path.parent)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["contig_id", "genus"])
        for contig in sorted(mapping):
            writer.writerow([contig, mapping[contig]])


def shell_assignments(
    values: dict[str, str | int | Path | Sequence[str | int | Path]],
) -> str:
    assignments = []
    for key, value in values.items():
        if isinstance(value, (list, tuple)):
            items = " ".join(shlex.quote(str(item)) for item in value)
            assignments.append(f"{key}=({items})")
        else:
            assignments.append(f"{key}={shlex.quote(str(value))}")
    return "\n".join(assignments)
