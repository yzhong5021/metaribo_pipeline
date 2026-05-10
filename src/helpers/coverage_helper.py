"""Coverage and metagene helpers used by the coverage notebook."""

from __future__ import annotations

from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam


POINT_METHODS = {
    "5prime": "5p",
    "5p": "5p",
    "3prime": "3p",
    "3p": "3p",
    "median": "center",
    "center": "center",
}


def _cds_id_and_rows(master_df: pd.DataFrame) -> tuple[pd.Index, pd.DataFrame]:
    df = master_df.copy()
    if "CDS_id" in df.columns:
        cds_ids = df["CDS_id"].astype(str)
        df = df.set_index(cds_ids, drop=False)
    else:
        cds_ids = df.index.astype(str)
        df = df.copy()
        df["CDS_id"] = cds_ids
    return pd.Index(cds_ids, name="CDS_id"), df


def _master_records(master_df: pd.DataFrame) -> dict[str, list[dict]]:
    _cds_ids, df = _cds_id_and_rows(master_df)
    required = {"CDS_id", "contig_id", "start", "end", "strand"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"master_df is missing required columns: {', '.join(sorted(missing))}")

    records: dict[str, list[dict]] = defaultdict(list)
    for _, row in df.iterrows():
        start_1 = int(row["start"])
        end_1 = int(row["end"])
        if end_1 < start_1:
            continue
        record = {
            "CDS_id": str(row["CDS_id"]),
            "contig_id": str(row["contig_id"]),
            "start0": start_1 - 1,
            "end0": end_1,
            "strand": str(row["strand"]),
        }
        records[record["contig_id"]].append(record)
    for contig in records:
        records[contig].sort(key=lambda r: (r["start0"], r["end0"], r["CDS_id"]))
    return records


def _valid_alignment(aln: pysam.AlignedSegment) -> bool:
    return not (
        aln.is_unmapped
        or aln.is_secondary
        or aln.is_supplementary
        or aln.is_duplicate
        or aln.is_qcfail
        or aln.reference_start is None
        or aln.reference_end is None
    )


def _read_length(aln: pysam.AlignedSegment) -> int | None:
    qlen = aln.query_length
    return qlen if qlen is not None else aln.infer_query_length(always=False)


def _same_strand(aln: pysam.AlignedSegment, cds_strand: str) -> bool:
    if cds_strand not in {"+", "-"}:
        return True
    return aln.is_reverse == (cds_strand == "-")


def _read_point(aln: pysam.AlignedSegment, method: str) -> int:
    method = POINT_METHODS[method]
    if method == "center":
        return (aln.reference_start + aln.reference_end - 1) // 2
    if method == "5p":
        return aln.reference_end - 1 if aln.is_reverse else aln.reference_start
    return aln.reference_start if aln.is_reverse else aln.reference_end - 1


def _overlap_len(aln: pysam.AlignedSegment, record: dict) -> int:
    return max(0, min(aln.reference_end, record["end0"]) - max(aln.reference_start, record["start0"]))


def bin_ribo_counts(bam_path, master_df, method="overlap"):
    # Bin each usable ribo read to one CDS; overlap is the default assignment method.
    method = method.lower()
    if method not in {"overlap", *POINT_METHODS}:
        raise ValueError("method must be one of: overlap, median, center, 5prime, 5p, 3prime, 3p")

    cds_ids, _df = _cds_id_and_rows(master_df)
    counts = pd.Series(0, index=cds_ids, name="binned_counts", dtype="int64")
    records_by_contig = _master_records(master_df)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for aln in bam.fetch(until_eof=True):
            if not _valid_alignment(aln) or aln.reference_id < 0:
                continue
            contig = bam.get_reference_name(aln.reference_id)
            records = records_by_contig.get(contig, [])
            if not records:
                continue

            best_cds = None
            best_score = 0
            if method == "overlap":
                for record in records:
                    if not _same_strand(aln, record["strand"]):
                        continue
                    score = _overlap_len(aln, record)
                    if score > best_score:
                        best_cds = record["CDS_id"]
                        best_score = score
            else:
                point = _read_point(aln, method)
                for record in records:
                    if not _same_strand(aln, record["strand"]):
                        continue
                    if record["start0"] <= point < record["end0"]:
                        best_cds = record["CDS_id"]
                        break

            if best_cds is not None:
                counts.loc[best_cds] += 1

    return counts


def build_metagene(bam_path, master_df, anchor, window, length=None):
    # Build start and stop metagene profiles around Bakta CDS boundaries.
    anchor = POINT_METHODS.get(str(anchor).lower())
    if anchor is None:
        raise ValueError("anchor must be one of: median, center, 5prime, 5p, 3prime, 3p")
    window = int(window)
    length = None if length is None else int(length)
    L = window * 2

    metagene_start = np.zeros(L, dtype=float)
    metagene_stop = np.zeros(L, dtype=float)
    n_used_start = 0
    n_used_stop = 0
    records_by_contig = _master_records(master_df)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        refs = set(bam.references)
        for contig, records in records_by_contig.items():
            if contig not in refs:
                continue
            contig_len = bam.get_reference_length(contig)
            for record in records:
                start0 = record["start0"]
                end0 = record["end0"]
                if end0 <= start0:
                    continue

                is_minus = record["strand"] == "-"
                anchor_start = end0 - 1 if is_minus else start0
                anchor_stop = start0 + 2 if is_minus else end0 - 3
                edge_pad = window + ((length - 1) if length is not None else 0)

                has_start = edge_pad <= anchor_start < contig_len - edge_pad
                has_stop = edge_pad <= anchor_stop < contig_len - edge_pad
                if not has_start and not has_stop:
                    continue

                v_start = np.zeros(L, dtype=float)
                v_stop = np.zeros(L, dtype=float)
                fetch_start = max(0, min(anchor_start, anchor_stop) - edge_pad)
                fetch_end = min(contig_len, max(anchor_start, anchor_stop) + edge_pad + 1)

                for aln in bam.fetch(contig, fetch_start, fetch_end):
                    if not _valid_alignment(aln) or not _same_strand(aln, record["strand"]):
                        continue
                    qlen = _read_length(aln)
                    if length is not None and qlen != length:
                        continue

                    p = _read_point(aln, anchor)
                    if has_start:
                        rel = p - anchor_start
                        if is_minus:
                            rel = -rel
                        if -window <= rel < window:
                            v_start[rel + window] += 1

                    if has_stop:
                        rel = p - anchor_stop
                        if is_minus:
                            rel = -rel
                        if -window <= rel < window:
                            v_stop[rel + window] += 1

                if has_start:
                    metagene_start += v_start
                    n_used_start += 1
                if has_stop:
                    metagene_stop += v_stop
                    n_used_stop += 1

    positions = np.arange(-window, window)
    return pd.DataFrame(
        {
            "position": positions,
            "metagene_start": metagene_start,
            "metagene_stop": metagene_stop,
            "metagene_start_mean": metagene_start / n_used_start if n_used_start else metagene_start,
            "metagene_stop_mean": metagene_stop / n_used_stop if n_used_stop else metagene_stop,
            "n_used_start": n_used_start,
            "n_used_stop": n_used_stop,
        }
    )


def plot_metagene(metagene_df, anchor, window, title):
    """
    Plots the start and stop metagene means within a specified window.
    """
    if isinstance(metagene_df, dict):
        x = np.arange(-window, window)
        start = np.asarray(metagene_df["metagene_start_mean"])
        stop = np.asarray(metagene_df["metagene_stop_mean"])
    else:
        x = metagene_df["position"].to_numpy()
        start = metagene_df["metagene_start_mean"].to_numpy()
        stop = metagene_df["metagene_stop_mean"].to_numpy()

    fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharey=True)
    for ax, y, label, color in [
        (axes[0], start, "Start", "tab:blue"),
        (axes[1], stop, "Stop", "tab:red"),
    ]:
        ax.plot(x, y, color=color)
        ax.axvline(0, color="black", linewidth=0.8, linestyle="--")
        ax.set_title(f"{label}")
        ax.set_xlabel("Distance (nt)")
        ax.set_xlim(-window, window - 1)
    axes[0].set_ylabel("Mean ribo-seq counts per CDS")
    fig.suptitle(str(title).strip())
    fig.tight_layout()
    return fig
