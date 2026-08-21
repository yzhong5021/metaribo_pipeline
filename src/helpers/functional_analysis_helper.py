"""Helpers for the functional analysis notebook."""

from __future__ import annotations

from collections import Counter
import textwrap

import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
import pysam
from scipy.stats import binomtest, false_discovery_control


def _numeric(series, fill=0.0):
    return pd.to_numeric(series, errors="coerce").fillna(fill).astype(float)


def _signal_column(df: pd.DataFrame) -> pd.Series:
    if "normalized_footprints" in df.columns:
        return _numeric(df["normalized_footprints"])
    if "normalized_ribo" in df.columns:
        return _numeric(df["normalized_ribo"])
    return normalize_ribo(df)


def _clean_ko(value) -> list[str]:
    if pd.isna(value) or str(value).strip() in {"", "-", "nan"}:
        return []
    return [ko.replace("ko:", "").strip() for ko in str(value).split(",") if ko.replace("ko:", "").strip()]


def _load_mapping(mapping_path, target_col="kofam") -> pd.DataFrame:
    mapping = pd.read_csv(mapping_path, sep="\t", header=None, names=["ko_raw", target_col])
    mapping = mapping.dropna(subset=["ko_raw", target_col]).copy()
    mapping["ko_clean"] = mapping["ko_raw"].astype(str).str.replace("ko:", "", regex=False).str.strip()
    mapping[target_col] = mapping[target_col].astype(str).str.strip()
    return mapping[(mapping["ko_clean"] != "") & (mapping[target_col] != "")][["ko_clean", target_col]]


def _ko_long(df: pd.DataFrame) -> pd.DataFrame:
    ids = df["CDS_id"].astype(str) if "CDS_id" in df.columns else pd.Index(df.index).astype(str)
    ko_df = pd.DataFrame({"CDS_id": ids, "kofam_ko": df["kofam_ko"].values})
    ko_df["gene_idx"] = np.arange(len(ko_df), dtype=np.int32)
    ko_df["ko_clean"] = ko_df["kofam_ko"].apply(_clean_ko)
    ko_df = ko_df.explode("ko_clean")
    return ko_df[ko_df["ko_clean"].notna() & (ko_df["ko_clean"] != "")].copy()


def normalize_ribo(df):
    # Normalize ribo counts by gene length, total ribo signal, and metagenome coverage.
    if "binned_counts" in df.columns:
        counts = _numeric(df["binned_counts"])
    elif "num_footprints" in df.columns:
        counts = _numeric(df["num_footprints"])
    else:
        counts = pd.Series(0.0, index=df.index)

    gene_length = _numeric(df["gene_length"]).replace(0, np.nan)
    total_counts = counts.sum()
    rpkm = counts / (gene_length / 1000.0)
    if total_counts > 0:
        rpkm = rpkm / (total_counts / 1_000_000.0)

    coverage = _numeric(df["coverage"]).replace(0, np.nan)
    normalized = (rpkm / coverage).replace([np.inf, -np.inf], np.nan).fillna(0.0)
    normalized.name = "normalized_ribo"
    return normalized


def _taxon_props(df: pd.DataFrame, value_col: str, cutoff: float) -> pd.Series:
    work = df.copy()
    work["genus"] = work.get("genus", "Unknown")
    work["genus"] = work["genus"].replace("", np.nan).fillna("Unknown")
    values = _numeric(work[value_col]) if value_col in work.columns else _signal_column(work)
    grouped = values.groupby(work["genus"]).sum().sort_values(ascending=False)
    total = grouped.sum()
    props = grouped / total if total > 0 else grouped
    if props.empty:
        return pd.Series({"Other": 0.0})
    keep = props[props >= cutoff]
    other = props[props < cutoff].sum()
    if other > 0:
        keep.loc["Other"] = other
    return keep


def _plot_taxon_stack(props: pd.Series, title: str):
    fig, ax = plt.subplots(figsize=(8, 8))
    bottom = 0.0
    cmap = plt.colormaps.get_cmap("tab20")
    labels = [f"{idx} ({val:.3f})" for idx, val in props.items()]
    for i, (_idx, val) in enumerate(props.items()):
        color = "black" if _idx == "Other" else cmap(i % 20)
        ax.bar(0, val, bottom=bottom, color=color, width=0.5)
        bottom += val
    handles = [plt.Rectangle((0, 0), 1, 1, color=("black" if idx == "Other" else cmap(i % 20))) for i, idx in enumerate(props.index)]
    ax.set_title(title)
    ax.set_xticks([])
    ax.set_ylabel("Proportion")
    ax.set_ylim(0, 1)
    ax.legend(handles, labels, title="Genus", bbox_to_anchor=(1.05, 1), loc="upper left")
    fig.tight_layout()
    return fig


def plot_taxonomy(df, cutoff):
    # Compare genus distribution by metagenome coverage and normalized ribo signal.
    work = df.copy()
    if "normalized_ribo" not in work.columns:
        work["normalized_ribo"] = normalize_ribo(work)
    meta_props = _taxon_props(work, "coverage", cutoff)
    ribo_props = _taxon_props(work, "normalized_ribo", cutoff)
    return (
        _plot_taxon_stack(meta_props, "Metagenome Coverage Genus Distribution"),
        _plot_taxon_stack(ribo_props, "Ribo-seq Genus Distribution"),
    )


def taxonomy_read_length(df, bam_path, top_n=10):
    # Count ribo read lengths by genus using the contig-to-genus calls in the master sheet.
    try:
        import seaborn as sns
    except ImportError:
        sns = None

    contig_genus = (
        df[["contig_id", "genus"]]
        .dropna(subset=["contig_id"])
        .assign(genus=lambda x: x["genus"].replace("", np.nan).fillna("Unknown"))
        .drop_duplicates("contig_id")
        .set_index("contig_id")["genus"]
        .to_dict()
    )
    counts = Counter()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for aln in bam.fetch(until_eof=True):
            if aln.is_unmapped or aln.reference_id < 0:
                continue
            contig = bam.get_reference_name(aln.reference_id)
            genus = contig_genus.get(contig)
            if genus is None:
                continue
            read_len = aln.infer_query_length(always=True)
            if read_len is not None:
                counts[(read_len, genus)] += 1

    rows = [{"length": length, "genus": genus, "count": count} for (length, genus), count in counts.items()]
    if rows:
        readlens = pd.DataFrame(rows).pivot_table(index="length", columns="genus", values="count", fill_value=0).sort_index()
    else:
        readlens = pd.DataFrame(index=pd.Index([], name="length"))

    top_genera = readlens.sum().sort_values(ascending=False).head(top_n).index
    plot_data = readlens[top_genera].T if len(top_genera) else readlens.T
    plot_norm = plot_data.div(plot_data.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)

    fig, ax = plt.subplots(figsize=(10, max(3, 0.45 * max(1, len(plot_norm)))))
    if sns is not None:
        sns.heatmap(plot_norm, cmap="viridis", norm=PowerNorm(gamma=0.5), ax=ax, cbar_kws={"label": "Row-normalized proportion"})
    else:
        image = ax.imshow(plot_norm.to_numpy(), aspect="auto", cmap="viridis", norm=PowerNorm(gamma=0.5))
        ax.set_xticks(range(len(plot_norm.columns)))
        ax.set_xticklabels(plot_norm.columns)
        ax.set_yticks(range(len(plot_norm.index)))
        ax.set_yticklabels(plot_norm.index)
        fig.colorbar(image, ax=ax, label="Row-normalized proportion")
    for row_idx, row_label in enumerate(plot_norm.index):
        if plot_norm.shape[1] == 0:
            continue
        max_col = plot_norm.loc[row_label].idxmax()
        x = plot_norm.columns.get_loc(max_col)
        ax.add_patch(Rectangle((x, row_idx), 1, 1, fill=False, edgecolor="white", linewidth=0.5))
    ax.set_title(f"Read Length Distribution, Top {top_n} Genera")
    ax.set_xlabel("Read length")
    ax.set_ylabel("Genus")
    fig.tight_layout()
    return fig


def process_footprints(df, mapping_path):
    # Map KO-level ribo signal to functional labels before plotting pathways.
    mapping = _load_mapping(mapping_path, target_col="kofam").set_index("ko_clean")["kofam"]
    work = df.copy()
    signal = _signal_column(work).reset_index(drop=True)
    genus = work.get("genus", pd.Series("Unknown", index=work.index)).reset_index(drop=True)
    genus = genus.replace("", np.nan).fillna("Unknown")

    ko_long = _ko_long(work)
    ko_long["normalized_footprints"] = signal.iloc[ko_long["gene_idx"].to_numpy()].to_numpy()
    ko_long["genus"] = genus.iloc[ko_long["gene_idx"].to_numpy()].to_numpy()
    exploded = ko_long
    exploded = exploded.join(mapping, on="ko_clean", how="left")
    exploded = exploded[(exploded["kofam"].notna()) & (exploded["kofam"] != "Function unknown")]
    exploded = exploded.drop_duplicates(["gene_idx", "kofam"])

    function_totals = exploded.groupby("kofam", dropna=False)["normalized_footprints"].sum().reset_index()
    function_totals.sort_values(by="normalized_footprints", ascending=False, inplace=True)

    function_by_genus = exploded.groupby(["kofam", "genus"], dropna=False)["normalized_footprints"].sum().unstack(fill_value=0).reset_index()
    function_by_genus = function_by_genus.assign(total=function_by_genus.drop(columns="kofam").sum(axis=1))
    function_by_genus.sort_values(by="total", ascending=False, inplace=True)
    function_by_genus.drop(columns="total", inplace=True)
    return function_totals, function_by_genus


def _format_stacked_barplot(ax, plot_data, title, xlabel, fontsize):
    ax.set_title(title, fontsize=18, pad=20)
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel("")
    pathways = [textwrap.fill(str(name), width=60, max_lines=2, placeholder="...") for name in plot_data.index]
    ax.set_yticklabels(pathways, fontsize=fontsize, va="center")
    ax.tick_params(axis="x", which="major", labelsize=8, rotation=30)
    ax.tick_params(axis="y", which="major", labelsize=fontsize, length=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.legend(bbox_to_anchor=(0.5, -0.1), ncols=4, loc="upper center", fontsize=14)


def plot_functional(df, mapping_path, top_n, top_genera, title, fontsize):
    _function_totals, function_by_genus = process_footprints(df, mapping_path)
    if function_by_genus.empty:
        fig, ax = plt.subplots(figsize=(12, 10))
        ax.set_title(title)
        return fig

    aggregated_df = function_by_genus.set_index("kofam")
    aggregated_df["total_score"] = aggregated_df.sum(axis=1)
    top_pathways_df = aggregated_df.sort_values(by="total_score", ascending=False).head(top_n).drop(columns=["total_score"])
    genera = top_pathways_df.sum().sort_values(ascending=False).head(top_genera).index
    plot_data = top_pathways_df[genera].iloc[::-1]

    fig, ax = plt.subplots(figsize=(12, 10))
    plot_data.plot(kind="barh", stacked=True, ax=ax, colormap="tab20", width=0.9)
    _format_stacked_barplot(ax, plot_data, str(title), "Relative Translational Activity", fontsize)
    fig.tight_layout()
    return fig


def _kegg_assignments(df, mapping_path):
    ko_long = _ko_long(df)
    mapping = _load_mapping(mapping_path, target_col="kegg_category")
    assignments = ko_long.merge(mapping, on="ko_clean", how="inner")
    return assignments[["gene_idx", "kegg_category"]].drop_duplicates()


def test_kegg_categories(df, mapping_path):
    if "binned_counts" in df.columns:
        observed_reads = _numeric(df["binned_counts"]).to_numpy(dtype=np.int64)
    elif "num_footprints" in df.columns:
        observed_reads = _numeric(df["num_footprints"]).to_numpy(dtype=np.int64)
    else:
        observed_reads = np.zeros(len(df), dtype=np.int64)

    gene_length = _numeric(df["gene_length"]).to_numpy(dtype=np.float64)
    coverage = _numeric(df["coverage"]).to_numpy(dtype=np.float64)
    genomic_opportunity = gene_length * coverage
    assignments = _kegg_assignments(df, mapping_path)
    eligible = genomic_opportunity > 0
    assignments = assignments[eligible[assignments["gene_idx"].to_numpy()]]
    if assignments.empty:
        return pd.DataFrame()

    total_reads = int(observed_reads[eligible].sum())
    if total_reads == 0:
        raise ValueError("No Ribo-seq reads are assigned to CDSs")

    gene_indices = assignments["gene_idx"].to_numpy(dtype=np.int32)
    stats = (
        assignments.assign(
            observed_reads=observed_reads[gene_indices],
            genomic_opportunity=genomic_opportunity[gene_indices],
        )
        .groupby("kegg_category", sort=False, as_index=False)
        .agg(
            observed_reads=("observed_reads", "sum"),
            genomic_opportunity=("genomic_opportunity", "sum"),
        )
    )
    stats["expected_probability"] = (
        stats["genomic_opportunity"] / genomic_opportunity[eligible].sum()
    ).clip(0.0, 1.0)
    stats["expected_reads"] = total_reads * stats["expected_probability"]
    with np.errstate(divide="ignore"):
        stats["log2_effect"] = np.log2(
            (stats["observed_reads"] / total_reads) / stats["expected_probability"]
        )
    stats["p_two_sided"] = [
        binomtest(int(reads), total_reads, float(probability), alternative="two-sided").pvalue
        for reads, probability in zip(stats["observed_reads"], stats["expected_probability"])
    ]
    stats["q_two_sided"] = false_discovery_control(stats["p_two_sided"])
    return stats.sort_values("q_two_sided", ascending=True)


def plot_kegg_test(stats, log2_thr=1.0, sig_thr=2.5, max_labels=20):
    fig, ax = plt.subplots()
    if stats.empty:
        ax.set_title("KEGG Category Enrichment")
        return fig
    vol = stats.copy()
    vol["neglog10_q"] = -np.log10(np.clip(vol["q_two_sided"], 1e-300, 1))
    vol["direction"] = np.where(vol["log2_effect"] >= 0, "Enriched", "Depleted")
    finite_effects = vol.loc[np.isfinite(vol["log2_effect"]), "log2_effect"]
    left_edge = finite_effects.min() - 1 if not finite_effects.empty else -log2_thr - 1
    vol["plot_effect"] = vol["log2_effect"].replace(-np.inf, left_edge)
    for direction, color in [("Enriched", "tab:red"), ("Depleted", "tab:blue")]:
        sub = vol[vol["direction"] == direction]
        ax.scatter(sub["plot_effect"], sub["neglog10_q"], s=20, alpha=0.6, c=color, label=direction)
    ax.axhline(sig_thr, color="black", linestyle="--", linewidth=1, label=f"-log10(q) = {sig_thr}")
    ax.axvline(log2_thr, color="gray", linestyle=":", linewidth=1)
    ax.axvline(-log2_thr, color="gray", linestyle=":", linewidth=1)
    label_df = (
        vol[
            (vol["neglog10_q"] > sig_thr)
            & (vol["log2_effect"].abs() > log2_thr)
        ]
        .sort_values("neglog10_q", ascending=False)
        .head(max_labels)
    )
    for _, row in label_df.iterrows():
        x = row["plot_effect"]
        y = row["neglog10_q"]
        ax.text(
            x,
            y,
            str(row["kegg_category"])[:60],
            fontsize=8,
        )
    ax.set_xlabel("log2(Observed / Expected)")
    ax.set_ylabel("-log10(q-value)")
    ax.set_title("KEGG Category Enrichment")
    ax.legend()
    fig.tight_layout()
    return fig
