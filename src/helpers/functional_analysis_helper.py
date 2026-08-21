"""Functional analysis helpers used by the notebook."""

from __future__ import annotations

from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
import textwrap

import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
import pysam


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
    mapping["ko_clean"] = mapping["ko_raw"].astype(str).str.replace("ko:", "", regex=False).str.strip()
    mapping[target_col] = mapping[target_col].astype(str).str.strip()
    return mapping[["ko_clean", target_col]].dropna()


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


def compute_gene_signal(reads_vec, gene_length, meta_cov=None, meta_cov_pseudocount=1.0):
    raw = np.asarray(reads_vec, dtype=float) / np.asarray(gene_length, dtype=float)
    if meta_cov is not None:
        raw = raw / (np.asarray(meta_cov, dtype=float) + meta_cov_pseudocount)
    raw = np.nan_to_num(raw, nan=0.0, posinf=0.0, neginf=0.0)
    scale = raw.sum()
    return (raw / scale) * 1_000_000 if scale > 0 else np.zeros_like(raw)


def bh_fdr(pvals):
    pvals = np.asarray(pvals, dtype=float)
    n = pvals.size
    order = np.argsort(pvals)
    ranked = pvals[order]
    q_ranked = ranked * n / (np.arange(1, n + 1))
    q_ranked = np.minimum.accumulate(q_ranked[::-1])[::-1]
    q_ranked = np.clip(q_ranked, 0, 1)
    qvals = np.empty_like(q_ranked)
    qvals[order] = q_ranked
    return qvals


def _permute_chunk(
    n_chunk,
    seed,
    total_reads,
    prob,
    gene_length,
    meta_cov,
    gene_indices,
    category_indices,
    n_categories,
):
    rng = np.random.default_rng(seed)
    out = np.zeros((n_chunk, n_categories), dtype=np.float64)
    for i in range(n_chunk):
        perm_reads = rng.multinomial(total_reads, prob)
        perm_norm = compute_gene_signal(perm_reads, gene_length, meta_cov=meta_cov)
        out[i, :] = np.bincount(
            category_indices,
            weights=perm_norm[gene_indices],
            minlength=n_categories,
        )
    return out


def run_permutations_parallel(
    n_perm,
    total_reads,
    prob,
    gene_length,
    meta_cov,
    gene_indices,
    category_indices,
    n_categories,
    seed=42,
    n_threads=1,
):
    n_threads = max(1, int(n_threads))
    n_perm = int(n_perm)
    if n_perm <= 0:
        return np.zeros((0, n_categories), dtype=np.float64)
    if n_threads == 1:
        return _permute_chunk(
            n_perm,
            seed,
            total_reads,
            prob,
            gene_length,
            meta_cov,
            gene_indices,
            category_indices,
            n_categories,
        )

    chunk_sizes = [n_perm // n_threads] * n_threads
    for i in range(n_perm % n_threads):
        chunk_sizes[i] += 1
    chunk_sizes = [c for c in chunk_sizes if c > 0]
    child_seeds = np.random.SeedSequence(seed).spawn(len(chunk_sizes))

    blocks = [None] * len(chunk_sizes)
    with ThreadPoolExecutor(max_workers=n_threads) as ex:
        futures = {
            ex.submit(
                _permute_chunk,
                chunk_n,
                int(child_seed.generate_state(1)[0]),
                total_reads,
                prob,
                gene_length,
                meta_cov,
                gene_indices,
                category_indices,
                n_categories,
            ): idx
            for idx, (chunk_n, child_seed) in enumerate(zip(chunk_sizes, child_seeds))
        }
        for fut in as_completed(futures):
            blocks[futures[fut]] = fut.result()
    return np.vstack(blocks)


def _category_assignments(df, mapping_path=None, target_col="pathway"):
    ko_long = _ko_long(df)
    if mapping_path is None:
        assignments = ko_long.rename(columns={"ko_clean": target_col})
    else:
        mapping = _load_mapping(mapping_path, target_col=target_col)
        assignments = ko_long.merge(mapping, on="ko_clean", how="inner")
    assignments[target_col] = assignments[target_col].astype(str)
    return assignments[["gene_idx", target_col]].drop_duplicates()


def null_functional_test(
    df,
    mapping_path=None,
    n_perm=10000,
    seed=42,
    n_threads=1,
    target_col="pathway",
    min_observed_reads=2,
):
    # Permute ribo counts across genes to see which functions are higher or lower than expected.
    work = df.copy()
    if "binned_counts" in work.columns:
        observed_reads = _numeric(work["binned_counts"]).to_numpy(dtype=np.int64)
    elif "num_footprints" in work.columns:
        observed_reads = _numeric(work["num_footprints"]).to_numpy(dtype=np.int64)
    else:
        observed_reads = np.zeros(len(work), dtype=np.int64)

    gene_length = _numeric(work["gene_length"]).replace(0, np.nan).fillna(1).to_numpy(dtype=np.float64)
    meta_cov = _numeric(work["coverage"]).to_numpy(dtype=np.float64) if "coverage" in work.columns else None
    assignments = _category_assignments(work, mapping_path=mapping_path, target_col=target_col)
    if assignments.empty:
        categories = pd.Index([], name=target_col)
        return {"stats": pd.DataFrame(), "perm": np.zeros((int(n_perm), 0)), "targets": categories}

    category_reads = (
        assignments.assign(
            observed_reads=observed_reads[assignments["gene_idx"].to_numpy()]
        )
        .groupby(target_col, sort=False)["observed_reads"]
        .sum()
    )
    category_reads = category_reads[category_reads >= int(min_observed_reads)]
    categories = pd.Index(category_reads.index, name=target_col)
    if categories.empty:
        return {"stats": pd.DataFrame(), "perm": np.zeros((int(n_perm), 0)), "targets": categories}

    assignments = assignments[assignments[target_col].isin(categories)]
    gene_indices = assignments["gene_idx"].to_numpy(dtype=np.int32)
    category_indices = categories.get_indexer(assignments[target_col]).astype(np.int32)
    observed_category_reads = category_reads.to_numpy(dtype=np.int64)

    obs_norm = compute_gene_signal(observed_reads, gene_length, meta_cov=meta_cov)
    observed_category = np.bincount(
        category_indices,
        weights=obs_norm[gene_indices],
        minlength=len(categories),
    )
    total_reads = int(observed_reads.sum())
    prob_weights = gene_length * (meta_cov + 1.0 if meta_cov is not None else 1.0)
    prob = prob_weights / prob_weights.sum() if prob_weights.sum() > 0 else np.repeat(1 / len(gene_length), len(gene_length))

    permuted_categories = run_permutations_parallel(
        n_perm=n_perm,
        total_reads=total_reads,
        prob=prob,
        gene_length=gene_length,
        meta_cov=meta_cov,
        gene_indices=gene_indices,
        category_indices=category_indices,
        n_categories=len(categories),
        seed=seed,
        n_threads=n_threads,
    )
    perm_mean = permuted_categories.mean(axis=0) if len(permuted_categories) else np.zeros(len(categories))
    perm_sd = permuted_categories.std(axis=0, ddof=1) if len(permuted_categories) > 1 else np.zeros(len(categories))
    denom = max(1, int(n_perm))
    p_enrich = (1 + (permuted_categories >= observed_category).sum(axis=0)) / (1 + denom)
    p_deplete = (1 + (permuted_categories <= observed_category).sum(axis=0)) / (1 + denom)
    abs_obs_dev = np.abs(observed_category - perm_mean)
    abs_perm_dev = np.abs(permuted_categories - perm_mean[None, :]) if len(permuted_categories) else np.zeros((0, len(categories)))
    p_two_sided = (1 + (abs_perm_dev >= abs_obs_dev[None, :]).sum(axis=0)) / (1 + denom)
    eps = 1e-12

    stats = pd.DataFrame(
        {
            target_col: categories,
            "target": categories,
            "target_display": categories,
            "observed_reads": observed_category_reads,
            "observed_norm": observed_category,
            "perm_mean_norm": perm_mean,
            "perm_sd_norm": perm_sd,
            "log2_fc_obs_vs_null": np.log2((observed_category + eps) / (perm_mean + eps)),
            "z_score": (observed_category - perm_mean) / (perm_sd + eps),
            "p_enrich": p_enrich,
            "p_deplete": p_deplete,
            "p_two_sided": p_two_sided,
            "q_enrich": bh_fdr(p_enrich),
            "q_deplete": bh_fdr(p_deplete),
            "q_two_sided": bh_fdr(p_two_sided),
        }
    ).sort_values("q_two_sided", ascending=True)
    stats["neglog10_q"] = -np.log10(np.clip(stats["q_two_sided"], 1e-300, 1))
    return {
        "stats": stats,
        "perm": permuted_categories,
        "targets": categories,
        "observed": observed_category,
        "target_col": target_col,
    }


def plot_functional_test(null_results, log2_thr=1.0, sig_thr=2.5, max_labels=20):
    stats = null_results["stats"] if isinstance(null_results, dict) else null_results
    fig, ax = plt.subplots()
    if stats.empty:
        ax.set_title("Pathway Permutation Significance")
        return fig
    vol = stats.copy()
    vol["direction"] = np.where(vol["log2_fc_obs_vs_null"] >= 0, "Enriched", "Depleted")
    for direction, color in [("Enriched", "tab:red"), ("Depleted", "tab:blue")]:
        sub = vol[vol["direction"] == direction]
        ax.scatter(sub["log2_fc_obs_vs_null"], sub["neglog10_q"], s=20, alpha=0.6, c=color, label=direction)
    ax.axhline(sig_thr, color="black", linestyle="--", linewidth=1, label=f"-log10(q) = {sig_thr}")
    ax.axvline(log2_thr, color="gray", linestyle=":", linewidth=1)
    ax.axvline(-log2_thr, color="gray", linestyle=":", linewidth=1)
    label_df = (
        vol[
            (vol["neglog10_q"] > sig_thr)
            & (vol["log2_fc_obs_vs_null"].abs() > log2_thr)
        ]
        .sort_values("neglog10_q", ascending=False)
        .head(max_labels)
    )
    for _, row in label_df.iterrows():
        x = row["log2_fc_obs_vs_null"]
        y = row["neglog10_q"]
        ax.text(
            x,
            y,
            str(row.get("target_display", row.get("target", "")))[:60],
            fontsize=8,
        )
    ax.set_xlabel("log2(Observed / Null mean)")
    ax.set_ylabel("-log10(FDR q-value, two-sided)")
    ax.set_title("Pathway Permutation Significance")
    ax.legend()
    fig.tight_layout()
    return fig


def run_mapping_permutation(df, mapping_path, target_col, label, n_perm=1000, seed=42, n_threads=1):
    result = null_functional_test(df, mapping_path=mapping_path, n_perm=n_perm, seed=seed, n_threads=n_threads, target_col=target_col)
    result["stats"]["mapping"] = label
    return result


def run_multi_mapping_permutations(df, mapping_specs, module_name_path=None, n_perm=1000, seed=42, n_threads=1):
    # mapping_specs entries are (label, mapping_path, target_col).
    module_name_map = {}
    if module_name_path is not None:
        module_names = pd.read_csv(module_name_path, sep="\t", header=None, names=["module_clean", "module_name"])
        module_names["module_clean"] = module_names["module_clean"].astype(str).str.replace("md:", "", regex=False).str.strip()
        module_name_map = dict(zip(module_names["module_clean"], module_names["module_name"]))

    mapping_workers = min(len(mapping_specs), max(1, int(n_threads))) or 1
    perm_threads = max(1, int(n_threads) // mapping_workers)

    def _run_one(idx, spec):
        label, path, target_col = spec
        result = run_mapping_permutation(df, path, target_col, label, n_perm=n_perm, seed=seed + idx * 1000, n_threads=perm_threads)
        stats = result["stats"]
        if label == "KO_MODULE" and module_name_map:
            stats["module_clean"] = stats["target"].astype(str).str.replace("md:", "", regex=False).str.strip()
            stats["module_name"] = stats["module_clean"].map(module_name_map)
            stats["target_display"] = stats.apply(
                lambda row: f"{row['target']} | {row['module_name']}" if pd.notna(row["module_name"]) else row["target"],
                axis=1,
            )
        return label, result

    if mapping_workers == 1:
        return dict(_run_one(i, spec) for i, spec in enumerate(mapping_specs))

    results = {}
    with ThreadPoolExecutor(max_workers=mapping_workers) as ex:
        futures = {ex.submit(_run_one, i, spec): spec[0] for i, spec in enumerate(mapping_specs)}
        for fut in as_completed(futures):
            label, result = fut.result()
            results[label] = result
    return results
