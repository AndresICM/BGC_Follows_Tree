#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import numpy as np
import pandas as pd
from scipy import sparse

GENOME_TABLE = "20260610_Taxonomy_HighQActinos.csv"
GENOME_COL   = "Genome"

# --- Distance input options 
NPY_DIST     = "pcoa_out/jaccard_distance.npy"
NPY_ORDER    = "pcoa_out/genome_order.txt"
# square CSV/TSV distance matrix exported from R (rownames=genomes)
CSV_DIST     = "20260611_jac_matrix.csv"          # set to "" to skip
#rebuild from the BGC table if no distance file is found
BGC_TABLE    = "bgc_table.tsv"
GCF_COL      = "gcf_id"
MIN_GENOMES_PER_GCF = 1                   # set 2 to drop singleton GCFs


FACTORS      = ["Class", "Order", "Family", "Genus", "Clade"]  # ranks to test
PERMUTATIONS = 999
MIN_GROUP_N  = 25            # drop groups with fewer than this many genomes
RANDOM_SEED  = 42
MTC_METHOD   = "holm" #choose holm or bonferroni
RESULTS_TSV  = "permanova_results_python.tsv"
RESULTS_TXT  = "permanova_results_python.txt"

# Genus-balanced analyses retain genera with at least this many genomes,
# then randomly sample exactly this number from each retained genus.
GENOMES_PER_GENUS = 25

# Independent balanced samples.
# Start with 3 to check runtime; increase to 20 or more for the final analysis.
N_BALANCED_REPLICATES = 20

EXCLUDE_GENERA = ["Streptomyces", "Rhodococcus"]

# Minimum number of genomes per tested taxonomic group AFTER sampling.
# Use the same threshold for all scenarios.
# Must not exceed GENOMES_PER_GENUS when testing Genus.
MIN_GROUP_N = 1

RESULTS_TSV = "permanova_sampling_sensitivity.tsv"
RESULTS_TXT = "permanova_sampling_sensitivity.txt"
SAMPLED_GENOMES_TSV = "permanova_sampling_genomes.tsv"

def log(m): print(f"[PERMANOVA] {m}", flush=True)


def load_distance():
    # A1: numpy
    if os.path.exists(NPY_DIST) and os.path.exists(NPY_ORDER):
        log(f"Loading numpy distance: {NPY_DIST}")
        dist = np.load(NPY_DIST).astype(np.float32)
        with open(NPY_ORDER) as fh:
            ids = [l.strip() for l in fh if l.strip()]
        return dist, ids
    # A2: CSV/TSV square matrix from R
    if CSV_DIST and os.path.exists(CSV_DIST):
        log(f"Loading square distance matrix: {CSV_DIST}")
        sep = "\t" if CSV_DIST.endswith((".tsv", ".txt")) else ","
        dm = pd.read_csv(CSV_DIST, sep=sep, index_col=0)
        ids = [str(x) for x in dm.index]
        return dm.values.astype(np.float32), ids
    # B: rebuild from BGC table
    log("No distance file found -> rebuilding Jaccard from BGC table (float32).")
    return rebuild_distance()


def rebuild_distance():
    df = pd.read_csv(BGC_TABLE, sep="\t", usecols=[GENOME_COL, GCF_COL], dtype=str)
    df = df.dropna(subset=[GENOME_COL, GCF_COL]).drop_duplicates([GENOME_COL, GCF_COL])
    gcodes, gidx = pd.factorize(df[GENOME_COL], sort=True)
    fcodes, fidx = pd.factorize(df[GCF_COL], sort=True)
    mat = sparse.csr_matrix((np.ones(len(df), np.float32), (gcodes, fcodes)),
                            shape=(len(gidx), len(fidx)), dtype=np.float32)
    mat.data[:] = 1
    if MIN_GENOMES_PER_GCF > 1:
        keep = np.asarray((mat > 0).sum(axis=0)).ravel() >= MIN_GENOMES_PER_GCF
        mat = mat[:, keep]
    nonempty = np.asarray((mat > 0).sum(axis=1)).ravel() > 0
    mat, gidx = mat[nonempty, :], gidx[nonempty]

    inter = (mat @ mat.T).toarray().astype(np.float32)
    a = np.asarray(mat.sum(axis=1)).ravel().astype(np.float32)
    union = a[:, None] + a[None, :] - inter
    with np.errstate(divide="ignore", invalid="ignore"):
        sim = np.where(union > 0, inter / union, 0.0).astype(np.float32)
    dist = 1.0 - sim
    np.fill_diagonal(dist, 0.0)
    dist = 0.5 * (dist + dist.T)
    np.clip(dist, 0.0, 1.0, out=dist)
    # cache so you never rebuild again
    os.makedirs("pcoa_out", exist_ok=True)
    np.save(NPY_DIST, dist)
    with open(NPY_ORDER, "w") as fh:
        fh.write("\n".join(str(g) for g in gidx))
    log("Cached rebuilt distance to pcoa_out/ for future runs.")
    return dist, [str(g) for g in gidx]



def run_one(dist, ids, meta, factor, permutations, min_group_n, seed
):
    from skbio.stats.distance import DistanceMatrix, permanova, permdisp

    if factor not in meta.columns:
        log(f"  SKIP '{factor}': not in metadata columns."); return None

    grouping = meta.set_index(GENOME_COL).reindex(ids)[factor]
    valid = grouping.notna() & (grouping.astype(str).str.strip() != "")
    counts = grouping[valid].value_counts()
    big = counts[counts >= min_group_n].index
    valid = valid & grouping.isin(big)

    idx = np.where(valid.values)[0]
    if len(idx) < 3 or grouping[valid].nunique() < 2:
        log(f"  SKIP '{factor}': too few samples/groups after filtering."); return None

    sub_ids = [ids[i] for i in idx]
    sub_dist = dist[np.ix_(idx, idx)]
    sub_grp = grouping.iloc[idx].astype(str).values
    log(f"  '{factor}': {len(sub_ids):,} genomes, {len(set(sub_grp)):,} groups (>= {min_group_n}).")

    dm = DistanceMatrix(sub_dist, ids=sub_ids)
    pmv = permanova(dm,
    grouping=list(sub_grp),
    permutations=permutations,
    seed=seed,)
    try:
        pdisp = permdisp(dm, grouping=list(sub_grp), permutations=permutations, seed=seed)
        dF, dp = pdisp["test statistic"], pdisp["p-value"]
    except Exception as e:
        dF, dp = np.nan, np.nan
        log(f"  PERMDISP failed for '{factor}': {e}")

    n, k = len(sub_ids), len(set(sub_grp))
    F = pmv["test statistic"]
    r2 = (F * (k - 1)) / (F * (k - 1) + (n - k)) if np.isfinite(F) else np.nan
    log(f"    PERMANOVA pseudo-F={F:.3f}  R2={r2:.3f}  p={pmv['p-value']}")
    return {"factor": factor, "n_genomes": n, "n_groups": k,
            "permanova_pseudoF": F, "permanova_R2": r2, "permanova_p": pmv["p-value"],
            "permdisp_F": dF, "permdisp_p": dp, "permutations": permutations}


def main():
    # Sampling settings
    MAX_GENOMES_PER_GENUS = 25
    N_REPLICATES = 20

    factors = FACTORS

    # Load the existing distance matrix.
    dist, ids = load_distance()

    if dist.shape != (len(ids), len(ids)):
        raise ValueError(
            "Distance matrix dimensions do not match genome IDs."
        )

    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate genome IDs in the distance matrix.")

    log(f"Distance matrix: {dist.shape}, {len(ids):,} genomes")

    # Load and clean metadata.
    meta = pd.read_csv(GENOME_TABLE, sep=",", dtype=str)

    required = [GENOME_COL, "Genus"]
    missing = [column for column in required if column not in meta.columns]

    if missing:
        raise ValueError(f"Missing metadata columns: {missing}")

    for column in set([GENOME_COL, "Genus"] + factors):
        if column in meta.columns:
            meta[column] = meta[column].astype("string").str.strip()
            meta[column] = meta[column].replace("", pd.NA)

    if meta[GENOME_COL].isna().any():
        raise ValueError("Metadata contains missing genome IDs.")

    if meta[GENOME_COL].duplicated().any():
        raise ValueError("Metadata contains duplicate genome IDs.")

    # Retain matched genomes with known genus assignments.
    pool = (
        meta.loc[
            meta[GENOME_COL].isin(ids) & meta["Genus"].notna()
        ]
        .sort_values(GENOME_COL)
        .reset_index(drop=True)
    )

    if pool.empty:
        raise ValueError(
            "No matched genomes with genus assignments were found."
        )

    genus_counts = pool.groupby("Genus")[GENOME_COL].nunique()

    expected_genomes = int(
        genus_counts.clip(upper=MAX_GENOMES_PER_GENUS).sum()
    )

    log(
        f"Matched input: {len(pool):,} genomes; "
        f"{len(genus_counts):,} genera"
    )
    log(
        f"Genus-capped sampling: maximum "
        f"{MAX_GENOMES_PER_GENUS} genomes per genus"
    )
    log(
        f"Each replicate: {expected_genomes:,} genomes; "
        f"{len(genus_counts):,} genera"
    )
    log(
        f"Genera above the cap: "
        f"{int((genus_counts > MAX_GENOMES_PER_GENUS).sum()):,}"
    )
    log(f"Minimum tested group size: {MIN_GROUP_N}")

    rows = []
    sampled_records = []

    for replicate in range(1, N_REPLICATES + 1):
        seed = RANDOM_SEED + replicate
        rng = np.random.default_rng(seed)

        selected_parts = []

        for genus, group in pool.groupby("Genus", sort=True):
            if len(group) <= MAX_GENOMES_PER_GENUS:
                selected_parts.append(group)
            else:
                selected_parts.append(
                    group.sample(
                        n=MAX_GENOMES_PER_GENUS,
                        replace=False,
                        random_state=int(
                            rng.integers(0, 2**31 - 1)
                        ),
                    )
                )

        selected = pd.concat(selected_parts, ignore_index=True)

        log(
            f"\nReplicate {replicate}/{N_REPLICATES}: "
            f"{len(selected):,} genomes; "
            f"{selected['Genus'].nunique():,} genera"
        )

        # Preserve the original distance-matrix ordering.
        selected_set = set(selected[GENOME_COL])

        idx = np.array(
            [
                i for i, genome in enumerate(ids)
                if genome in selected_set
            ],
            dtype=int,
        )

        sub_ids = [ids[i] for i in idx]
        sub_dist = dist[np.ix_(idx, idx)]

        membership = selected[[GENOME_COL, "Genus"]].copy()
        membership["replicate"] = replicate
        sampled_records.append(membership)

        replicate_rows = []

        for factor in factors:
            result = run_one(
                sub_dist,
                sub_ids,
                selected,
                factor,
                PERMUTATIONS,
                MIN_GROUP_N,
                seed=seed,
            )

            if result is None:
                continue

            result.update({
                "scenario": "genus_capped",
                "replicate": replicate,
                "sampling_pool_genomes": len(selected),
                "sampling_pool_genera": selected["Genus"].nunique(),
                "maximum_genomes_per_genus": MAX_GENOMES_PER_GENUS,
                "minimum_tested_group_size": MIN_GROUP_N,
            })

            replicate_rows.append(result)

        if replicate_rows:
            block = pd.DataFrame(replicate_rows)

            # Correct across tested ranks within this replicate.
            if MTC_METHOD:
                from statsmodels.stats.multitest import multipletests

                method = MTC_METHOD.lower()

                if method not in ["holm", "bonferroni"]:
                    raise ValueError(
                        f"Unknown MTC_METHOD: {MTC_METHOD}"
                    )

                block["permanova_p_adj"] = multipletests(
                    block["permanova_p"].to_numpy(),
                    method=method,
                )[1]

                block["p_adj_method"] = method

            rows.extend(block.to_dict("records"))

            # Save completed results after every replicate.
            pd.DataFrame(rows).to_csv(
                RESULTS_TSV, sep="\t", index=False
            )

        pd.concat(sampled_records, ignore_index=True).to_csv(
            SAMPLED_GENOMES_TSV, sep="\t", index=False
        )

        del sub_dist

    if not rows:
        log("No valid tests produced.")
        return

    out = pd.DataFrame(rows)

    # Summarize effect sizes across sampling replicates.
    summary = (
        out.groupby("factor", sort=False)
        .agg(
            N_replicates=("replicate", "nunique"),
            N_genomes_min=("n_genomes", "min"),
            N_genomes_max=("n_genomes", "max"),
            N_groups_min=("n_groups", "min"),
            N_groups_max=("n_groups", "max"),
            R2_median=("permanova_R2", "median"),
            R2_min=("permanova_R2", "min"),
            R2_max=("permanova_R2", "max"),
            PERMANOVA_P_max=("permanova_p", "max"),
            PERMDISP_P_min=("permdisp_p", "min"),
            PERMDISP_successful_tests=("permdisp_p", "count"),
        )
        .reset_index()
    )

    if "permanova_p_adj" in out.columns:
        adjusted = (
            out.groupby("factor")["permanova_p_adj"]
            .max()
            .rename("PERMANOVA_P_adjusted_max")
            .reset_index()
        )

        summary = summary.merge(adjusted, on="factor", how="left")

    out.to_csv(RESULTS_TSV, sep="\t", index=False)

    with open(RESULTS_TXT, "w") as fh:
        fh.write(
            "PERMANOVA of GCF composition: genus-capped sampling\n"
        )
        fh.write(
            f"Maximum genomes per genus: {MAX_GENOMES_PER_GENUS}\n"
            f"Sampling replicates: {N_REPLICATES}\n"
            f"Minimum tested group size: {MIN_GROUP_N}\n\n"
        )

        fh.write("SUMMARY ACROSS REPLICATES\n")
        fh.write(summary.to_string(index=False))

        fh.write("\n\nALL RESULTS\n")
        fh.write(out.to_string(index=False))

        fh.write(
            "\n\nNotes:\n"
            "- Genera below the sampling cap retain all genomes.\n"
            "- Rank-specific filtering can reduce the tested sample.\n"
            "- Failed PERMDISP tests are recorded as missing values.\n"
            "- Sampling replicates are not independent biological replicates.\n"
        )

    log(f"Saved detailed results -> {RESULTS_TSV}")
    log(f"Saved summary and details -> {RESULTS_TXT}")
    log(f"Saved sampled genome IDs -> {SAMPLED_GENOMES_TSV}")

    print("\nSUMMARY ACROSS REPLICATES")
    print(summary.to_string(index=False))

if __name__ == "__main__":
    main()
