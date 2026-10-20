#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import numpy as np
import pandas as pd
from scipy import sparse

# ------------------------------------------------------------------ CONFIG ----
# Metadata table (the ~30k-genome table with taxonomy). REQUIRED.
GENOME_TABLE = "20260610_Taxonomy_HighQActinos.csv"
GENOME_COL   = "Genome"

# --- Distance input options (the script uses the FIRST one it finds) ---------
# Option A1: numpy distance saved by 01_pcoa_gcf.py
NPY_DIST     = "pcoa_out/jaccard_distance.npy"
NPY_ORDER    = "pcoa_out/genome_order.txt"
# Option A2: a square CSV/TSV distance matrix exported from R (rownames=genomes)
#   In R:  write.csv(as.matrix(jac), "jac_matrix.csv")
CSV_DIST     = "20260611_jac_matrix.csv"          # set to "" to skip
# Option B: rebuild from the BGC table if no distance file is found
BGC_TABLE    = "bgc_table.tsv"
GCF_COL      = "gcf_id"
MIN_GENOMES_PER_GCF = 1                   # set 2 to drop singleton GCFs

# --- Test settings -----------------------------------------------------------
FACTORS      = ["Class", "Order", "Family", "Genus", "Clade"]  # ranks to test
PERMUTATIONS = 999
MIN_GROUP_N  = 25            # drop groups with fewer than this many genomes
RANDOM_SEED  = 42
MTC_METHOD   = "holm" #choose holm or bonferroni
RESULTS_TSV  = "permanova_results_python.tsv"
RESULTS_TXT  = "permanova_results_python.txt"
# ------------------------------------------------------------------------------


def log(m): print(f"[PERMANOVA] {m}", flush=True)


# ---------- load or build the distance matrix --------------------------------
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


# ---------- the test ----------------------------------------------------------
def run_one(dist, ids, meta, factor, permutations, min_group_n):
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
    pmv = permanova(dm, grouping=list(sub_grp), permutations=permutations)
    try:
        pdisp = permdisp(dm, grouping=list(sub_grp), permutations=permutations)
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
    np.random.seed(RANDOM_SEED)
    factors = [sys.argv[1]] if len(sys.argv) > 1 else FACTORS

    dist, ids = load_distance()
    log(f"Distance matrix: {dist.shape}, {len(ids):,} genomes, dtype={dist.dtype}")

    meta = pd.read_csv(GENOME_TABLE, sep=",", dtype=str)
    log(f"Metadata: {meta.shape[0]:,} genomes, columns include {list(meta.columns)[:8]}...")

    rows = [r for f in factors
            if (r := run_one(dist, ids, meta, f, PERMUTATIONS, MIN_GROUP_N))]

    if not rows:
        log("No valid tests produced."); return

    out = pd.DataFrame(rows)

#Corrections

    if MTC_METHOD and len(out) > 1:
        m = len(out)
        raw = out["permanova_p"].astype(float).values
        if MTC_METHOD.lower() == "bonferroni":
            adj = np.minimum(raw * m, 1.0)
        elif MTC_METHOD.lower() == "holm":
            order = np.argsort(raw)                 # smallest p first
            adj_sorted = np.empty(m, float)
            running = 0.0
            for rank, i in enumerate(order):        # step-down
                val = (m - rank) * raw[i]
                running = max(running, val)         # enforce monotonicity
                adj_sorted[i] = min(running, 1.0)
            adj = adj_sorted
        else:
            raise ValueError(f"Unknown MTC_METHOD: {MTC_METHOD}")
        out["permanova_p_adj"] = adj
        out["p_adj_method"]    = MTC_METHOD.lower()
        log(f"Applied {MTC_METHOD} correction across {m} factors -> column 'permanova_p_adj'.")
    elif MTC_METHOD:
        log("Only one factor tested; no multiple-comparison correction applied.")
    out.to_csv(RESULTS_TSV, sep="\t", index=False)
    with open(RESULTS_TXT, "w") as fh:
        fh.write("PERMANOVA of GCF composition (Jaccard) by taxonomic rank\n")
        fh.write("=" * 60 + "\n\n")
        fh.write(out.to_string(index=False))
        fh.write("\n\nInterpretation:\n")
        fh.write(" - permanova_p significant + permdisp_p NOT significant\n")
        fh.write("   => real difference in GCF COMPOSITION (location).\n")
        fh.write(" - higher permanova_R2 = that rank explains more variance.\n")
    log(f"Saved -> {RESULTS_TSV} and {RESULTS_TXT}")
    print("\n" + out.to_string(index=False))


if __name__ == "__main__":
    main()
