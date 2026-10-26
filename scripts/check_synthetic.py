"""Sanity checks on the synthetic PRISM-like dataset.

1. Standard scanpy workflow -> UMAP coloured by cell type and visit.
2. Pseudobulk (mean counts per cell, per sample x cell type), averaged over
   patients within drug class -> empirical trajectory shapes (OTS), compared
   with the planted ground truth. This is the "empirical means" estimator that
   LongOmicsTraj accepts; a mixed model is the better estimator for real data.

Usage: python scripts/check_synthetic.py --data data/synthetic --figdir docs/figures
"""

import argparse
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

MIN_CELLS = 30  # minimum median cells per sample for a cell type to be scored


def pseudobulk_ots(adata, thr):
    rows = []
    X = adata.X
    for (sample, ct), idx in adata.obs.groupby(["sample", "cell_type"], observed=True).indices.items():
        o = adata.obs.iloc[idx[0]]
        mean = np.asarray(X[idx].mean(axis=0)).ravel()
        rows.append(pd.DataFrame({"gene": adata.var_names, "mean": mean, "sample": sample,
                                  "cell_type": ct, "patient": o["patient"],
                                  "drug_class": o["drug_class"], "visit_month": o["visit_month"],
                                  "n_cells": len(idx)}))
    pb = pd.concat(rows)
    vmm = (pb.groupby(["drug_class", "cell_type", "gene", "visit_month"], observed=True)["mean"].mean()
             .unstack("visit_month"))
    vmm = np.log2(vmm + 1e-3)
    d1, d2 = vmm[1] - vmm[0], vmm[6] - vmm[1]
    cls = lambda d: np.where(d > thr, "up", np.where(d < -thr, "down", "flat"))
    vmm["est_OTS"] = pd.Series(cls(d1), index=vmm.index) + "_" + pd.Series(cls(d2), index=vmm.index)
    return vmm.reset_index()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/synthetic")
    ap.add_argument("--figdir", default="docs/figures")
    args = ap.parse_args()
    data, figdir = Path(args.data), Path(args.figdir)
    figdir.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(data / "prism_like_synthetic.h5ad")
    thr = adata.uns["synthetic"]["ots_threshold_log2"]
    print(adata)

    # 1. standard workflow
    a = adata.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, n_top_genes=800)
    sc.pp.pca(a, n_comps=30, mask_var="highly_variable")
    sc.pp.neighbors(a, n_neighbors=15)
    sc.tl.umap(a, random_state=0)
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    sc.pl.umap(a, color="cell_type", ax=axes[0], show=False, title="Cell type (planted)")
    sc.pl.umap(a, color="visit", ax=axes[1], show=False, title="Visit")
    sc.pl.umap(a, color="cell_state", ax=axes[2], show=False, title="Cell state", legend_loc=None)
    fig.suptitle("SYNTHETIC PRISM-like PBMC data (not real measurements)", fontsize=13)
    fig.tight_layout()
    fig.savefig(figdir / "synthetic_umap.png", dpi=130)

    # 2. topology recovery
    est = pseudobulk_ots(adata, thr)
    truth = pd.read_csv(data / "ground_truth_gene_topology.csv")
    m = truth.merge(est[["drug_class", "cell_type", "gene", "est_OTS"]],
                    on=["drug_class", "cell_type", "gene"])
    # score only genes the ground truth marks as expressed, in cell types with enough cells
    ncell = adata.obs.groupby(["sample", "cell_type"], observed=True).size().groupby("cell_type", observed=True).median()
    ok_ct = ncell[ncell >= MIN_CELLS].index
    print(f"Cell types scored (median >= {MIN_CELLS} cells/sample): {len(ok_ct)}/{len(ncell)}")
    m = m[m["expressed"] & m["cell_type"].isin(ok_ct)]
    m["correct"] = m["true_OTS"] == m["est_OTS"]
    planted = m[m["true_OTS"] != "flat_flat"]
    print(f"\nExpressed gene x cell type x drug class rows: {len(m)}")
    print(f"Overall OTS agreement: {m['correct'].mean():.3f}")
    print(f"Planted (non-flat) recovery: {planted['correct'].mean():.3f} ({planted['correct'].sum()}/{len(planted)})")
    print("False non-flat calls among true flat_flat:",
          int(((m['true_OTS'] == 'flat_flat') & (m['est_OTS'] != 'flat_flat')).sum()))
    print("\nRecovery by planted shape:")
    print(planted.groupby(["drug_class", "true_OTS"])["correct"].agg(["mean", "size"]).round(2).to_string())

    # composition check: classical monocyte IL1B+ fraction by visit
    cm = adata.obs[adata.obs["cell_type"] == "Classical monocyte"]
    frac = (cm.assign(il1b=cm["cell_state"] == "IL1B+ CM")
              .groupby(["drug_class", "visit_month"], observed=True)["il1b"].mean().unstack())
    print("\nObserved IL1B+ fraction of classical monocytes:\n", frac.round(2).to_string())


if __name__ == "__main__":
    main()
