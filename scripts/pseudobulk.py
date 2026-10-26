"""Pseudobulk an AnnData object by sample x cell type.

--mode sum   (recommended) summed UMI counts per gene and sample, for a count
             model (edgeR/limma-voom, scripts/normalise_pseudobulk.R).
--mode mean  log2(mean UMI per cell + 1e-3): the scale of the synthetic ground
             truth; quick, but no library-size normalisation or count model.

Writes one genes x samples CSV per cell type, a sample table and an index.

Usage: python scripts/pseudobulk.py --h5ad data/synthetic/prism_like_synthetic.h5ad \
         --out results/synthetic/pseudobulk --min-cells 20
"""

import argparse
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


def slug(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", s).strip("_")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", default="data/synthetic/prism_like_synthetic.h5ad")
    ap.add_argument("--out", default="results/synthetic/pseudobulk")
    ap.add_argument("--min-cells", type=int, default=20, help="drop sample x cell type groups with fewer cells")
    ap.add_argument("--offset", type=float, default=1e-3)
    ap.add_argument("--mode", choices=["sum", "mean"], default="sum")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(args.h5ad)
    obs = adata.obs

    meta_cols = [c for c in ["sample", "patient", "drug", "drug_class", "visit", "visit_month",
                             "processing_batch"] if c in obs]
    samples = (obs[meta_cols]
               .drop_duplicates().set_index("sample").sort_index())
    samples.to_csv(out / "samples.csv")

    index = []
    for ct, idx_ct in obs.groupby("cell_type", observed=True).indices.items():
        cols, names, ncell = [], [], {}
        for sample, idx in obs.iloc[idx_ct].groupby("sample", observed=True).indices.items():
            rows = idx_ct[idx]
            ncell[sample] = len(rows)
            if len(rows) < args.min_cells:
                continue
            agg = adata.X[rows].sum(axis=0) if args.mode == "sum" else adata.X[rows].mean(axis=0)
            cols.append(np.asarray(agg).ravel())
            names.append(sample)
        if not cols:
            continue
        vals = np.vstack(cols).T
        if args.mode == "mean":
            vals = np.log2(vals + args.offset)
        else:
            vals = vals.astype(np.int64)
        mat = pd.DataFrame(vals, index=adata.var_names, columns=names)
        f = out / f"{slug(ct)}.csv"
        mat.to_csv(f)
        index.append({"cell_type": ct, "file": f.name, "n_samples": len(names), "mode": args.mode,
                      "median_cells": int(np.median(list(ncell.values())))})
    pd.DataFrame(index).to_csv(out / "index.csv", index=False)
    print(pd.DataFrame(index).to_string(index=False))


if __name__ == "__main__":
    main()
