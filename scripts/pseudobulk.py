"""Pseudobulk an AnnData object by sample x cell type for LongOmicsTraj.

Value per gene and sample: log2(mean UMI per cell + 1e-3), i.e. the same scale
as the synthetic ground truth. Writes one genes x samples CSV per cell type and
one sample metadata table.

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
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(args.h5ad)
    obs = adata.obs

    samples = (obs[["sample", "patient", "drug", "drug_class", "visit", "visit_month"]]
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
            cols.append(np.asarray(adata.X[rows].mean(axis=0)).ravel())
            names.append(sample)
        if not cols:
            continue
        mat = pd.DataFrame(np.log2(np.vstack(cols).T + args.offset), index=adata.var_names, columns=names)
        f = out / f"{slug(ct)}.csv"
        mat.to_csv(f)
        index.append({"cell_type": ct, "file": f.name, "n_samples": len(names),
                      "median_cells": int(np.median(list(ncell.values())))})
    pd.DataFrame(index).to_csv(out / "index.csv", index=False)
    print(pd.DataFrame(index).to_string(index=False))


if __name__ == "__main__":
    main()
