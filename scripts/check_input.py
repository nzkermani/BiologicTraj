"""Check that an .h5ad file is ready for the BiologicTraj pipeline.

Run this before scripts/run_pipeline.sh on real data. It reads the file, checks
the counts, gene names and cell metadata the pipeline needs, and writes a short
report of counts per sample and cell type. The report holds no cell-level data.

Required in adata.obs (one row per cell):
  sample        sample ID, one per patient x visit (e.g. P1_M0)
  patient       pseudonymised patient ID
  drug_class    e.g. anti-IL-5, anti-IL-4Ra (constant within a patient)
  visit         M0, M1 or M6
  visit_month   0, 1 or 6 (numeric, matching visit)
  cell_type     cell type label
Recommended (per sample, repeated on every cell of that sample):
  drug, processing_batch, processing_date, hours_to_freeze, sequencing_batch
adata.X must hold raw UMI counts; adata.var_names must be gene symbols (FOS, CD3E, ...).

Usage:
  python scripts/check_input.py --h5ad data/prism/prism.h5ad --out results/prism/input_check
Exit code 0 = ready (warnings may remain), 1 = must be fixed first.
"""

import argparse
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

REQUIRED = ["sample", "patient", "drug_class", "visit", "visit_month", "cell_type"]
RECOMMENDED = ["drug", "processing_batch", "processing_date", "hours_to_freeze", "sequencing_batch"]
VISITS = {"M0": 0, "M1": 1, "M6": 6}
STRESS_GENES = ["HSPA1A", "HSPA1B", "DNAJB1", "HSPH1", "HSP90AA1", "HSPA8",
                "FOS", "FOSB", "JUN", "JUNB", "IER2", "EGR1", "DUSP1", "ZFP36"]
MIN_CELLS = 20  # pseudobulk drops sample x cell type groups below this


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    a = ad.read_h5ad(args.h5ad)
    errors, warnings, notes = [], [], []
    notes.append(f"{a.n_obs:,} cells x {a.n_vars:,} genes")

    # counts
    block = a.X[: min(2000, a.n_obs)]
    vals = block.data if sp.issparse(block) else np.asarray(block).ravel()
    if vals.size and (vals.min() < 0 or not np.allclose(vals, np.round(vals))):
        hint = " (adata.raw is present: save adata.raw.to_adata() instead)" if a.raw is not None else ""
        errors.append("adata.X is not raw integer counts" + hint)

    # gene names
    ensembl = np.mean([bool(re.match(r"^ENS[A-Z]*G\d+", g)) for g in a.var_names[:500]])
    if ensembl > 0.5:
        errors.append("var_names look like Ensembl IDs; gene symbols are needed (map them first)")
    if not a.var_names.is_unique:
        errors.append("var_names are not unique (run adata.var_names_make_unique())")
    present = [g for g in STRESS_GENES if g in a.var_names]
    if len(present) < 10:
        warnings.append(f"only {len(present)}/{len(STRESS_GENES)} stress-check genes present: {', '.join(present)}")

    # metadata
    obs = a.obs
    missing = [c for c in REQUIRED if c not in obs]
    if missing:
        errors.append("missing required obs columns: " + ", ".join(missing))
    absent = [c for c in RECOMMENDED if c not in obs]
    if absent:
        warnings.append("recommended obs columns absent (the stress check is weaker without them): " + ", ".join(absent))

    if not missing:
        bad_visit = sorted(set(obs["visit"].astype(str)) - set(VISITS))
        if bad_visit:
            errors.append(f"visit must be M0, M1 or M6; found {bad_visit}")
        vm = obs[["visit", "visit_month"]].astype(str).drop_duplicates()
        wrong = vm[[str(VISITS.get(v, "?")) != str(m).replace(".0", "") for v, m in zip(vm["visit"], vm["visit_month"])]]
        if len(wrong):
            errors.append("visit and visit_month disagree: " + "; ".join(f"{v}/{m}" for v, m in wrong.values))
        per_sample = obs.groupby("sample", observed=True)[["patient", "visit", "drug_class"]].nunique()
        if (per_sample > 1).any().any():
            errors.append("some samples have more than one patient, visit or drug_class")
        if (obs.groupby("patient", observed=True)["drug_class"].nunique() > 1).any():
            errors.append("drug_class changes within a patient")
        design = obs[["patient", "drug_class", "visit"]].drop_duplicates()
        incomplete = design.groupby("patient", observed=True)["visit"].nunique()
        incomplete = incomplete[incomplete < 3]
        if len(incomplete):
            warnings.append("patients without all three visits: " + ", ".join(map(str, incomplete.index)))
        notes.append(f"{obs['patient'].nunique()} patients, {obs['sample'].nunique()} samples, "
                     f"{obs['cell_type'].nunique()} cell types")
        notes.append("patients per drug class: " + ", ".join(
            f"{k} {v}" for k, v in design.groupby("drug_class", observed=True)["patient"].nunique().items()))

        n = obs.groupby(["cell_type", "sample"], observed=True).size().unstack(fill_value=0)
        n.to_csv(out / "cells_per_celltype_sample.csv")
        thin = (n < MIN_CELLS).sum(axis=1)
        thin = thin[thin > 0]
        if len(thin):
            warnings.append(f"cell types with samples below {MIN_CELLS} cells (dropped there by pseudobulk): "
                            + ", ".join(f"{ct} ({k} samples)" for ct, k in thin.items()))

    status = "NOT READY" if errors else "READY"
    lines = [f"BiologicTraj input check: {status}", f"file: {args.h5ad}", ""]
    lines += ["Summary"] + [f"  {x}" for x in notes] + [""]
    lines += ["Errors (must fix)"] + [f"  - {x}" for x in errors or ["none"]] + [""]
    lines += ["Warnings"] + [f"  - {x}" for x in warnings or ["none"]]
    text = "\n".join(lines)
    print(text)
    (out / "input_check.txt").write_text(text + "\n")
    raise SystemExit(1 if errors else 0)


if __name__ == "__main__":
    main()
