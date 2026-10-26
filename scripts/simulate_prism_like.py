"""Simulate a synthetic PBMC scRNA-seq dataset with the PRISM biologics design.

SYNTHETIC DATA ONLY. Nothing here is measured. The design (patients, drugs,
visits, PBMC cell types) mirrors Park et al., Exp Mol Med 2024, so the
BiologicTraj pipeline can be built and tested before the real data arrive.
The planted effects are loosely inspired by the published findings but their
sizes are invented; they exist so every analysis step has a known answer.

Design
------
8 patients: mepolizumab (P1, P2), reslizumab (P3, P4) -> anti-IL-5;
            dupilumab (P5-P8) -> anti-IL-4Ra
3 visits per patient: 0, 1 and 6 months  -> 24 samples
~3,700 cells per sample by default (PRISM: 89,728 cells / 24 samples)

Planted ground truth (log-scale effects relative to baseline)
--------------------------------------------------------------
NF-kB/AP-1 module, all cell types, both drugs   flat_down  (late suppression)
CD4 T cells shift from state S0 to S1 by 6 m    GIMAP genes flat_up
Classical monocytes shift IL1B+ -> S100A+ by 6 m (composition)
Dupilumab: STAT6 targets in B cells             down_flat  (early, sustained)
Dupilumab: interferon module in monocytes       up_down    (transient)
Anti-IL-5: CLC / IL5RA / HDC in basophils       down_flat; basophils fall
Anti-IL-5: cytotoxic module in CD8 T cells      down_up    (dip, recovery)
Everything else                                  flat_flat

The true topology of every gene x cell type x drug class is computed from the
noise-free generative means and written alongside the data.

Usage
-----
python scripts/simulate_prism_like.py --out data/synthetic --seed 42
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

VISITS = [0, 1, 6]  # months
PATIENTS = {
    "P1": "mepolizumab", "P2": "mepolizumab",
    "P3": "reslizumab", "P4": "reslizumab",
    "P5": "dupilumab", "P6": "dupilumab", "P7": "dupilumab", "P8": "dupilumab",
}
DRUG_CLASS = {"mepolizumab": "anti-IL-5", "reslizumab": "anti-IL-5", "dupilumab": "anti-IL-4Ra"}

# ---------------------------------------------------------------- genes
HOUSEKEEPING = ["ACTB", "GAPDH", "B2M", "MALAT1", "TMSB4X", "RPL13", "RPS27",
                "HLA-A", "HLA-B", "PTPRC", "EEF1A1", "TPT1"]
MARKERS = {
    "CD4 T naive": ["CD3D", "CD3E", "CD4", "IL7R", "CCR7", "SELL", "LEF1", "TCF7"],
    "CD4 T central memory": ["CD3D", "CD3E", "CD4", "IL7R", "CCR7", "S100A4", "AQP3"],
    "CD4 T effector memory": ["CD3D", "CD3E", "CD4", "IL7R", "GZMK", "CCL5", "S100A4"],
    "Treg": ["CD3D", "CD4", "FOXP3", "IL2RA", "CTLA4"],
    "CD8 T": ["CD3D", "CD3E", "CD8A", "CD8B", "CCL5"],
    "NK": ["GNLY", "NKG7", "KLRD1", "KLRF1", "FCGR3A"],
    "B naive": ["MS4A1", "CD79A", "CD79B", "CD19", "IGHD", "IGHM", "TCL1A", "BANK1"],
    "B memory": ["MS4A1", "CD79A", "CD79B", "CD19", "CD27", "BANK1"],
    "Classical monocyte": ["CD14", "LYZ", "VCAN", "FCN1", "S100A8", "S100A9"],
    "Non-classical monocyte": ["FCGR3A", "MS4A7", "CDKN1C", "LST1", "LYZ"],
    "cDC": ["FCER1A", "CD1C", "CLEC10A", "HLA-DQA1"],
    "pDC": ["LILRA4", "IL3RA", "CLEC4C", "TCF4"],
    "Basophil": ["HDC", "GATA2", "CPA3", "MS4A2", "ENPP3"],
    "Platelet": ["PPBP", "PF4", "GP9", "ITGA2B"],
}
MODULES = {
    "nfkb_ap1": ["FOS", "FOSB", "JUN", "JUNB", "JUND", "DUSP1", "DUSP2", "NFKBIA",
                 "NFKBIZ", "ZFP36", "TNFAIP3", "CD69", "PPP1R15A", "SOCS3", "BTG2",
                 "IER2", "KLF6", "EGR1"],
    "gimap": ["GIMAP1", "GIMAP4", "GIMAP5", "GIMAP7"],
    "il1b_cm": ["IL1B", "CXCL8", "CCL3", "CCL4", "TNF"],
    "s100_cm": ["S100A12"],  # S100A8/9 are already CM markers
    "stat6_b": ["FCER2", "IL4R", "IGHE", "CISH", "SOCS1"],
    "ifn": ["ISG15", "IFI6", "IFI44L", "MX1", "IFIT1", "IFIT3", "OAS1"],
    "il5_baso": ["CLC", "IL5RA"],
    "cytotoxic": ["GZMB", "GZMH", "PRF1"],
    "mhc2": ["HLA-DRA", "HLA-DRB1", "CD74"],
}
N_FILLER = 1500  # synthetic background genes, always flat
EXPRESSED_MIN = 0.5  # mean UMI per cell; below this a gene is not scored in the ground truth

# CD4 cell types that carry the S0 -> S1 state switch
CD4_SWITCH = ["CD4 T naive", "CD4 T central memory", "CD4 T effector memory"]

BASE_PROPORTIONS = {
    "CD4 T naive": 0.15, "CD4 T central memory": 0.12, "CD4 T effector memory": 0.06,
    "Treg": 0.03, "CD8 T": 0.12, "NK": 0.10, "B naive": 0.07, "B memory": 0.04,
    "Classical monocyte": 0.18, "Non-classical monocyte": 0.04, "cDC": 0.015,
    "pDC": 0.008, "Basophil": 0.012, "Platelet": 0.01,
}

# fraction of CD4 cells in state S1, and of classical monocytes that are IL1B+, by visit
S1_FRACTION = {0: 0.10, 1: 0.12, 6: 0.70}
IL1B_CM_FRACTION = {0: 0.60, 1: 0.58, 6: 0.20}


def build_gene_table(rng: np.random.Generator) -> pd.DataFrame:
    rows = []
    seen = set()

    def add(g, module):
        if g not in seen:
            seen.add(g)
            rows.append({"gene": g, "module": module, "synthetic_name": False})

    for g in HOUSEKEEPING:
        add(g, "housekeeping")
    for genes in MARKERS.values():
        for g in genes:
            add(g, "marker")
    for m, genes in MODULES.items():
        for g in genes:
            add(g, m)
    for i in range(N_FILLER):
        rows.append({"gene": f"SYNG{i + 1:04d}", "module": "background", "synthetic_name": True})
    return pd.DataFrame(rows).set_index("gene")


def baseline_log_means(genes: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    """Log mean UMI per gene for each cell type at baseline (state S0 / IL1B+)."""
    cts = list(BASE_PROPORTIONS)
    mu = pd.DataFrame(-3.0, index=genes.index, columns=cts)
    bg = genes.index[genes["module"] == "background"]
    mu.loc[bg, :] = rng.normal(-1.5, 1.0, size=len(bg))[:, None] + rng.normal(0, 0.3, size=(len(bg), len(cts)))
    mu.loc[HOUSEKEEPING, :] = rng.uniform(1.5, 4.0, size=len(HOUSEKEEPING))[:, None]
    for ct, genes_ct in MARKERS.items():
        mu.loc[genes_ct, ct] = rng.uniform(1.0, 2.8, size=len(genes_ct))
    # AP-1/NF-kB module: expressed everywhere, higher in myeloid
    mu.loc[MODULES["nfkb_ap1"], :] = rng.uniform(0.0, 1.2, size=(len(MODULES["nfkb_ap1"]), len(cts)))
    mu.loc[MODULES["nfkb_ap1"], ["Classical monocyte", "Non-classical monocyte"]] += 0.8
    mu.loc[MODULES["gimap"], [c for c in cts if "T" in c or c == "NK"]] = 0.8
    mu.loc[MODULES["il1b_cm"], "Classical monocyte"] = 0.5  # raised further in IL1B+ state
    mu.loc[MODULES["s100_cm"], "Classical monocyte"] = 1.0
    mu.loc[MODULES["stat6_b"], ["B naive", "B memory"]] = 1.2
    mu.loc[MODULES["ifn"], :] = 0.2
    mu.loc[MODULES["il5_baso"], "Basophil"] = 2.0
    mu.loc[MODULES["cytotoxic"], ["CD8 T", "NK"]] = 2.0
    mu.loc[MODULES["mhc2"], ["B naive", "B memory", "Classical monocyte",
                             "Non-classical monocyte", "cDC", "pDC"]] = 2.5
    return mu


def state_effects(mu: pd.DataFrame) -> dict:
    """Log-scale shifts for substates relative to the baseline profile."""
    s1 = pd.Series(0.0, index=mu.index)
    s1[MODULES["gimap"]] = 0.9
    s1[["FOS", "FOSB", "JUN", "JUNB", "DUSP1", "CD69"]] = -0.9
    il1b = pd.Series(0.0, index=mu.index)
    il1b[MODULES["il1b_cm"]] = 1.8
    il1b[MODULES["s100_cm"] + ["S100A8", "S100A9"]] = -0.8
    s100 = pd.Series(0.0, index=mu.index)
    s100[MODULES["s100_cm"] + ["S100A8", "S100A9"]] = 0.8
    return {"S1": s1, "IL1B+": il1b, "S100A+": s100}


def treatment_effect(genes: pd.Index, cell_type: str, drug_class: str, visit: int) -> pd.Series:
    """Log-scale treatment effect on gene means for one cell type, drug class and visit."""
    eff = pd.Series(0.0, index=genes)
    v = {0: 0, 1: 1, 6: 2}[visit]
    # late NF-kB/AP-1 suppression, all cells, both drugs
    eff[MODULES["nfkb_ap1"]] += [0.0, 0.0, -1.0][v]
    if drug_class == "anti-IL-4Ra":
        if cell_type in ("B naive", "B memory"):
            eff[MODULES["stat6_b"]] += [0.0, -1.2, -1.2][v]
        if cell_type in ("Classical monocyte", "Non-classical monocyte"):
            eff[MODULES["ifn"]] += [0.0, 1.2, 0.0][v]
    else:  # anti-IL-5
        if cell_type == "Basophil":
            eff[MODULES["il5_baso"] + ["HDC"]] += [0.0, -1.0, -1.0][v]
        if cell_type == "CD8 T":
            eff[MODULES["cytotoxic"]] += [0.0, -0.9, 0.0][v]
    return eff


def visit_proportions(drug_class: str, visit: int) -> pd.Series:
    p = pd.Series(BASE_PROPORTIONS, dtype=float)
    if drug_class == "anti-IL-5" and visit > 0:
        p["Basophil"] *= 0.35
    return p / p.sum()


def classify(delta: np.ndarray, thr: float) -> np.ndarray:
    return np.where(delta > thr, "up", np.where(delta < -thr, "down", "flat"))


def simulate(seed: int = 42, cells_per_sample: int = 3700, theta: float = 6.0,
             patient_sd: float = 0.15, ots_threshold: float = 0.35):
    rng = np.random.default_rng(seed)
    genes = build_gene_table(rng)
    mu0 = baseline_log_means(genes, rng)
    states = state_effects(mu0)
    gidx = genes.index

    patient_gene = {p: rng.normal(0, patient_sd, size=len(gidx)) for p in PATIENTS}
    patient_comp = {p: rng.dirichlet(np.full(len(BASE_PROPORTIONS), 200.0)) for p in PATIENTS}

    blocks, obs_rows, truth_rows = [], [], []
    for patient, drug in PATIENTS.items():
        dclass = DRUG_CLASS[drug]
        for visit in VISITS:
            sample = f"{patient}_M{visit}"
            props = visit_proportions(dclass, visit).values * patient_comp[patient] * len(BASE_PROPORTIONS)
            props = props / props.sum()
            n = rng.poisson(cells_per_sample)
            cts = rng.choice(list(BASE_PROPORTIONS), size=n, p=props)
            for ct in BASE_PROPORTIONS:
                k = int((cts == ct).sum())
                if k == 0:
                    continue
                base = mu0[ct].values + patient_gene[patient] + treatment_effect(gidx, ct, dclass, visit).values
                if ct in CD4_SWITCH:
                    is_s1 = rng.random(k) < S1_FRACTION[visit]
                    sub = np.where(is_s1, "S1", "S0")
                elif ct == "Classical monocyte":
                    is_il1b = rng.random(k) < IL1B_CM_FRACTION[visit]
                    sub = np.where(is_il1b, "IL1B+", "S100A+")
                else:
                    sub = np.full(k, "")
                logmu = np.tile(base, (k, 1))
                for s, shift in states.items():
                    m = sub == s
                    if m.any():
                        logmu[m] += shift.values
                size = rng.lognormal(0.0, 0.3, size=(k, 1))
                mean = np.exp(logmu) * size
                counts = rng.negative_binomial(theta, theta / (theta + mean))
                blocks.append(sp.csr_matrix(counts.astype(np.int32)))
                label = ct if ct not in ("Classical monocyte",) else None
                for s in sub:
                    obs_rows.append({
                        "patient": patient, "drug": drug, "drug_class": dclass,
                        "visit": f"M{visit}", "visit_month": visit, "sample": sample,
                        "cell_type": ct,
                        "cell_state": (f"CD4 {s}" if ct in CD4_SWITCH else
                                       (f"{s} CM" if ct == "Classical monocyte" else ct)),
                    })

    X = sp.vstack(blocks).tocsr()
    obs = pd.DataFrame(obs_rows)
    obs.index = [f"cell{i:06d}" for i in range(len(obs))]
    for c in ["patient", "drug", "drug_class", "visit", "sample", "cell_type", "cell_state"]:
        obs[c] = obs[c].astype("category")
    adata = ad.AnnData(X=X, obs=obs, var=genes.copy())
    adata.uns["synthetic"] = {
        "description": "SYNTHETIC data mimicking the PRISM scRNA-seq design (Park et al. 2024). Not real measurements.",
        "seed": seed, "cells_per_sample": cells_per_sample, "nb_theta": theta,
        "patient_sd": patient_sd, "ots_threshold_log2": ots_threshold,
    }

    # ---------------- ground truth from noise-free means (pseudobulk scale)
    for dclass in ["anti-IL-5", "anti-IL-4Ra"]:
        pts = [p for p, d in PATIENTS.items() if DRUG_CLASS[d] == dclass]
        for ct in BASE_PROPORTIONS:
            means = []
            for visit in VISITS:
                acc = np.zeros(len(gidx))
                for p in pts:
                    base = mu0[ct].values + patient_gene[p] + treatment_effect(gidx, ct, dclass, visit).values
                    if ct in CD4_SWITCH:
                        f = S1_FRACTION[visit]
                        m = (1 - f) * np.exp(base) + f * np.exp(base + states["S1"].values)
                    elif ct == "Classical monocyte":
                        f = IL1B_CM_FRACTION[visit]
                        m = f * np.exp(base + states["IL1B+"].values) + (1 - f) * np.exp(base + states["S100A+"].values)
                    else:
                        m = np.exp(base)
                    acc += m
                means.append(np.log2(acc / len(pts) + 1e-3))
            means = np.vstack(means).T  # genes x visits
            d1, d2 = means[:, 1] - means[:, 0], means[:, 2] - means[:, 1]
            ots = np.char.add(np.char.add(classify(d1, ots_threshold), "_"), classify(d2, ots_threshold))
            expressed = means.max(axis=1) > np.log2(EXPRESSED_MIN)  # score only genes with real signal in this cell type
            for g, o, e, row in zip(gidx, ots, expressed, means):
                truth_rows.append({"gene": g, "module": genes.loc[g, "module"], "cell_type": ct,
                                   "drug_class": dclass, "true_OTS": o, "expressed": bool(e),
                                   "log2_mean_M0": row[0], "log2_mean_M1": row[1], "log2_mean_M6": row[2]})
    truth = pd.DataFrame(truth_rows)

    comp_rows = []
    for dclass in ["anti-IL-5", "anti-IL-4Ra"]:
        for visit in VISITS:
            p = visit_proportions(dclass, visit)
            for ct, v in p.items():
                comp_rows.append({"drug_class": dclass, "visit_month": visit, "cell_type": ct, "expected_fraction": v})
            for v_name, frac in (("CD4 S1 (within CD4 switching types)", S1_FRACTION[visit]),
                                 ("IL1B+ CM (within classical monocytes)", IL1B_CM_FRACTION[visit])):
                comp_rows.append({"drug_class": dclass, "visit_month": visit, "cell_type": v_name, "expected_fraction": frac})
    composition = pd.DataFrame(comp_rows)
    return adata, truth, composition


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="data/synthetic")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--cells-per-sample", type=int, default=3700)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    adata, truth, composition = simulate(seed=args.seed, cells_per_sample=args.cells_per_sample)
    adata.write_h5ad(out / "prism_like_synthetic.h5ad", compression="gzip")
    truth.to_csv(out / "ground_truth_gene_topology.csv", index=False)
    composition.to_csv(out / "ground_truth_composition.csv", index=False)
    print(f"cells={adata.n_obs} genes={adata.n_vars} samples={adata.obs['sample'].nunique()}")
    print(truth.query("expressed and true_OTS != 'flat_flat'")
          .groupby(["drug_class", "true_OTS"]).size().to_string())


if __name__ == "__main__":
    main()
