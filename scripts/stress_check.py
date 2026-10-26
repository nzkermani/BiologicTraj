"""Check whether a visit effect looks like a sample-processing stress artefact.

Processing and dissociation stress raises immediate-early genes (FOS, JUN,
EGR1, DUSP1, ...) together with heat-shock genes (HSPA1A, HSPA1B, DNAJB1, ...)
in every cell type (van den Brink et al., Nat Methods 2017). A biologic acting
on NF-kB/AP-1 signalling would not be expected to move heat-shock genes. So:

  immediate-early genes change, heat-shock genes do not    -> consistent with biology
  both change together, same direction, across cell types  -> suspect processing artefact

For each sample x cell type the script computes a heat-shock score and an
immediate-early score (mean log-normalised expression of each gene set), then
for each visit transition and drug class the within-patient change in both
scores. A transition is flagged when the heat-shock change exceeds
--min-change in the same direction as the immediate-early change in at least
--min-fraction of cell types. Scores are also correlated across samples
after removing each patient's mean.

This is a screen, not proof: always compare against processing dates,
time-to-freeze and sequencing batch, which the real data need.

Usage: python scripts/stress_check.py --h5ad data/synthetic/prism_like_synthetic.h5ad \
         --out results/synthetic/stress_check
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

HEAT_SHOCK = ["HSPA1A", "HSPA1B", "DNAJB1", "HSPH1", "HSP90AA1", "HSPA8"]
IEG = ["FOS", "FOSB", "JUN", "JUNB", "IER2", "EGR1", "DUSP1", "ZFP36"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", default="data/synthetic/prism_like_synthetic.h5ad")
    ap.add_argument("--out", default="results/synthetic/stress_check")
    ap.add_argument("--min-cells", type=int, default=20)
    ap.add_argument("--min-change", type=float, default=0.15, help="log-scale change in heat-shock score to count")
    ap.add_argument("--min-fraction", type=float, default=0.75, help="fraction of cell types needed to flag")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    a = ad.read_h5ad(args.h5ad)
    hs = [g for g in HEAT_SHOCK if g in a.var_names]
    ieg = [g for g in IEG if g in a.var_names]
    if len(hs) < 3 or len(ieg) < 3:
        raise SystemExit(f"too few marker genes present (heat-shock {len(hs)}, immediate-early {len(ieg)})")
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    X = a.X
    a.obs["hs_score"] = np.asarray(X[:, [a.var_names.get_loc(g) for g in hs]].mean(axis=1)).ravel()
    a.obs["ieg_score"] = np.asarray(X[:, [a.var_names.get_loc(g) for g in ieg]].mean(axis=1)).ravel()

    keys = ["sample", "patient", "drug_class", "visit_month", "cell_type"]
    extra = [c for c in ("processing_batch",) if c in a.obs]
    g = a.obs.groupby(keys + extra, observed=True)
    per = g[["hs_score", "ieg_score"]].mean()
    per["n_cells"] = g.size()
    per = per[per["n_cells"] >= args.min_cells].reset_index()
    per.to_csv(out / "stress_scores_per_sample_celltype.csv", index=False)

    # within-patient change per transition, then mean over patients
    visits = sorted(per["visit_month"].unique())
    rows = []
    for (dc, ct), d in per.groupby(["drug_class", "cell_type"]):
        wide = d.pivot_table(index="patient", columns="visit_month", values=["hs_score", "ieg_score"])
        for v0, v1 in zip(visits[:-1], visits[1:]):
            if (("hs_score", v0) not in wide) or (("hs_score", v1) not in wide):
                continue
            dhs = (wide[("hs_score", v1)] - wide[("hs_score", v0)]).dropna()
            dieg = (wide[("ieg_score", v1)] - wide[("ieg_score", v0)]).dropna()
            if len(dhs) == 0:
                continue
            rows.append({"drug_class": dc, "cell_type": ct, "transition": f"M{v0}->M{v1}",
                         "n_patients": len(dhs), "delta_heat_shock": dhs.mean(), "delta_ieg": dieg.mean()})
    tr = pd.DataFrame(rows)
    tr["co_moving"] = (tr["delta_heat_shock"].abs() >= args.min_change) & \
                      (np.sign(tr["delta_heat_shock"]) == np.sign(tr["delta_ieg"]))
    tr.to_csv(out / "stress_transitions.csv", index=False)

    flags = (tr.groupby(["drug_class", "transition"])
               .agg(cell_types=("co_moving", "size"), co_moving=("co_moving", "sum"),
                    mean_delta_heat_shock=("delta_heat_shock", "mean"), mean_delta_ieg=("delta_ieg", "mean"))
               .reset_index())
    flags["fraction_co_moving"] = flags["co_moving"] / flags["cell_types"]
    flags["flag"] = np.where(flags["fraction_co_moving"] >= args.min_fraction,
                             "SUSPECT processing artefact", "ok")

    # correlation of heat-shock and immediate-early scores across samples, patient means removed
    cors = []
    for ct, d in per.groupby("cell_type"):
        d = d.copy()
        for c in ("hs_score", "ieg_score"):
            d[c] = d[c] - d.groupby("patient")[c].transform("mean")
        if len(d) >= 6:
            cors.append({"cell_type": ct, "r_heatshock_vs_ieg": d["hs_score"].corr(d["ieg_score"]), "n_samples": len(d)})
    cors = pd.DataFrame(cors)
    flags.to_csv(out / "stress_flags.csv", index=False)
    cors.to_csv(out / "stress_correlations.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    for ax, col, title in ((axes[0], "hs_score", "Heat-shock score"), (axes[1], "ieg_score", "Immediate-early score")):
        s = per.groupby(["sample", "patient", "visit_month"], observed=True)[col].mean().reset_index()
        for p, d in s.groupby("patient"):
            ax.plot(d["visit_month"], d[col], marker="o", color="0.6", lw=1)
        m = s.groupby("visit_month")[col].mean()
        ax.plot(m.index, m.values, marker="o", color="#B45309", lw=3, label="mean")
        ax.set_xticks(visits); ax.set_xlabel("Month"); ax.set_title(title + " per sample (mean over cell types)")
        ax.legend(frameon=False)
    fig.suptitle(f"Processing-stress check: {Path(args.h5ad).name}")
    fig.tight_layout()
    fig.savefig(out / "stress_check.png", dpi=120)

    pd.set_option("display.width", 200)
    print(flags.round(3).to_string(index=False))
    print(f"\nMedian within-patient correlation of heat-shock vs immediate-early scores: "
          f"{cors['r_heatshock_vs_ieg'].median():.2f}")


if __name__ == "__main__":
    main()
