"""Score LongOmicsTraj trajectory shapes (OTS) against the synthetic ground truth.

Scored rows: gene x cell type x drug class marked `expressed` in the ground
truth and analysed by LongOmicsTraj. For every estimator x threshold reports
recovery of planted (non-flat) shapes and false non-flat calls among true
flat_flat rows.

With --transitions (limma output from scripts/normalise_pseudobulk.R) each
row is also scored after a support filter: a transition called up/down is set
to flat unless its BH-adjusted p-value is below --alpha. This is reported as
threshold "<threshold>+support".

Comparator-aware queries (from run_longomicstraj.R) are scored against the
genes the truth says should match, with and without the support filter.

Usage: python scripts/score_longomicstraj.py --res results/synthetic/longomicstraj_voom \
         --truth-file data/synthetic/ground_truth_gene_topology.csv \
         --transitions results/synthetic/pseudobulk_voom/transitions.csv
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

QUERY_TRUTH = {
    "shared delayed suppression": {"anti-IL-5": "flat_down", "anti-IL-4Ra": "flat_down"},
    "early sustained fall, anti-IL-4Ra only": {"anti-IL-4Ra": "down_flat"},
    "transient rise, anti-IL-4Ra only": {"anti-IL-4Ra": "up_down", "anti-IL-5": "flat_flat"},
    "dip and recovery, anti-IL-5 only": {"anti-IL-5": "down_up", "anti-IL-4Ra": "flat_flat"},
}
KEY = ["gene", "cell_type", "drug_class"]


def apply_support(ots: pd.DataFrame, trans: pd.DataFrame, alpha: float) -> pd.DataFrame:
    p = trans.pivot_table(index=KEY, columns="transition", values="adj_p").reset_index()
    o = ots.merge(p, on=KEY, how="left")
    parts = o["est_OTS"].str.split("_", expand=True)
    for j, tcol in enumerate(["M0->M1", "M1->M6"]):
        unsupported = ~(o[tcol] < alpha)  # NaN (no test) counts as unsupported
        parts[j] = np.where(unsupported, "flat", parts[j])
    o["est_OTS"] = parts[0] + "_" + parts[1]
    o["threshold"] = o["threshold"] + "+support"
    return o[ots.columns]


def score(m: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (est, thr), g in m.groupby(["estimator", "threshold"]):
        p, f = g[g["planted"]], g[~g["planted"]]
        rows.append({"estimator": est, "threshold": thr, "rows_scored": len(g),
                     "planted_recovered": f"{p['correct'].sum()}/{len(p)}",
                     "planted_recovery": round(p["correct"].mean(), 3) if len(p) else np.nan,
                     "false_nonflat": int((f["est_OTS"] != "flat_flat").sum()),
                     "false_nonflat_rate": round((f["est_OTS"] != "flat_flat").mean(), 3)})
    return pd.DataFrame(rows)


def score_queries(ots_q: pd.DataFrame, truth: pd.DataFrame, label: str) -> pd.DataFrame:
    """Re-run the four comparator queries on an OTS table (one estimator x threshold)."""
    w = ots_q.pivot_table(index=["cell_type", "gene"], columns="drug_class", values="est_OTS", aggfunc="first")
    t = truth.pivot_table(index=["cell_type", "gene"], columns="drug_class", values="true_OTS", aggfunc="first")
    rows = []
    for q, cond in QUERY_TRUTH.items():
        got, want = w.copy(), t.copy()
        for grp, shape in cond.items():
            got = got[got.get(grp) == shape] if grp in got else got.iloc[0:0]
            want = want[want.get(grp) == shape] if grp in want else want.iloc[0:0]
        if q == "early sustained fall, anti-IL-4Ra only" and "anti-IL-5" in w:
            got = got[got["anti-IL-5"] != "down_flat"]
        gs, ws = set(got.index), set(want.index)
        tp = len(gs & ws)
        rows.append({"filter": label, "query": q, "truth_matches": len(ws), "returned": len(gs),
                     "true_positives": tp, "precision": round(tp / len(gs), 3) if gs else np.nan,
                     "recall": round(tp / len(ws), 3) if ws else np.nan})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res", default="results/synthetic/longomicstraj_voom")
    ap.add_argument("--truth-file", default="data/synthetic/ground_truth_gene_topology.csv")
    ap.add_argument("--transitions", default=None)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--query-threshold", default="fixed_0.35")
    args = ap.parse_args()
    res = Path(args.res)

    ots = pd.read_csv(res / "ots_all.csv").rename(
        columns={"object_id": "gene", "group": "drug_class", "topology_label": "est_OTS"})
    ots = ots[KEY + ["estimator", "threshold", "est_OTS"]]
    if args.transitions:
        ots = pd.concat([ots, apply_support(ots, pd.read_csv(args.transitions), args.alpha)])
    truth = pd.read_csv(args.truth_file)
    truth = truth[truth["expressed"]]
    m = truth.merge(ots, on=KEY)
    m["correct"] = m["true_OTS"] == m["est_OTS"]
    m["planted"] = m["true_OTS"] != "flat_flat"

    summary = score(m).sort_values(["estimator", "threshold"])
    by_shape = (m[m["planted"]].groupby(["estimator", "threshold", "drug_class", "true_OTS"])["correct"]
                .agg(recovered="sum", total="size").reset_index())
    stress = None
    if "stress_gene" in m:
        s = m[m["stress_gene"] & (m["true_OTS"] == "flat_flat")]
        stress = (s.groupby(["estimator", "threshold"])
                   .apply(lambda g: pd.Series({"stress_genes_truly_flat": len(g),
                                               "called_changing": int((g["est_OTS"] != "flat_flat").sum())}),
                          include_groups=False).reset_index())

    qs = []
    for thr in [args.query_threshold, args.query_threshold + "+support"]:
        sub = ots[(ots["estimator"] == "lmm") & (ots["threshold"] == thr)]
        if len(sub):
            qs.append(score_queries(sub, truth, f"lmm {thr}"))
    qsum = pd.concat(qs) if qs else pd.DataFrame()

    analysed = ots[KEY].drop_duplicates()
    planted_truth = truth[truth["true_OTS"] != "flat_flat"]
    na = planted_truth.merge(analysed, how="left", indicator=True)
    na = na[na["_merge"] == "left_only"][["cell_type", "drug_class", "gene", "true_OTS"]]

    summary.to_csv(res / "score_summary.csv", index=False)
    by_shape.to_csv(res / "score_by_shape.csv", index=False)
    qsum.to_csv(res / "score_queries.csv", index=False)
    if stress is not None:
        stress.to_csv(res / "score_stress_genes.csv", index=False)
    pd.set_option("display.width", 220)
    print(summary.to_string(index=False))
    print(f"\nPlanted rows not analysed: {len(na)}"
          + (f" ({', '.join(sorted(set(na['cell_type'] + ' / ' + na['drug_class'])))})" if len(na) else ""))
    if stress is not None and len(stress):
        print("\nStress genes that are biologically flat but called changing:")
        print(stress.to_string(index=False))
    print("\nQueries:")
    print(qsum.to_string(index=False))


if __name__ == "__main__":
    main()
