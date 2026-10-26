"""Score LongOmicsTraj trajectory shapes (OTS) against the synthetic ground truth.

Scored rows: gene x cell type x drug class marked `expressed` in the ground
truth and analysed by LongOmicsTraj. Reports, per estimator x threshold:
recovery of planted (non-flat) shapes, false non-flat calls among true
flat_flat rows, and recovery per planted shape. Also scores the
comparator-aware queries against the genes the truth says should match.

Usage: python scripts/score_longomicstraj.py --res results/synthetic/longomicstraj --truth data/synthetic
"""

import argparse
from pathlib import Path

import pandas as pd

QUERY_TRUTH = {
    "shared delayed suppression": {"anti-IL-5": "flat_down", "anti-IL-4Ra": "flat_down"},
    "early sustained fall, anti-IL-4Ra only": {"anti-IL-4Ra": "down_flat"},
    "transient rise, anti-IL-4Ra only": {"anti-IL-4Ra": "up_down", "anti-IL-5": "flat_flat"},
    "dip and recovery, anti-IL-5 only": {"anti-IL-5": "down_up", "anti-IL-4Ra": "flat_flat"},
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--res", default="results/synthetic/longomicstraj")
    ap.add_argument("--truth", default="data/synthetic")
    args = ap.parse_args()
    res, tdir = Path(args.res), Path(args.truth)

    ots = pd.read_csv(res / "ots_all.csv").rename(
        columns={"object_id": "gene", "group": "drug_class", "topology_label": "est_OTS"})
    truth = pd.read_csv(tdir / "ground_truth_gene_topology.csv")
    m = truth[truth["expressed"]].merge(ots, on=["gene", "cell_type", "drug_class"])
    m["correct"] = m["true_OTS"] == m["est_OTS"]
    m["planted"] = m["true_OTS"] != "flat_flat"

    rows = []
    for (est, thr), g in m.groupby(["estimator", "threshold"]):
        p, f = g[g["planted"]], g[~g["planted"]]
        rows.append({"estimator": est, "threshold": thr, "rows_scored": len(g),
                     "planted_recovered": f"{p['correct'].sum()}/{len(p)}",
                     "planted_recovery": round(p["correct"].mean(), 3),
                     "false_nonflat": int((f["est_OTS"] != "flat_flat").sum()),
                     "false_nonflat_rate": round((f["est_OTS"] != "flat_flat").mean(), 3),
                     "overall_agreement": round(g["correct"].mean(), 3)})
    summary = pd.DataFrame(rows)
    by_shape = (m[m["planted"]].groupby(["estimator", "threshold", "drug_class", "true_OTS"])["correct"]
                .agg(recovered="sum", total="size").reset_index())
    missed = m[m["planted"] & ~m["correct"]][["estimator", "threshold", "cell_type", "drug_class",
                                               "gene", "true_OTS", "est_OTS"]]

    # planted rows not analysed at all (e.g. a group dropped for missing visits)
    analysed = ots[["gene", "cell_type", "drug_class"]].drop_duplicates()
    planted_truth = truth[truth["expressed"] & (truth["true_OTS"] != "flat_flat")]
    not_analysed = planted_truth.merge(analysed, how="left", indicator=True)
    not_analysed = not_analysed[not_analysed["_merge"] == "left_only"][["cell_type", "drug_class", "gene", "true_OTS"]]

    # comparator queries vs truth (LMM + fixed threshold)
    qfile = res / "queries.csv"
    qrows = []
    if qfile.exists() and qfile.stat().st_size > 0:
        q = pd.read_csv(qfile)
        t = truth[truth["expressed"]].pivot_table(index=["cell_type", "gene"], columns="drug_class",
                                                   values="true_OTS", aggfunc="first")
        for label, cond in QUERY_TRUTH.items():
            want = t.copy()
            for grp, shape in cond.items():
                want = want[want[grp] == shape]
            want_set = set(want.index)
            got = q[q["query"] == label]
            got_set = set(zip(got["cell_type"], got["gene"]))
            tp = len(want_set & got_set)
            qrows.append({"query": label, "truth_matches": len(want_set), "returned": len(got_set),
                          "true_positives": tp,
                          "precision": round(tp / len(got_set), 3) if got_set else None,
                          "recall": round(tp / len(want_set), 3) if want_set else None})
    qsum = pd.DataFrame(qrows)

    summary.to_csv(res / "score_summary.csv", index=False)
    by_shape.to_csv(res / "score_by_shape.csv", index=False)
    missed.to_csv(res / "missed_planted.csv", index=False)
    qsum.to_csv(res / "score_queries.csv", index=False)
    pd.set_option("display.width", 200)
    print(summary.to_string(index=False))
    print("\nPlanted rows not analysed:", len(not_analysed))
    if len(not_analysed):
        print(not_analysed.to_string(index=False))
    print("\nBy shape (lmm, fixed):")
    print(by_shape.query("estimator == 'lmm' and threshold == 'fixed'").to_string(index=False))
    print("\nQueries (lmm, fixed):")
    print(qsum.to_string(index=False))


if __name__ == "__main__":
    main()
