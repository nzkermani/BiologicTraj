"""Aim 1: do cell embeddings capture treatment, drug and cell type, without
washing out patients? Compares scGPT zero-shot embeddings with a PCA baseline
on the same cells.

Embeddings
  pca          log-normalised expression, 2,000 highly variable genes (or all if
               fewer), 50 principal components. Always computed.
  scgpt:<dir>  scGPT zero-shot cell embeddings from a pretrained checkpoint
               folder (best_model.pt, vocab.json, args.json), via
               scgpt.tasks.embed_data on CPU. Can be given more than once.

Metrics (all on the same stratified subsample)
  cell_type_knn_acc   k-NN accuracy for cell type, leave-one-patient-out
  visit_knn_bacc      within each cell type, k-NN balanced accuracy for
                      6 months vs earlier, leave-one-patient-out; mean over cell types
  drug_knn_bacc_B     in B cells at 6 months (the only cell type where the planted
                      drug difference persists to 6 months), k-NN balanced accuracy
                      for drug class, leave-one-patient-out. With 4 patients per
                      drug, leaving one out tilts training against the held-out
                      class, so this metric is biased below 0.5 when there is no signal.
  same_patient_frac   fraction of each cell's 15 nearest neighbours (same cell
                      type) from the same patient; lower = patients better mixed.
                      Chance level is reported alongside.

Leave-one-patient-out: the neighbours of a held-out patient's cells come only
from other patients, so a method cannot score by recognising the patient.

Usage:
  python scripts/embed_and_evaluate.py --h5ad data/synthetic/prism_like_synthetic.h5ad \
     --out results/synthetic/aim1 --n-cells 5000 \
     --model scgpt:models/scGPT_blood --model scgpt:models/scGPT_human
"""

import argparse
import json
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import balanced_accuracy_score
from sklearn.neighbors import KNeighborsClassifier, NearestNeighbors

K = 15


def stratified_subsample(a, n, seed=0):
    rng = np.random.default_rng(seed)
    per = max(1, n // a.obs["sample"].nunique())
    idx = []
    for _, ix in a.obs.groupby("sample", observed=True).indices.items():
        idx.extend(rng.choice(ix, size=min(per, len(ix)), replace=False))
    return a[np.sort(idx)].copy()


def pca_embedding(a):
    b = a.copy()
    sc.pp.normalize_total(b, target_sum=1e4)
    sc.pp.log1p(b)
    n_hvg = min(2000, b.n_vars)
    sc.pp.highly_variable_genes(b, n_top_genes=n_hvg)
    sc.pp.pca(b, n_comps=50, mask_var="highly_variable")
    return b.obsm["X_pca"]


def scgpt_embedding(a, model_dir):
    from scgpt.tasks import embed_data
    b = a.copy()
    b.var["feature_name"] = b.var_names
    out = embed_data(b, model_dir=model_dir, gene_col="feature_name", max_length=1200,
                     batch_size=32, device="cpu", use_fast_transformer=False, return_new_adata=True)
    vocab = json.loads((Path(model_dir) / "vocab.json").read_text())
    in_vocab = int(np.isin(a.var_names, list(vocab)).sum())
    return np.asarray(out.X), in_vocab


def lopo_knn(X, y, groups, mask=None):
    """Leave-one-patient-out k-NN predictions."""
    pred = np.empty(len(y), dtype=object)
    idx_all = np.arange(len(y)) if mask is None else np.where(mask)[0]
    for g in np.unique(groups[idx_all]):
        te = idx_all[groups[idx_all] == g]
        tr = idx_all[groups[idx_all] != g]
        if len(np.unique(y[tr])) < 2 or len(tr) < K:
            continue
        clf = KNeighborsClassifier(n_neighbors=K).fit(X[tr], y[tr])
        pred[te] = clf.predict(X[te])
    ok = np.array([p is not None for p in pred[idx_all]])
    return idx_all[ok], pred[idx_all][ok]


def evaluate(X, obs):
    ct = obs["cell_type"].astype(str).values
    pt = obs["patient"].astype(str).values
    late = np.where(obs["visit_month"].values == 6, "M6", "M0/M1")
    drug = obs["drug_class"].astype(str).values
    res = {}
    i, p = lopo_knn(X, ct, pt)
    res["cell_type_knn_acc"] = float((p == ct[i]).mean())

    vis, drg, same, chance = [], [], [], []
    for c in np.unique(ct):
        m = ct == c
        if m.sum() < 60:
            continue
        i, p = lopo_knn(X, late, pt, m)
        if len(i) and len(np.unique(late[i])) == 2:
            vis.append(balanced_accuracy_score(late[i], p))
        m6 = m & (late == "M6")
        if c in ("B naive", "B memory") and m6.sum() >= 40:
            i, p = lopo_knn(X, drug, pt, m6)
            if len(i) and len(np.unique(drug[i])) == 2:
                drg.append(balanced_accuracy_score(drug[i], p))
        Xi, pi = X[m], pt[m]
        nn = NearestNeighbors(n_neighbors=K + 1).fit(Xi)
        _, nb = nn.kneighbors(Xi)
        same.append(float((pi[nb[:, 1:]] == pi[:, None]).mean()))
        f = pd.Series(pi).value_counts(normalize=True)
        chance.append(float((f ** 2).sum()))
    res["visit_knn_bacc"] = float(np.mean(vis))
    res["drug_knn_bacc_B"] = float(np.mean(drg)) if drg else np.nan
    res["same_patient_frac"] = float(np.mean(same))
    res["same_patient_chance"] = float(np.mean(chance))
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", default="data/synthetic/prism_like_synthetic.h5ad")
    ap.add_argument("--out", default="results/synthetic/aim1")
    ap.add_argument("--n-cells", type=int, default=5000)
    ap.add_argument("--model", action="append", default=[], help="scgpt:<checkpoint dir>; repeatable")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    a = stratified_subsample(ad.read_h5ad(args.h5ad), args.n_cells, args.seed)
    print(f"subsample: {a.n_obs} cells, {a.obs['sample'].nunique()} samples, {a.n_vars} genes")
    rows, embs = [], {}

    t = time.time()
    embs["PCA (50 PCs)"] = pca_embedding(a)
    rows.append({"method": "PCA (50 PCs)", "seconds": round(time.time() - t, 1), "genes_used": a.n_vars,
                 **evaluate(embs["PCA (50 PCs)"], a.obs)})
    for spec in args.model:
        kind, path = spec.split(":", 1)
        if kind != "scgpt":
            raise SystemExit(f"unknown model kind {kind}")
        name = f"scGPT zero-shot ({Path(path).name})"
        t = time.time()
        X, in_vocab = scgpt_embedding(a, path)
        embs[name] = X
        rows.append({"method": name, "seconds": round(time.time() - t, 1), "genes_used": in_vocab,
                     **evaluate(X, a.obs)})
        print(f"{name}: {in_vocab}/{a.n_vars} genes in vocabulary, {rows[-1]['seconds']} s")

    res = pd.DataFrame(rows)
    res.to_csv(out / "aim1_metrics.csv", index=False)
    np.savez_compressed(out / "embeddings.npz", **{k.replace(" ", "_"): v for k, v in embs.items()},
                        obs_index=a.obs_names.values)
    a.obs.to_csv(out / "subsample_obs.csv")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(len(embs), 3, figsize=(15, 4.6 * len(embs)), squeeze=False)
    for r, (name, X) in enumerate(embs.items()):
        b = ad.AnnData(X=np.asarray(X, dtype=np.float32), obs=a.obs.copy())
        sc.pp.neighbors(b, n_neighbors=K, use_rep="X")
        sc.tl.umap(b, random_state=0)
        for c, col in enumerate(["cell_type", "visit", "patient"]):
            sc.pl.umap(b, color=col, ax=axes[r, c], show=False, title=f"{name}: {col}",
                       legend_loc="right margin" if c == 0 else "on data" if col == "visit" else None)
    fig.suptitle("Aim 1 on SYNTHETIC data (not real measurements)")
    fig.tight_layout()
    fig.savefig(out / "aim1_umaps.png", dpi=110)
    pd.set_option("display.width", 200)
    print(res.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
