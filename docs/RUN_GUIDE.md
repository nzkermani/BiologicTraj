# Running BiologicTraj on your own data

This guide is for a team that holds the data and runs the analysis on its own computers (the "code goes to the data" option). Individual cells never leave your site; you send back only the summary files listed in step 6.

Everything runs on a laptop-class CPU. Steps 1–5 take well under an hour for a PRISM-sized dataset (about 90,000 cells, 24 samples), apart from the optional mixed-model estimator (about 30 minutes) and Aim 1 with scGPT (about 30 minutes per checkpoint on all genes).

## 1. Get the code and software

```bash
git clone https://github.com/nzkermani/BiologicTraj.git
cd BiologicTraj

# Python 3.11
python3.11 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# R (4.5 or later): edgeR, limma, BiocParallel, LongOmicsTraj
Rscript scripts/install_R_packages.R      # also writes R_sessionInfo.txt
```

Optional, for Aim 1 (scGPT): follow the install order at the top of `requirements-aim1.txt`, then download the checkpoints as in the README (Environment).

Before running on real data, check everything works on the simulated data:

```bash
python scripts/simulate_prism_like.py --out data/synthetic --scenario drug
scripts/run_pipeline.sh data/synthetic/prism_like_synthetic.h5ad results/synthetic/drug \
  data/synthetic/ground_truth_gene_topology.csv empirical
```

The score file `results/synthetic/drug/longomicstraj/score_summary.csv` should match the README results table.

## 2. Prepare one .h5ad file

One AnnData file with all cells from all samples:

- `adata.X`: **raw UMI counts** (integers), not normalised values.
- `adata.var_names`: gene symbols (`FOS`, `CD3E`, ...), unique.
- `adata.obs`, one row per cell:

| Column | Required | Example | Notes |
| --- | --- | --- | --- |
| `sample` | yes | `P1_M0` | one per patient × visit |
| `patient` | yes | `P1` | pseudonymised |
| `drug_class` | yes | `anti-IL-5`, `anti-IL-4Ra` | same for all of a patient's samples |
| `visit` | yes | `M0`, `M1`, `M6` | |
| `visit_month` | yes | `0`, `1`, `6` | numeric, matching `visit` |
| `cell_type` | yes | `CD4 T naive` | your existing labels are fine to start |
| `drug` | recommended | `dupilumab` | |
| `processing_batch` | recommended | `B2` | |
| `processing_date` | recommended | `2022-03-14` | |
| `hours_to_freeze` | recommended | `3.5` | time from blood draw to freezing |
| `sequencing_batch` | recommended | `run04` | |

The processing columns are what lets the stress check tell a handling artefact from a drug effect, so please include them if you have them.

If the sample-level information is in a separate table, it can be joined onto the cells like this:

```python
import anndata as ad, pandas as pd
adata = ad.read_h5ad("your_cells.h5ad")
samples = pd.read_csv("samples.csv")          # one row per sample, columns as above
adata.obs = adata.obs.join(samples.set_index("sample"), on="sample")
adata.write_h5ad("data/prism/prism.h5ad")
```

## 3. Check the input

```bash
python scripts/check_input.py --h5ad data/prism/prism.h5ad --out results/prism/input_check
```

It prints `READY` or `NOT READY` with the reasons, and writes `input_check.txt` and a table of cells per cell type and sample. Fix every error before going on. Warnings are worth reading: a cell type with too few cells in some samples is dropped there, which is expected for depleted populations such as basophils under anti-IL-5.

## 4. Run the pipeline

```bash
scripts/run_pipeline.sh data/prism/prism.h5ad results/prism "" empirical
```

The empty `""` is the ground-truth file, which real data do not have, so the scoring step is skipped. Use `empirical,lmm` instead of `empirical` to add the mixed-model estimator.

## 5. Optional: Aim 1 (scGPT vs PCA)

```bash
python scripts/embed_and_evaluate.py --h5ad data/prism/prism.h5ad --out results/prism/aim1 \
  --n-cells 5000 --model scgpt:models/scGPT_blood --model scgpt:models/scGPT_human
```

## 6. What to send back

Send only these summary files, as agreed between the two teams. None of them holds individual cells.

| File or folder | What it holds |
| --- | --- |
| `results/prism/input_check/input_check.txt`, `cells_per_celltype_sample.csv` | cell counts per sample and cell type |
| `results/prism/stress_check/` | stress scores per sample and cell type, flags, figure |
| `results/prism/pseudobulk_voom/transitions.csv`, `index.csv`, `session.txt` | visit-to-visit test results per gene and cell type |
| `results/prism/longomicstraj/` (except per-gene files you have agreed not to share) | trajectory shape calls and queries |
| `results/prism/aim1/aim1_metrics.csv` | Aim 1 scores |
| `R_sessionInfo.txt` | software versions |

**Do not send** the `.h5ad` file, `results/prism/aim1/embeddings.npz` or `subsample_obs.csv` (cell-level). Per-sample summed counts (`pseudobulk_counts/`, `pseudobulk_voom/*.csv`) are summaries but still per person; send them only if your agreement allows.

## Problems

Open an issue on GitHub (without any data in it) or contact Dr Nazanin Zounemat-Kermani, National Heart and Lung Institute, Imperial College London. Fixes go into the GitHub repository, so both teams always run the same version: run `git pull` before each new run.
