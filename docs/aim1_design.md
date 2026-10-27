# Aim 1 (revised): representation, with the patient as the unit

**Status:** design proposal, not yet implemented. The current script (`scripts/embed_and_evaluate.py`) computes the cell-level checks in section 4; sections 2–3 are to be built.

## 1. Why revise

The first version of Aim 1 scored embeddings with cell-level k-NN: can a held-out patient's cells be classified as month 6 vs earlier, or by drug class? That is a useful sanity check but a weak main result:

- **Cells are not independent.** ~3,700 cells per sample, 24 samples, 8 patients. Cell-level accuracies hide that the effective sample size is 8.
- **Visit is confounded with processing.** In real PRISM data, a month-6 signal can come from handling or sequencing batch (see the `stress_batch` scenario). A visit classifier rewards that just as much as a drug effect.
- **"Is this cell from month 6?" is not a question anyone needs answered.** The questions are what changes, in which cells, consistently across patients, and whether it differs by drug.

## 2. Questions

| # | Question | Unit | Where it lives |
| --- | --- | --- | --- |
| Q1 | Does the embedding keep cell identity and remove technical batch? | cell, within sample | Aim 1 (this doc) |
| Q2 | Within each cell type, is there a shift from M0 to M6 that points the same way across patients? | patient × cell type | Aim 1 (this doc) |
| Q3 | Does that shift differ between anti-IL-5 and anti-IL-4Rα? | patient | Aim 1, exploratory |
| Q4 | Which cell types and states change in abundance? | sample | composition analysis (separate; see README Next steps) |

**Not answerable in PRISM:** responder vs non-responder. All 8 patients are good clinical responders. A response question needs a second cohort.

## 3. Methods

### Embeddings compared

| Method | Notes |
| --- | --- |
| PCA (50 PCs on HVGs) | baseline, already in the script |
| scVI | probabilistic baseline; `processing_batch` as batch covariate |
| scPoli | gives a per-sample embedding as well as per-cell; useful for Q2–Q3 |
| scGPT zero-shot, blood and whole-human checkpoints | already in the script |
| scGPT fine-tuned | HPC; only if zero-shot is competitive |
| PerturbGen | Aims 2–3; its embedding can be scored here too |

All methods get the same cells (stratified subsample, fixed seed) and the same gene set. A **gene-panel sensitivity run** repeats everything on the 125 planted genes only: on synthetic data the 1,500 background genes have real names but random expression, which handicaps models that read gene names (scGPT) more than PCA.

### Q1: identity vs batch (scIB-style)

Follow the scIB benchmark (Luecken et al., *Nat Methods* 2022): report a **biological conservation** score and a **batch removal** score side by side, not one number.

- Bio conservation: leave-one-patient-out k-NN cell-type accuracy (existing), ARI/NMI of Leiden clusters vs labels, cell-type silhouette.
- Batch removal: computed on **`processing_batch`**, within cell type (e.g. batch silhouette, iLISI).
- **Patient is not treated as batch.** Patient differences are partly biology (Aim 1 asks to preserve them). The same-patient neighbour fraction (existing) is reported descriptively, next to its chance level.

### Q2: consistent treatment shift (the main result)

1. For each method, average cell embeddings per **patient × cell type × visit** (pseudobulk in embedding space; cell types with ≥ 20 cells per sample).
2. Shift vector per patient and cell type: Δ = mean(M6) − mean(M0); likewise M1 − M0.
3. **Consistency:** mean pairwise cosine similarity of Δ across the 8 patients. High = the change points the same way in different people.
4. **Null:** permute visit labels within each patient (all 3! = 6 orderings per patient, so an exact or large Monte Carlo null), recompute. Report the observed value as a percentile of the null.
5. **Magnitude**, standardised by within-sample spread, so methods with different scales are comparable.

On synthetic data the expected answer is known: cell types with planted late effects should be consistent; untouched cell types should sit inside the null.

### Q3: drug class (exploratory)

- Leave-one-patient-out: classify each held-out patient's Δ (concatenated across cell types) as anti-IL-5 or anti-IL-4Rα, from the other 7 patients.
- With 4 vs 4 patients the minimum non-zero error is 1/8. Report with a label-permutation null (70 balanced relabelings) and treat as hypothesis-generating.

### Control: processing artefact

Run Q2 on both synthetic scenarios.

- In `drug`, the M6 shift is real.
- In `stress_batch`, it is an artefact and the embedding shift should **align with the stress score**.

For real data, report the correlation between each patient's Δ and their change in heat-shock/immediate-early score (from `stress_check.py`) and in processing metadata. A method that scores well on Q2 only because it encodes stress is not a better method.

## 4. Kept as sanity checks

From the current script, reported but not headline:

- leave-one-patient-out k-NN for visit (M6 vs M0/M1), balanced accuracy, per cell type;
- drug-class k-NN in B cells at M6;
- same-patient neighbour fraction vs chance.

## 5. What synthetic results can and cannot show

- **Can:** whether each metric recovers planted effects, rejects untouched cell types, and flags the stress artefact; whether code runs end to end.
- **Cannot:** rank scGPT against PCA for real biology. Background-gene expression is noise, so scGPT's pretrained gene knowledge has nothing true to match. Ranking needs a real public longitudinal PBMC dataset, then PRISM.

## 6. Build order

1. Add `processing_batch` batch-removal metrics and ARI/NMI to `embed_and_evaluate.py`.
2. New script `aim1_shift_consistency.py` for Q2–Q3 (reads `embeddings.npz` + `subsample_obs.csv`), including the stress control.
3. Add scVI and scPoli (CPU is fine at 5,000 cells; check timings).
4. Gene-panel sensitivity run (125 planted genes).
5. Choose a public longitudinal PBMC dataset for a real-data check before PRISM.

## References

- Luecken MD, Büttner M, Chaichoompu K, et al. Benchmarking atlas-level data integration in single-cell genomics. *Nat Methods* 2022;19:41–50.
- Lopez R, Regier J, Cole MB, Jordan MI, Yosef N. Deep generative modeling for single-cell transcriptomics. *Nat Methods* 2018;15:1053–1058. (scVI)
- De Donno C, Hediyeh-Zadeh S, Moinfar AA, et al. Population-level integration of single-cell datasets enables multi-scale analysis across samples. *Nat Methods* 2023;20:1683–1692. (scPoli)
- Büttner M, Ostner J, Müller CL, Theis FJ, Schubert B. scCODA is a Bayesian model for compositional single-cell data analysis. *Nat Commun* 2021;12:6876. (for Q4)
