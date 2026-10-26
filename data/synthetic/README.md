# Synthetic PRISM-like data

> **These data are simulated.** No value here was measured. Use them to build and test the pipeline, never to draw biological conclusions.

The dataset copies the design of the PRISM PBMC scRNA-seq study (Park et al., *Exp Mol Med* 2024) and plants known effects, so every analysis step can be checked against a known answer before the real data arrive.

## Generate

```bash
python scripts/simulate_prism_like.py --out data/synthetic --seed 42
python scripts/check_synthetic.py --data data/synthetic --figdir docs/figures
```

The background gene names come from `background_genes.txt` (regenerate with `Rscript scripts/export_background_genes.R`, which needs msigdbr 7.5.1). Generation takes about 20 seconds. `prism_like_synthetic.h5ad` (~73 MB) is not committed, because it is fully reproducible from the seed.

## Design

| | Synthetic | PRISM (published) |
| --- | --- | --- |
| Patients | 8: mepolizumab P1–P2, reslizumab P3–P4, dupilumab P5–P8 | same |
| Visits | 0, 1, 6 months (24 samples) | same |
| Cells | ~88,800 (~3,700 per sample) | 89,728 |
| Genes | 1,617 real human symbols: 117 immune/housekeeping genes carrying the planted biology + 1,500 background genes | 25,242 |
| Counts | negative binomial (θ = 6), patient effects (SD 0.15 log), cell size factors | — |

14 PBMC cell types, plus states: CD4 T cells S0/S1 and classical monocytes IL1B+/S100A+.

## Planted effects

| Effect | Where | Drug class | True shape |
| --- | --- | --- | --- |
| NF-κB/AP-1 genes fall late | all cell types | both | `flat_down` |
| CD4 T cells switch S0 → S1 (GIMAP genes up) | CD4 naive, TCM, EM | both | `flat_up` |
| IL1B+ → S100A+ classical monocytes | classical monocytes | both | IL1B module `flat_down`, S100A `flat_up` |
| STAT6 targets fall early (*FCER2*, *IL4R*, *IGHE*, *CISH*, *SOCS1*) | B cells | anti-IL-4Rα | `down_flat` |
| Interferon module rises transiently | monocytes | anti-IL-4Rα | `up_down` |
| *CLC*, *IL5RA*, *HDC* fall; basophils fall | basophils | anti-IL-5 | `down_flat` |
| Cytotoxic genes dip and recover | CD8 T | anti-IL-5 | `down_up` |
| Everything else | — | — | `flat_flat` |

Effect sizes are invented. The direction of the late NF-κB fall and the monocyte switch follows Park et al.; the drug-specific effects are made up to give each trajectory shape at least one test case.

## Ground truth files

- `ground_truth_gene_topology.csv`: one row per gene × cell type × drug class, with noise-free log2 means at each visit and the true shape (`true_OTS`). The shapes use a log2 change threshold of 0.35. Only rows with `expressed = True` (noise-free mean ≥ 0.5 UMI per cell) should be scored.
- `ground_truth_composition.csv`: expected cell-type fractions by drug class and visit, plus the CD4 S1 and IL1B+ CM fractions.

## Baseline check (seed 42)

Pseudobulk means per sample and cell type, averaged within drug class (`check_synthetic.py`):

- planted shapes recovered: 571/574 (99.5%)
- false non-flat calls among true `flat_flat` rows: 428 of 9,755 (4.4%)

This is the bar the LongOmicsTraj mixed-model estimator should meet or beat. The false calls are deliberate noise and a useful test of threshold choice.

## Limitations

- Real PBMC data have continuous states, ambient RNA, doublets and batch effects; clusters here are far cleaner.
- Background genes have real names (a random draw from MSigDB Hallmark and Reactome genes; see `scripts/export_background_genes.R` and `background_genes.txt`), so scGPT and PerturbGen will tokenise them. But their expression is random noise, not their real biology: a foundation model's prior about, say, *ABCB6* is not reflected in these data. Only the 117 planted genes behave biologically, and only in the ways listed above.
- There is no ground truth for in silico perturbation (Aim 3) beyond the direction of planted drug effects.
