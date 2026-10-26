# Synthetic PRISM-like data

> **These data are simulated.** No value here was measured. Use them to build and test the pipeline, never to draw biological conclusions.

The dataset copies the design of the PRISM PBMC scRNA-seq study (Park et al., *Exp Mol Med* 2024) and plants known effects, so every analysis step can be checked against a known answer before the real data arrive.

## Generate

```bash
python scripts/simulate_prism_like.py --out data/synthetic --scenario drug           # ~20 s
python scripts/simulate_prism_like.py --out data/synthetic --scenario stress_batch   # ~20 s
python scripts/check_synthetic.py --data data/synthetic --figdir docs/figures       # UMAP + quick check
```

Background gene names come from `background_genes.txt` (regenerate with `Rscript scripts/export_background_genes.R`, which needs msigdbr 7.5.1). The `.h5ad` files (~73 MB each) are not committed; they are fully reproducible from the seed.

## Design

| | Synthetic | PRISM (published) |
| --- | --- | --- |
| Patients | 8: mepolizumab P1–P2, reslizumab P3–P4, dupilumab P5–P8 | same |
| Visits | 0, 1, 6 months (24 samples) | same |
| Cells | ~88,800 (~3,700 per sample) | 89,728 |
| Genes | 1,625 real human symbols: 125 immune, housekeeping and stress genes carrying the planted biology + 1,500 background genes | 25,242 |
| Counts | negative binomial (θ = 6), patient effects (SD 0.15 log), cell size factors, per-sample processing stress | — |

14 PBMC cell types, plus states: CD4 T cells S0/S1 and classical monocytes IL1B+/S100A+. Cell metadata include `processing_batch` and the simulated `stress_level`; real data will need the equivalent (processing date, time-to-freeze, sequencing batch).

## Two scenarios

| Scenario | Late fall of NF-κB/AP-1 genes at 6 months is… | Processing stress | Files |
| --- | --- | --- | --- |
| `drug` (default) | a **drug effect** (biology) | same level at every visit (0.5 ± 0.1) | `prism_like_synthetic.h5ad`, `ground_truth_*.csv` |
| `stress_batch` | a **processing artefact**: 0- and 1-month samples in a high-stress batch, 6-month samples in a low-stress batch; no NF-κB drug effect | 1.0 at months 0 and 1, 0.0 at month 6 (± 0.1) | `*_stress_batch.*` |

Stress raises immediate-early genes (*FOS, FOSB, JUN, JUNB, IER2, EGR1, DUSP1, ZFP36, NFKBIA, CD69*) **and** heat-shock genes (*HSPA1A, HSPA1B, DNAJB1, HSPH1, HSP90AA1, HSPA8*; plus *ATF3, IER3*) in every cell type, by 1.2 log units at stress level 1 (gene sets after van den Brink et al., *Nat Methods* 2017). The stress-batch scenario exists because the published 6-month signal is made of exactly these immediate-early genes: the pipeline must be able to tell a drug effect from a handling artefact. Heat-shock genes are the internal control.

## Planted effects

| Effect | Where | Drug class | True shape | Scenario |
| --- | --- | --- | --- | --- |
| NF-κB/AP-1 genes fall late | all cell types | both | `flat_down` | `drug` only |
| CD4 T cells switch S0 → S1 (GIMAP up, some AP-1 down) | CD4 naive, TCM, EM | both | `flat_up` / `flat_down` | both |
| IL1B+ → S100A+ classical monocytes | classical monocytes | both | IL1B module `flat_down`, S100A `flat_up` | both |
| STAT6 targets fall early (*FCER2*, *IL4R*, *IGHE*, *CISH*, *SOCS1*) | B cells | anti-IL-4Rα | `down_flat` | both |
| Interferon module rises transiently | monocytes | anti-IL-4Rα | `up_down` | both |
| *CLC*, *IL5RA*, *HDC* fall; basophils fall | basophils | anti-IL-5 | `down_flat` | both |
| Cytotoxic genes dip and recover | CD8 T | anti-IL-5 | `down_up` | both |
| Processing-stress genes fall at 6 months | all cell types | both | **artefact**: biologically `flat_flat`, observed `flat_down` | `stress_batch` only |
| Everything else | — | — | `flat_flat` | both |

Effect sizes are invented. The direction of the late NF-κB fall and the monocyte switch follows Park et al.; the drug-specific effects are made up to give each trajectory shape at least one test case.

## Ground truth files

- `ground_truth_gene_topology[_stress_batch].csv`: one row per gene × cell type × drug class, with noise-free log2 means at each visit and
  - `true_OTS`: the **biological** shape (stress held constant), threshold 0.35 on the log2 scale;
  - `expected_observed_OTS`: the shape the data should show, including any processing artefact;
  - `stress_gene`: whether the gene is raised by processing stress;
  - `expressed`: score only rows where this is `True` (noise-free mean ≥ 0.5 UMI per cell).
- `ground_truth_composition[_stress_batch].csv`: expected cell-type fractions by drug class and visit, plus the CD4 S1 and IL1B+ CM fractions.

## Limitations

- Real PBMC data have continuous states, ambient RNA, doublets and unstructured batch effects; clusters here are far cleaner, and the pipeline uses the true cell labels rather than annotating cells itself.
- The scoring threshold of the ground truth (0.35) is one of the thresholds the pipeline uses; results at 0.25 and 0.5 are reported to show sensitivity to a mismatched choice.
- Background genes have real names (a random draw from MSigDB Hallmark and Reactome genes; see `scripts/export_background_genes.R` and `background_genes.txt`), so scGPT and PerturbGen will tokenise them. But their expression is random noise, not their real biology. Only the 125 planted genes behave biologically, and only in the ways listed above.
- In real data a biologic could, in principle, move heat-shock genes, and a cold-chain artefact could move immediate-early genes without heat-shock genes; the stress check is a screen, not proof.
- There is no ground truth for in silico perturbation (Aim 3) beyond the direction of planted drug effects.
