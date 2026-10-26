# Count-model pseudobulk: edgeR TMM normalisation + limma-voom with patient blocking.
#
# Input : summed-count pseudobulk from scripts/pseudobulk.py --mode sum
# Output: (1) voom logCPM per cell type, same layout as the input folder, to be
#             passed to scripts/run_longomicstraj.R as its visit-level input;
#         (2) transitions.csv: for every gene x cell type x drug class, the
#             log2 fold change and BH-adjusted p-value for M0->M1 and M1->M6,
#             from limma with duplicateCorrelation(block = patient).
#             These give each LongOmicsTraj transition statistical support.
#
# Groups (drug classes) that lack any visit in a cell type are dropped for that
# cell type, as in run_longomicstraj.R.
#
# Usage: Rscript scripts/normalise_pseudobulk.R results/synthetic/pseudobulk_counts results/synthetic/pseudobulk_voom
suppressMessages({
  library(edgeR)
  library(limma)
})

args    <- commandArgs(trailingOnly = TRUE)
in_dir  <- if (length(args) >= 1) args[1] else "results/synthetic/pseudobulk_counts"
out_dir <- if (length(args) >= 2) args[2] else "results/synthetic/pseudobulk_voom"
visits  <- c("M0", "M1", "M6")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

samples <- read.csv(file.path(in_dir, "samples.csv"))
index   <- read.csv(file.path(in_dir, "index.csv"))
file.copy(file.path(in_dir, "samples.csv"), file.path(out_dir, "samples.csv"), overwrite = TRUE)

safe <- function(x) gsub("[^A-Za-z0-9]", "", x)  # syntactic level names for contrasts
out_index <- list()
trans <- list()
for (i in seq_len(nrow(index))) {
  ct <- index$cell_type[i]
  cnt <- as.matrix(read.csv(file.path(in_dir, index$file[i]), row.names = 1, check.names = FALSE))
  s <- samples[match(colnames(cnt), samples$sample), ]
  ok <- names(which(tapply(s$visit, s$drug_class, function(x) all(visits %in% x))))
  keep <- s$drug_class %in% ok
  if (!any(keep)) next
  cnt <- cnt[, keep, drop = FALSE]; s <- s[keep, ]
  grp <- factor(paste(safe(s$drug_class), s$visit, sep = "_"))
  design <- model.matrix(~ 0 + grp); colnames(design) <- levels(grp)

  y <- DGEList(cnt, group = grp)
  k <- filterByExpr(y, group = grp)
  y <- normLibSizes(y[k, , keep.lib.sizes = FALSE], method = "TMM")
  v <- voom(y, design)
  dc <- duplicateCorrelation(v, design, block = s$patient)
  v <- voom(y, design, block = s$patient, correlation = dc$consensus.correlation)
  dc <- duplicateCorrelation(v, design, block = s$patient)
  fit <- lmFit(v, design, block = s$patient, correlation = dc$consensus.correlation)

  for (g in ok) {
    gs <- safe(g)
    cm <- makeContrasts(contrasts = c(sprintf("%s_M1-%s_M0", gs, gs), sprintf("%s_M6-%s_M1", gs, gs)),
                        levels = design)
    f2 <- eBayes(contrasts.fit(fit, cm))
    for (j in 1:2) {
      tt <- topTable(f2, coef = j, number = Inf, sort.by = "none")
      trans[[length(trans) + 1]] <- data.frame(
        cell_type = ct, drug_class = g, transition = c("M0->M1", "M1->M6")[j],
        gene = rownames(tt), logFC = tt$logFC, adj_p = tt$adj.P.Val)
    }
  }
  f <- file.path(out_dir, index$file[i])
  write.csv(v$E, f)
  out_index[[length(out_index) + 1]] <- data.frame(
    cell_type = ct, file = index$file[i], n_samples = ncol(cnt), n_genes = nrow(v$E),
    patient_correlation = round(dc$consensus.correlation, 3), groups = paste(ok, collapse = ";"))
  message(sprintf("%s: %d genes kept, patient correlation %.2f", ct, nrow(v$E), dc$consensus.correlation))
}
write.csv(do.call(rbind, out_index), file.path(out_dir, "index.csv"), row.names = FALSE)
write.csv(do.call(rbind, trans), file.path(out_dir, "transitions.csv"), row.names = FALSE)
writeLines(sprintf("edgeR %s, limma %s, R %s", packageVersion("edgeR"), packageVersion("limma"), getRversion()),
           file.path(out_dir, "session.txt"))
