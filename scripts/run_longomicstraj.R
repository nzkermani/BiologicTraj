# Run LongOmicsTraj on pseudobulk data, one object per cell type.
#
# Features = genes; group = drug class; subject = patient; visits = M0, M1, M6.
# Estimators: empirical visit means and a linear mixed model
#   y ~ group * visit + (1 | subject_id)
# Thresholds: "fixed" at delta_min (default 0.35, the synthetic ground-truth
# threshold) and the package default "hybrid".
#
# Usage: Rscript scripts/run_longomicstraj.R results/synthetic/pseudobulk results/synthetic/longomicstraj [delta_min] [n_cores]
suppressMessages({
  library(LongOmicsTraj)
  library(BiocParallel)
})

args      <- commandArgs(trailingOnly = TRUE)
pb_dir    <- if (length(args) >= 1) args[1] else "results/synthetic/pseudobulk"
out_dir   <- if (length(args) >= 2) args[2] else "results/synthetic/longomicstraj"
delta_min <- if (length(args) >= 3) as.numeric(args[3]) else 0.35
n_cores   <- if (length(args) >= 4) as.integer(args[4]) else max(1L, parallel::detectCores() - 1L)
min_expr  <- log2(0.1)  # analyse genes reaching 0.1 UMI per cell in at least one sample
visits    <- c("M0", "M1", "M6")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

samples <- read.csv(file.path(pb_dir, "samples.csv"))
index   <- read.csv(file.path(pb_dir, "index.csv"))
bp      <- if (n_cores > 1) MulticoreParam(workers = n_cores) else SerialParam()

run_one <- function(lot, estimator, threshold) {
  lot <- if (estimator == "empirical") {
    lot_from_empirical(lot, assay = "transcriptomics", visits = visits, overwrite = TRUE, verbose = FALSE)
  } else {
    lot_from_lmm(lot, assay = "transcriptomics", formula_fixed = y ~ group * visit + (1 | subject_id),
                 visits = visits, BPPARAM = bp, overwrite = TRUE, verbose = FALSE)
  }
  lot <- lot_compute_deltas(lot, assay = "transcriptomics", visits = visits, overwrite = TRUE)
  lot <- if (threshold == "fixed") {
    lot_compute_thresholds(lot, assay = "transcriptomics", method = "fixed", delta_min = delta_min,
                           visits = visits, overwrite = TRUE, verbose = FALSE)
  } else {
    lot_compute_thresholds(lot, assay = "transcriptomics", method = "hybrid",
                           visits = visits, overwrite = TRUE, verbose = FALSE)
  }
  lot <- lot_map_to_ots(lot, visits = visits, overwrite = TRUE)
  list(lot = lot, ots = as.data.frame(lot_ots(lot)))
}

all_ots <- list()
queries <- list()
timing  <- list()
for (i in seq_len(nrow(index))) {
  ct <- index$cell_type[i]
  m  <- as.matrix(read.csv(file.path(pb_dir, index$file[i]), row.names = 1, check.names = FALSE))
  m  <- m[apply(m, 1, max) > min_expr, , drop = FALSE]
  s  <- samples[match(colnames(m), samples$sample), ]
  meta <- data.frame(sample_id = colnames(m), subject_id = s$patient, visit = s$visit,
                     group = s$drug_class, row.names = colnames(m))
  # every group needs every visit for a complete visit mean matrix
  ok_groups <- names(which(tapply(meta$visit, meta$group, function(x) all(visits %in% x))))
  dropped <- setdiff(unique(meta$group), ok_groups)
  if (length(dropped)) {
    message(ct, ": dropping group(s) without all visits: ", paste(dropped, collapse = ", "))
    keep <- meta$group %in% ok_groups
    m <- m[, keep, drop = FALSE]; meta <- meta[keep, , drop = FALSE]
  }
  if (!length(ok_groups)) next
  base <- lot_create(expression = m, metadata = meta, visit_levels = visits,
                     time_values = c(M0 = 0, M1 = 1, M6 = 6))
  for (est in c("empirical", "lmm")) {
    for (thr in c("fixed", "hybrid")) {
      t0 <- Sys.time()
      res <- tryCatch(run_one(base, est, thr), error = function(e) { message(ct, " ", est, " ", thr, ": ", conditionMessage(e)); NULL })
      if (is.null(res)) next
      o <- res$ots
      o$cell_type <- ct; o$threshold <- thr
      all_ots[[length(all_ots) + 1]] <- o
      timing[[length(timing) + 1]] <- data.frame(cell_type = ct, estimator = est, threshold = thr,
                                                 n_genes = nrow(m), seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
      # comparator-aware queries (main analysis: LMM + fixed threshold)
      if (est == "lmm" && thr == "fixed" && length(ok_groups) == 2) {
        q <- function(cond, label) {
          f <- tryCatch(lot_query(res$lot, conditions = cond), error = function(e) character(0))
          ids <- if (is.data.frame(f)) unique(f$object_id) else unique(as.character(f))
          if (length(ids)) data.frame(cell_type = ct, query = label, gene = ids) else NULL
        }
        queries[[length(queries) + 1]] <- rbind(
          q(list(`anti-IL-5` = "flat_down", `anti-IL-4Ra` = "flat_down"), "shared delayed suppression"),
          q(list(`anti-IL-4Ra` = "down_flat", `anti-IL-5` = NOT("down_flat")), "early sustained fall, anti-IL-4Ra only"),
          q(list(`anti-IL-4Ra` = "up_down", `anti-IL-5` = "flat_flat"), "transient rise, anti-IL-4Ra only"),
          q(list(`anti-IL-5` = "down_up", `anti-IL-4Ra` = "flat_flat"), "dip and recovery, anti-IL-5 only")
        )
      }
    }
  }
  message(sprintf("%s: %d genes done", ct, nrow(m)))
}

ots <- do.call(rbind, all_ots)
write.csv(ots[, c("cell_type", "object_id", "group", "estimator", "threshold", "topology_label")],
          file.path(out_dir, "ots_all.csv"), row.names = FALSE)
write.csv(do.call(rbind, queries), file.path(out_dir, "queries.csv"), row.names = FALSE)
write.csv(do.call(rbind, timing), file.path(out_dir, "timing.csv"), row.names = FALSE)
writeLines(c(sprintf("LongOmicsTraj %s, R %s", packageVersion("LongOmicsTraj"), getRversion()),
             sprintf("delta_min = %s; min_expr = log2(0.1); cores = %d", delta_min, n_cores)),
           file.path(out_dir, "session.txt"))
cat("rows:", nrow(ots), "\n")
