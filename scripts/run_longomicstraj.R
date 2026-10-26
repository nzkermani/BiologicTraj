# Run LongOmicsTraj on pseudobulk data, one object per cell type.
#
# Features = genes; group = drug class; subject = patient; visits = M0, M1, M6.
# Estimators: empirical visit means and a linear mixed model
#   y ~ group * visit + (1 | subject_id)
# Thresholds: "fixed" at 0.25, 0.35 and 0.5 (the synthetic ground truth uses
# 0.35, so 0.25 and 0.5 test sensitivity to a mismatched threshold) and the
# package default "hybrid". Queries use the fixed threshold given as delta_min.
#
# Input: either voom logCPM (scripts/normalise_pseudobulk.R; recommended, use
# min_expr = -Inf because filterByExpr has already filtered) or log2 mean
# pseudobulk (scripts/pseudobulk.py --mode mean; min_expr = log2(0.1)).
#
# Usage: Rscript scripts/run_longomicstraj.R <input_dir> <out_dir> [delta_min] [n_cores] [min_expr]
suppressMessages({
  library(LongOmicsTraj)
  library(BiocParallel)
})

args      <- commandArgs(trailingOnly = TRUE)
pb_dir    <- if (length(args) >= 1) args[1] else "results/synthetic/pseudobulk"
out_dir   <- if (length(args) >= 2) args[2] else "results/synthetic/longomicstraj"
delta_min <- if (length(args) >= 3) as.numeric(args[3]) else 0.35
n_cores   <- if (length(args) >= 4) as.integer(args[4]) else max(1L, parallel::detectCores() - 1L)
min_expr  <- if (length(args) >= 5) as.numeric(args[5]) else log2(0.1)  # minimum max value to analyse a gene
fixed_thr <- sort(unique(c(0.25, 0.35, 0.5, delta_min)))
visits    <- c("M0", "M1", "M6")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

samples <- read.csv(file.path(pb_dir, "samples.csv"))
index   <- read.csv(file.path(pb_dir, "index.csv"))
bp      <- if (n_cores > 1) MulticoreParam(workers = n_cores) else SerialParam()

estimate <- function(lot, estimator, n_groups) {
  if (estimator == "empirical") {
    lot_from_empirical(lot, assay = "transcriptomics", visits = visits, overwrite = TRUE, verbose = FALSE)
  } else {
    # a group x visit interaction needs two groups; fall back to visit only
    f <- if (n_groups > 1) y ~ group * visit + (1 | subject_id) else y ~ visit + (1 | subject_id)
    lot_from_lmm(lot, assay = "transcriptomics", formula_fixed = f,
                 visits = visits, BPPARAM = bp, overwrite = TRUE, verbose = FALSE)
  }
}

run_one <- function(lot, threshold) {
  lot <- lot_compute_deltas(lot, assay = "transcriptomics", visits = visits, overwrite = TRUE)
  lot <- if (startsWith(threshold, "fixed")) {
    lot_compute_thresholds(lot, assay = "transcriptomics", method = "fixed",
                           delta_min = as.numeric(sub("fixed_", "", threshold)),
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
    t0 <- Sys.time()
    est_lot <- tryCatch(estimate(base, est, length(ok_groups)), error = function(e) { message(ct, " ", est, ": ", conditionMessage(e)); NULL })
    if (is.null(est_lot)) next
    for (thr in c(paste0("fixed_", fixed_thr), "hybrid")) {
      res <- tryCatch(run_one(est_lot, thr), error = function(e) { message(ct, " ", est, " ", thr, ": ", conditionMessage(e)); NULL })
      if (is.null(res)) next
      o <- res$ots
      o$cell_type <- ct; o$threshold <- thr
      all_ots[[length(all_ots) + 1]] <- o
      timing[[length(timing) + 1]] <- data.frame(cell_type = ct, estimator = est, threshold = thr,
                                                 n_genes = nrow(m), seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
      # comparator-aware queries (main analysis: LMM + fixed threshold)
      if (est == "lmm" && thr == paste0("fixed_", delta_min) && length(ok_groups) == 2) {
        q <- function(cond, label) {
          # lot_query returns an S4Vectors DataFrame of matching OTS rows
          f <- tryCatch(as.data.frame(lot_query(res$lot, conditions = cond)), error = function(e) NULL)
          ids <- if (!is.null(f) && "object_id" %in% names(f)) unique(as.character(f$object_id)) else character(0)
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
             sprintf("input = %s; query threshold = %s; fixed thresholds = %s; min_expr = %s; cores = %d",
                     pb_dir, delta_min, paste(fixed_thr, collapse = ","), min_expr, n_cores)),
           file.path(out_dir, "session.txt"))
cat("rows:", nrow(ots), "\n")
