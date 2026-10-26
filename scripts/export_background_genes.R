# Export a reproducible random set of real human gene symbols for the
# synthetic dataset's background genes.
#
# Source: msigdbr 7.5.1 (bundles MSigDB gene sets offline). Genes are drawn from
# the Hallmark and Reactome canonical-pathway sets, i.e. well-annotated
# protein-coding genes that single-cell foundation models will recognise.
#
# Usage: Rscript scripts/export_background_genes.R data/synthetic/background_genes.txt 3000
suppressMessages(library(msigdbr))
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[1] else "data/synthetic/background_genes.txt"
n   <- if (length(args) >= 2) as.integer(args[2]) else 3000L

h <- msigdbr(species = "Homo sapiens", category = "H")
r <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
pool <- unique(c(h$gene_symbol, r$gene_symbol))
# plain gene symbols only: no mitochondrial, ribosomal-protein, histone or orf-style names
pool <- pool[!grepl("^(MT-|RPL|RPS|MRP[LS]|H[1-4]|H2A|H2B)|orf|-AS", pool)]
set.seed(20261026)
writeLines(sort(sample(pool, min(n, length(pool)))), out)
cat(sprintf("pool: %d symbols; wrote %d to %s (msigdbr %s)\n",
            length(pool), min(n, length(pool)), out, as.character(packageVersion("msigdbr"))))
