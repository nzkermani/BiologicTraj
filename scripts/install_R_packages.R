# Install the R packages the pipeline needs (steps 3 and 4).
# The results in this repository were produced under R 4.5.3.
# Usage: Rscript scripts/install_R_packages.R
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("edgeR", "limma", "BiocParallel"), update = FALSE, ask = FALSE)
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
# LongOmicsTraj's DESCRIPTION asks for R >= 4.6; under R 4.5 it installs if that
# requirement is relaxed (see README, Environment).
remotes::install_github("nzkermani/LongOmicsTraj", upgrade = "never")

# Record versions next to the results (send this file back with the outputs)
writeLines(capture.output(sessionInfo()), "R_sessionInfo.txt")
cat("Installed. Versions written to R_sessionInfo.txt\n")
