if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
library(BiocManager)

if (!require("TxDb.Hsapiens.UCSC.hg38.refGene", quietly = TRUE)) {
  BiocManager::install("TxDb.Hsapiens.UCSC.hg38.refGene")
}

if (!require("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}

library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(org.Hs.eg.db)

if (requireNamespace("ArchR", quietly = TRUE) & dir.exists("PBMCSmall")) {
  library(ArchR)
  # ArchR is installed and the test project is downloaded, 
  # proceed with optional tests
  withr::local_options(list(timeout = 600))
  capture.output(proj <- ArchR::getTestProject(), type = "message")
  outdir <- dirname(ArchR::getOutputDirectory(proj))
  
  # Remove this directory after all tests
  withr::defer(unlink(ArchR::getOutputDirectory(proj), recursive = TRUE), teardown_env())
  # Remove additional empty directories
  withr::defer(unlink(file.path(outdir, "__MACOSX"), recursive = TRUE), teardown_env())
  withr::defer(unlink(file.path(outdir, "ArchRLogs"), recursive = TRUE), teardown_env())
} else { 
  message("Skipping tests involving ArchR project")
}

