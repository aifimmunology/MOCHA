library(TxDb.Hsapiens.UCSC.hg38.refGene) 
library(org.Hs.eg.db)

if (!requireNamespace("ArchR", quietly = TRUE)) {
  message("Skipping tests involving downloading ArchR project")
} else {
  # ArchR is installed and the test project is downloaded, 
  # proceed with optional tests
  capture.output(proj <- ArchR::getTestProject(), type = "message")
  outdir <- dirname(getOutputDirectory(proj))
  
  # Remove this directory after all tests
  # withr::defer(unlink(getOutputDirectory(proj), recursive = TRUE), teardown_env())
  # 
  # Remove additional empty directories
  withr::defer(unlink(file.path(outdir, "__MACOSX"), recursive = TRUE), teardown_env())
  withr::defer(unlink(file.path(outdir, "ArchRLogs"), recursive = TRUE), teardown_env())
}

