library(ArchR)
library(TxDb.Hsapiens.UCSC.hg38.refGene) 
library(org.Hs.eg.db)
# Download test project before all tests
capture.output(proj <- ArchR::getTestProject(), type = "message")
outdir <- dirname(getOutputDirectory(proj))

# Remove this directory after all tests
withr::defer(unlink(getOutputDirectory(proj), recursive = TRUE), teardown_env())
# Remove additional empty directories
withr::defer(unlink(file.path(outdir, "__MACOSX"), recursive = TRUE), teardown_env())
withr::defer(unlink(file.path(outdir, "ArchRLogs"), recursive = TRUE), teardown_env())
