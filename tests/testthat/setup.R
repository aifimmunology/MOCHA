if (requireNamespace("TxDb.Hsapiens.UCSC.hg38.refGene", quietly = TRUE)) {
  library(TxDb.Hsapiens.UCSC.hg38.refGene)
}

if (requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)) {
  library(BSgenome.Hsapiens.UCSC.hg19)
}

if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  library(org.Hs.eg.db)
}

if (requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)) {
  library(BSgenome.Hsapiens.UCSC.hg19)
}

TxDb <- TxDb.Hsapiens.UCSC.hg38.refGene
Org <- org.Hs.eg.db
genome <- BSgenome.Hsapiens.UCSC.hg19

if (requireNamespace("ArchR", quietly = TRUE) && dir.exists("PBMCSmall")) {
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
}
