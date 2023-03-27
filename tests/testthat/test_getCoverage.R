if (requireNamespace("TxDb.Hsapiens.UCSC.hg38.refGene", quietly = TRUE)) {
  
  test_that("getCoverage works on a 1 sample test dataset", {
    capture.output(
      covFiles <- getCoverage(
        popFrags = MOCHA::exampleFragments,
        normFactor = c(1,1),
        filterEmpty = FALSE,
        cl = 1, TxDb = TxDb.Hsapiens.UCSC.hg38.refGene
      )
    )
    expect_snapshot(covFiles)
  })
  
  test_that("getCoverage works on a 1 sample test dataset with filterEmpty=TRUE", {
    capture.output(
      covFiles <- getCoverage(
        popFrags = MOCHA::exampleFragments,
        normFactor = c(1,1),
        filterEmpty = TRUE,
        cl = 1, TxDb = TxDb.Hsapiens.UCSC.hg38.refGene
      )
    )
    expect_snapshot(covFiles)
  })
  
  test_that("getSpecificCoverage works on a 1 sample test dataset", {
    capture.output(
      covFiles <- getCoverage(
        popFrags = MOCHA::exampleFragments,
        normFactor = c(1,1),
        filterEmpty = TRUE,
        cl = 1, TxDb = TxDb.Hsapiens.UCSC.hg38.refGene
      )
    )
    regions <- StringsToGRanges(c("chr1:565291-569412", "chr2:243031253-243034339"))
    seqinfo(regions) <- seqinfo(covFiles[[1]])
    counts <- getSpecificCoverage(
      covFiles, 
      regions, 
      numCores = 1
    ) 
    expect_snapshot(counts)
  })
    
}