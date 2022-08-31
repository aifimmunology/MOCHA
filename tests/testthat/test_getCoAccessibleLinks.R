test_that("FindCoAccessibleLinks works on a 1 sample test dataset", {
  
  # Check for existence of ArchR test data first:
  capture.output(ArchR::addArchRThreads(threads = 10), type = "message")
  withr::local_options(timeout = 600) # 10 minute timeout for download
  capture.output(ArchR::addArchRVerbose(verbose = FALSE), type = "message")
  
  # This only downloads the test project the first time
  # subsequent runs load the object itself
  capture.output(testProj <- ArchR::getTestProject(), type = "message")
  
  cellPopulations <- c("C1", "C2")
  library(TxDb.Hsapiens.UCSC.hg38.refGene)
  library(org.Hs.eg.db)
  TxDb <- TxDb.Hsapiens.UCSC.hg38.refGene
  Org <- org.Hs.eg.db
  
  capture.output(
    tileResults <- scMACS::callOpenTiles( 
      testProj,
      cellPopLabel = "Clusters",
      cellPopulations = cellPopulations,
      TxDb = TxDb,
      Org = Org
    ), type = "message"
  )
  
  capture.output(
    SampleTileMatrix <- scMACS::getSampleTileMatrix( 
          tileResults,
          cellPopulations = cellPopulations,
          threshold = 0
    )
  )
  
  cellPopulation <- "C1"
  regions <- scMACS::StringsToGRanges(c(
      "chr1:102368000-102368499",
      "chr1:101873000-101873499"
  ))
  links <- scMACS::getCoAccessibleLinks(SampleTileMatrix, 
                                        cellPopulation, 
                                        regions, 
                                        verbose=TRUE)
  
  expect_false(any((links$Peak1 == links$Peak2)))
  expect_snapshot_output(
    links,
    variant = "basic"
  )
  expect_snapshot_output(
    nrow(links),
    variant = "nrows"
  )
  
  
})