test_that("FindCoAccessibleLinks works on a 1 sample test dataset", {

  # Check for existence of ArchR test data first:
  capture.output(ArchR::addArchRThreads(threads = 10), type = "message")
  withr::local_options(timeout = 600) # 10 minute timeout for download
  capture.output(ArchR::addArchRVerbose(verbose = FALSE), type = "message")

  # This only downloads the test project the first time
  # subsequent runs load the object itself
  capture.output(testProj <- ArchR::getTestProject(), type = "message")

  cellPopulations <- c("C1", "C2")

  TxDb <- TxDb.Hsapiens.UCSC.hg38.refGene
  Org <- org.Hs.eg.db
  capture.output(
    tileResults <- MOCHA::callOpenTiles(
      testProj,
      TxDb = TxDb,
      Org = Org,
      cellPopLabel = "Clusters",
      cellPopulations = cellPopulations,
      numCores = 1
    ),
    type = "message"
  )

  capture.output(
    SampleTileMatrix <- MOCHA::getSampleTileMatrix(
      tileResults,
      cellPopulations = cellPopulations,
      threshold = 0
    )
  )

  cellPopulation <- "C1"
  regions <- MOCHA::StringsToGRanges(c(
    "chr1:102368000-102368499",
    "chr1:101873000-101873499"
  ))
  links <- MOCHA::getCoAccessibleLinks(SampleTileMatrix,
    cellPopulation,
    regions,
    verbose = FALSE
  )

  expect_false(any((links$Peak1 == links$Peak2)))
  expect_snapshot_output(
    links,
    variant = "tiles"
  ) # change variant to "basic" to compare with previous co_accessibility.R
  expect_snapshot_output(
    nrow(links),
    variant = "nrows"
  )
})
