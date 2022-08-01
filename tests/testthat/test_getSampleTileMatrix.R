test_that("getSampleTileMatrices works on a 1 sample test dataset", {
  
  # Check for existence of ArchR test data first:
  ArchR::addArchRThreads(threads = 10)
  withr::local_options(timeout = 600) # 10 minute timeout for download
  ArchR::addArchRVerbose(verbose = FALSE)
  
  # This only downloads the test project the first time
  # subsequent runs load the object itself
  testProj <- ArchR::getTestProject()
  
  cellPopulations <- c("C1", "C2")
  
  capture.output(
    tileResults <- scMACS::callOpenTiles(
      testProj,
      cellPopLabel = "Clusters",
      cellPopulations = cellPopulations,
      numCores = 10
    )
  )

  capture.output(
    expect_snapshot_value(
      scMACS::getSampleTileMatrix( 
        tileResults,
        cellPopulations = cellPopulations,
        threshold = 0
      ), 
      style = "json2"
    )
  )
  
})