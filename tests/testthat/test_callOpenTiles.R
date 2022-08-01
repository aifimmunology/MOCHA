test_that("We can call peaks by sample", {
  
  # Check for existence of ArchR test data first:
  ArchR::addArchRThreads(threads = 10)
  withr::local_options(timeout = 600) # 10 minute timeout for download
  ArchR::addArchRVerbose(verbose = FALSE)
  
  # This only downloads the test project the first time
  # subsequent runs load the object itself
  testProj <- ArchR::getTestProject()

  capture.output(expect_snapshot_value(
    scMACS::callOpenTiles(
      testProj,
      cellPopLabel = "Clusters",
      cellPopulations = c("C1", "C2"),
      numCores = 10
    ),
    style = "json2"
  ))
  
})