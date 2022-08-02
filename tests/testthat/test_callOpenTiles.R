test_that("We can call peaks by sample", {
  
  # Check for existence of ArchR test data first:
  capture.output(ArchR::addArchRThreads(threads = 10), type = "message")
  withr::local_options(timeout = 600) # 10 minute timeout for download
  capture.output(ArchR::addArchRVerbose(verbose = FALSE), type = "message")
  
  # This only downloads the test project the first time
  # subsequent runs load the object itself
  capture.output(testProj <- ArchR::getTestProject(), type = "message")

  capture.output(expect_snapshot_value(
    scMACS::callOpenTiles(
      testProj,
      cellPopLabel = "Clusters",
      cellPopulations = c("C1", "C2"),
      numCores = 10
    ),
    style = "json2"
  ), type = "message")
  
})