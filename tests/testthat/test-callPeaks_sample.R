test_that("We can call peaks by sample", {
  
  # testProj <- loadArchRProject("PBMCSmall", force = TRUE)
  withr::local_options(timeout = 600) # 10 minute timeout for download
  testProj <- getTestProject()
  
  peakCallResults <- scMACS::callPeaks_by_sample(
    testProj,
    metadf,
    cellSubsets = "C1",
    cellCol_label_name = "Clusters",
    cellType_sample_label_name = "Sample",
    returnAllPeaks = TRUE,
    numCores = 10
  )
  
})
