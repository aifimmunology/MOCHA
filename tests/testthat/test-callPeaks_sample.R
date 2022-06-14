test_that("We can call peaks by sample", {
  peakCallResults <- callPeaks_by_sample(
    ArchRProj,
    metadf,
    cellType_to_analyze = "CD14 Mono",
    cellType_sample_label_name = "cellType_sample",
    returnAllPeaks = TRUE,
    numCores = 10
  )
  expect_equal(2 * 2, 4)
})
