test_that("We can subset a tileResults object by celltypes", {
  
  capture.output(
    obj <- subsetMOCHAObject(
      MOCHA:::testTileResultsMultisample,
      subsetBy = "celltypes",
      groupList = c("C3"),
      na.rm = TRUE,
      subsetPeaks = FALSE,
      verbose = FALSE)
  )
   
  expect_snapshot_output(
    obj,
    variant = "tileResults"
  )
})

test_that("We can subset a sampleTileMatrix object by celltypes", {
  
  capture.output(
    SampleTileMatrix <- MOCHA::getSampleTileMatrix(
      MOCHA:::testTileResultsMultisample,
      cellPopulations = "all",
      threshold = 0
    )
  )
  
  capture.output(
    obj <- subsetMOCHAObject(
      SampleTileMatrix,
      subsetBy = "celltypes",
      groupList = c("C3"),
      na.rm = TRUE,
      subsetPeaks = FALSE,
      verbose = FALSE)
  )
   
  expect_snapshot_output(
    obj,
    variant = "sampleTileMatrix"
  )
})


test_that("We can subset a sampleTileMatrix object - and peaks - by celltypes", {
  
  capture.output(
    SampleTileMatrix <- MOCHA::getSampleTileMatrix(
      MOCHA:::testTileResultsMultisample,
      cellPopulations = "All",
      threshold = 0
    )
  )
  
  capture.output(
    obj <- subsetMOCHAObject(
      SampleTileMatrix,
      subsetBy = "celltypes",
      groupList = c("C3"),
      na.rm = TRUE,
      subsetPeaks = TRUE,
      verbose = FALSE)
  )
   
  expect_snapshot_output(
    obj,
    variant = "sampleTileMatrix_peaks"
  )
})
