test_that("getDifferentialAccessibleTiles works on a 3 sample test dataset", {
  
  cellPopulations = c("C3")
  capture.output(
    SampleTileMatrix <- MOCHA::getSampleTileMatrix(
      MOCHA:::testTileResultsMultisample,
      cellPopulations = cellPopulations,
      threshold = 0
    )
  )

  cellPopulation <- "C3"
  capture.output(
      differentials <- MOCHA::getDifferentialAccessibleTiles(
      SampleTileMatrix,
      cellPopulation,
      groupColumn = "Sample",
      foreground = "scATAC_BMMC_R1",
      background = "scATAC_CD34_BMMC_R1",
      outputGRanges = TRUE,
      numCores = detectCores(),
      verbose = FALSE
    )
  )         
  expect_snapshot_output(
    differentials,
    variant = "differentials"
  )
  expect_snapshot_output(
    nrow(differentials),
    variant = "nrows"
  )
})
