test_that("getBackGroundObj works on a 3-sample dataset", {
  
  cellPopulations = c("C2", "C3")
  capture.output(
    STObj <- MOCHA::getSampleTileMatrix(
      MOCHA:::testTileResultsMultisample,
      cellPopulations = cellPopulations,
      threshold = 0
    )
  )

  obj <- MOCHA::getBackGroundObj(STObj)
  
  expect_snapshot_output(
    obj,
    variant = "3sample"
  )
  
  }
)