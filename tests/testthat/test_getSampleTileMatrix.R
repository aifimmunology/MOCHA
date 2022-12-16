test_that("getSampleTileMatrices works on a 1 sample test dataset", {
  
  cellPopulations = c("C2", "C5")
  
  capture.output(
    tilemat <- MOCHA::getSampleTileMatrix(
      MOCHA:::testTileResults,
      cellPopulations = cellPopulations,
      threshold = 0
    ),
    type = "message"
  )
  
  expect_snapshot_value(
     dim(tilemat),
     style = "json2"
  )
  
  expect_equal(
    names(SummarizedExperiment::assays(tilemat)),
    c("C2", "C5")
  )
 
  # expect_snapshot_value(
  #   tilemat,
  #   style = "json2"
  # )
    
})
