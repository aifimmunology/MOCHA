test_that("getSampleTileMatrices works on a 1 sample test dataset", {
  
  cellPopulations = c("C1", "C2")
  
  capture.output(
    tilemat <- MOCHA::getSampleTileMatrix(
      MOCHA:::tileResults,
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
    names(assays(tilemat)),
    c("C1","C2")
  )
 
  # expect_snapshot_value(
  #   tilemat,
  #   style = "json2"
  # )
    
})
