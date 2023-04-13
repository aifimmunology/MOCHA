test_that("bulkLSI works on a 3 sample test dataset", {
  cellPopulations <- c("C2", "C3")
  capture.output(
    SampleTileMatrix <- MOCHA::getSampleTileMatrix(
      MOCHA:::testTileResultsMultisample,
      cellPopulations = cellPopulations,
      threshold = 0
    ),
    type = "message"
  )
  
  # Spoof data to remove NAs
  c2 <- SummarizedExperiment::assays(SampleTileMatrix)[[1]] 
  set.seed(1)
  c2[,2] <- runif(dim(c2)[1], min=1100, max=8700)
  set.seed(2)
  c2[,3] <- runif(dim(c2)[1], min=1100, max=8700)
  SummarizedExperiment::assays(SampleTileMatrix)[[1]] <- c2
  
  # Spoof data to remove NAs
  c3 <- SummarizedExperiment::assays(SampleTileMatrix)[[2]] 
  set.seed(3)
  c3[,2] <- runif(dim(c3)[1], min=1100, max=8700)
  set.seed(4)
  c3[,3] <- runif(dim(c3)[1], min=1100, max=8700)
  SummarizedExperiment::assays(SampleTileMatrix)[[2]] <- c3

  capture.output(
    LSIObj <- MOCHA::bulkLSI(SampleTileMatrix, cellType = 'all', componentNumber=2)
  )
  
  capture.output(
    UMAPvalues <- MOCHA::bulkUMAP(LSIObj, components = c(1:2), n_neighbors=4)
  )
  
  
  expect_snapshot(
    LSIObj
  )
  
  expect_snapshot(
    UMAPvalues
  )
  
})

test_that("bulkLSI errors NA columns", {
  cellPopulations <- c("C2", "C3")
  capture.output(
    SampleTileMatrix <- MOCHA::getSampleTileMatrix(
      MOCHA:::testTileResultsMultisample,
      cellPopulations = cellPopulations,
      threshold = 0
    ),
    type = "message"
  )

  expect_error(
    LSIse <- MOCHA::bulkLSI(SampleTileMatrix, cellType = 'all', componentNumber=2),
    "starting vector near the null space"
  )
})