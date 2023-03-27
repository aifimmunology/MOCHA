test_that("testCoAccessibleLinks works on a 1 sample test dataset", {
  cellPopulations <- c("C2", "C5")
  capture.output(
    SampleTileMatrix <- MOCHA::getSampleTileMatrix(
      MOCHA:::testTileResults,
      cellPopulations = cellPopulations,
      threshold = 0
    )
  )
  
  cellPopulation <- "C2"
  regions <- MOCHA::StringsToGRanges(c(
    "chr1:101873000-101873499"
  ))
  capture.output(
    links <- MOCHA::getCoAccessibleLinks(SampleTileMatrix,
                                         cellPopulation,
                                         regions,
                                         verbose = FALSE
    )
  )
  
  capture.output(
    results <- MOCHA::testCoAccessibilityRandom(
      SampleTileMatrix,
      tile1 = links$Tile1,
      tile2 = links$Tile2,
      numCores = 1,
      ZI = TRUE,
      backNumber = 1000,
      verbose = FALSE
    )
  )
  
  expect_snapshot(
    results,
    variant = "1sample"
  )
  
})

test_that("testCoAccessibleLinks works on a 3 sample test dataset", {
  cellPopulations <- c("C2", "C3")
  capture.output(
    SampleTileMatrix <- MOCHA::getSampleTileMatrix(
      MOCHA:::testTileResultsMultisample,
      cellPopulations = cellPopulations,
      threshold = 0
    )
  )
  
  cellPopulation <- "C2"
  regions <- MOCHA::StringsToGRanges(c(
    "chr1:101775000-101775499"
  ))
  
  capture.output(
    links <- MOCHA::getCoAccessibleLinks(SampleTileMatrix,
                                         cellPopulation,
                                         regions,
                                         verbose = FALSE
    )
  )
  
  capture.output(
    results <- MOCHA::testCoAccessibilityRandom( 
      SampleTileMatrix,
      tile1 = links$Tile1,
      tile2 = links$Tile2,
      numCores = 1,
      ZI = TRUE,
      backNumber = 1000,
      verbose = FALSE
    )
  )
  
  expect_snapshot(
    results,
    variant = "3sample"
  )

})
