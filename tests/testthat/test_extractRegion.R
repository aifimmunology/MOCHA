if (dir.exists("../../../HemeTutorial/MOCHA")) {
  test_that("extractRegion works on a 3 sample test dataset", {
    capture.output(
      SampleTileMatrix <- MOCHA::getSampleTileMatrix(
        MOCHA:::testTileResultsMultisample,
        cellPopulations = "all",
        threshold = 0
      )
    )
    
    # Since there are only cells for two samples in population "C3", 
    # remove the empty sample.
    # 
    # $CellCounts
    # cellTypeLabelList
    #                     C2  C3
    # scATAC_BMMC_R1      379  22
    # scATAC_CD34_BMMC_R1   0 302
    # scATAC_PBMC_R1        0   0
    # 
    SampleTileMatrix <- SampleTileMatrix[
      , SampleTileMatrix$Sample %in% c("scATAC_BMMC_R1", "scATAC_CD34_BMMC_R1")
    ]
    
    SampleTileMatrix@metadata$Directory <- "../../../HemeTutorial/MOCHA"
    cellPopulation <- "C3"
    capture.output(
      countSE <- MOCHA::extractRegion(
        SampleTileObj = SampleTileMatrix,
        cellPopulations = cellPopulation,
        region = "chr1:18137866-38139912",
        numCores = 30,
        sampleSpecific = FALSE
      )
    )
    expect_snapshot_output(
      countSE
    )
  })
  
  test_that("extractRegion errors when there is no fragment coverage for a cell population", {
    capture.output(
      SampleTileMatrix <- MOCHA::getSampleTileMatrix(
        MOCHA:::testTileResultsMultisample,
        cellPopulations = "all",
        threshold = 0
      )
    )
    
    SampleTileMatrix@metadata$Directory <- "../../../HemeTutorial/MOCHA"
    cellPopulation <- "C3"
    expect_error(
      countSE <- MOCHA::extractRegion(
        SampleTileObj = SampleTileMatrix,
        cellPopulations = cellPopulation,
        region = "chr1:18137866-38139912",
        numCores = 30,
        sampleSpecific = FALSE
      ),
      "There is no fragment coverage for cell population"
    )
  })
}

test_that("extractRegion errors when coverage files aren't saved locally", {
  capture.output(
    SampleTileMatrix <- MOCHA::getSampleTileMatrix(
      MOCHA:::testTileResultsMultisample,
      cellPopulations = "all",
      threshold = 0
    )
  )

  SampleTileMatrix@metadata$Directory <- "../../../HemeTutorial/idontexist"
  cellPopulation <- "C3"
  expect_error(
    MOCHA::extractRegion(
      SampleTileObj = SampleTileMatrix,
      cellPopulations = "ALL",
      region = "chr1:18137866-38139912",
      numCores = 30,
      sampleSpecific = FALSE
    ),
    "does not exist"
  )
})





