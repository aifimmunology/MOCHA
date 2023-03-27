if (dir.exists("../../../HemeTutorial/MOCHA")) {
  test_that("extractRegion works on a 3 sample test dataset", {
    capture.output(
      SampleTileMatrix <- MOCHA::getSampleTileMatrix(
        MOCHA:::testTileResultsMultisample,
        cellPopulations = "all",
        threshold = 0
      )
    )

    SampleTileMatrix@metadata$Directory <- "../../../HemeTutorial/MOCHA"
    cellPopulation <- "C3"
    capture.output(
      countSE <- MOCHA::extractRegion(
        SampleTileObj = SampleTileMatrix,
        cellPopulations = "ALL",
        region = "chr1:18137866-38139912",
        numCores = 30,
        sampleSpecific = FALSE
      )
    )
    expect_snapshot_output(
      countSE
    )
    expect_snapshot_output(
      length(differentials),
      variant = "length"
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

  SampleTileMatrix@metadata$Directory <- "../../../HemeTutorial/MOCHA"
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


