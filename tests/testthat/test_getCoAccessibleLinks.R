test_that("FindCoAccessibleLinks works on a 1 sample test dataset", {

  capture.output(
    SampleTileMatrix <- MOCHA::getSampleTileMatrix(
      MOCHA:::tileResults,
      cellPopulations = cellPopulations,
      threshold = 0
    )
  )

  cellPopulation <- "C1"
  regions <- MOCHA::StringsToGRanges(c(
    "chr1:102368000-102368499",
    "chr1:101873000-101873499"
  ))
  links <- MOCHA::getCoAccessibleLinks(SampleTileMatrix,
    cellPopulation,
    regions,
    verbose = FALSE
  )

  expect_false(any((links$Peak1 == links$Peak2)))
  expect_snapshot_output(
    links,
    variant = "tiles"
  ) # change variant to "basic" to compare with previous co_accessibility.R
  expect_snapshot_output(
    nrow(links),
    variant = "nrows"
  )
})
