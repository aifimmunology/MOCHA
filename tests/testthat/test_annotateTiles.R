test_that("annotateTiles works on a 1 sample test dataset", {
  cellPopulations <- c("C2", "C5")
  capture.output(
    SampleTileMatrix <- MOCHA::getSampleTileMatrix(
      MOCHA:::testTileResults,
      cellPopulations = cellPopulations,
      threshold = 0
    ),
    type = "message"
  )

  capture.output(
    STM <- MOCHA::annotateTiles(SampleTileMatrix,
      promoterRegion = c(2000, 100)
    ),
    type = "message"
  )

  expect_snapshot(
    rowRanges(STM),
    variant = "1sample"
  )
})

test_that("annotateTiles works on a 3 sample test dataset", {
  cellPopulations <- c("C2", "C3")
  capture.output(
    SampleTileMatrix <- MOCHA::getSampleTileMatrix(
      MOCHA:::testTileResultsMultisample,
      cellPopulations = cellPopulations,
      threshold = 0
    ),
    type = "message"
  )

  capture.output(
    STM <- MOCHA::annotateTiles(SampleTileMatrix,
      promoterRegion = c(2000, 100)
    )
  )

  expect_snapshot(
    rowRanges(STM),
    variant = "3sample"
  )
})

if (requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
  requireNamespace("TxDb.Hsapiens.UCSC.hg38.refGene", quietly = TRUE)) {
  test_that("annotateTiles works on a GRanges", {
    capture.output(
      ranges <- MOCHA::annotateTiles(rowRanges(SampleTileMatrix),
        TxDb = TxDb.Hsapiens.UCSC.hg38.refGene,
        Org = org.Hs.eg.db,
        promoterRegion = c(2000, 100)
      ),
      type = "message"
    )

    expect_snapshot(
      ranges,
      variant = "3sample"
    )
  })
}