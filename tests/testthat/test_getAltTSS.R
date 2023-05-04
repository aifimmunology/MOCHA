if (
  require("TxDb.Hsapiens.UCSC.hg38.refGene") &&
  require("org.Hs.eg.db")
) {
  test_that("getAltTSS works with different input formats", {

    testPeakDT <- MOCHA:::testPeaks
    set.seed(2023)
    testPeakDT$Log2FC_C <- sample(seq(-4, 4, .1), nrow(testPeakDT), replace=TRUE)

    testPeakDF <- as.data.frame(testPeakDT)

    testPeakGR <- MOCHA::differentialsToGRanges(
      testPeakDT, tileColumn = "tileID"
    )

    # This will drop these from metadata: tileID seqnames start end strand
    testPeakGR2 <- GenomicRanges::makeGRangesFromDataFrame(
      as.data.frame(testPeakDT),
      keep.extra.columns = TRUE
    )

    expect_true(all(
      nrow(testPeakDF) == nrow(testPeakDT),
      nrow(testPeakDF) == length(testPeakGR)
    ))

    # data.table
    tpeaksDT <- MOCHA::getAltTSS(
      completeDAPs = testPeakDT,
      returnAllTSS = FALSE,
      nuancedTSS = FALSE,
      nuancedTSSGap = 150,
      threshold = 0.2,
      TxDb = TxDb.Hsapiens.UCSC.hg38.refGene,
      Org = org.Hs.eg.db
    )

    # data.frame
    tpeaksDF <- MOCHA::getAltTSS(
      completeDAPs = testPeakDF,
      returnAllTSS = FALSE,
      nuancedTSS = FALSE,
      nuancedTSSGap = 150,
      threshold = 0.2,
      TxDb = TxDb.Hsapiens.UCSC.hg38.refGene,
      Org = org.Hs.eg.db
    )

    # GRanges
    tpeaksGR <- MOCHA::getAltTSS(
      completeDAPs = testPeakGR,
      returnAllTSS = FALSE,
      nuancedTSS = FALSE,
      nuancedTSSGap = 150,
      threshold = 0.2,
      TxDb = TxDb.Hsapiens.UCSC.hg38.refGene,
      Org = org.Hs.eg.db
    )

    expect_true(all(
      length(tpeaksDT) == length(tpeaksDF),
      length(tpeaksDT) == length(tpeaksGR)
    ))

    expect_snapshot_output(
      tpeaksDT
    )

  })

  test_that("getAltTSS works with returnAllTSS = TRUE", {

    testPeakDT <- MOCHA:::testPeaks
    set.seed(2023)
    testPeakDT$Log2FC_C <- sample(seq(-4, 4, .1), nrow(testPeakDT), replace=TRUE)

    tpeaks <- MOCHA::getAltTSS(
      completeDAPs = testPeakDT,
      returnAllTSS = TRUE,
      nuancedTSS = FALSE,
      nuancedTSSGap = 150,
      threshold = 0.2,
      TxDb = TxDb.Hsapiens.UCSC.hg38.refGene,
      Org = org.Hs.eg.db
    )

    expect_snapshot_output(
      tpeaks,
      variant = "returnAllTSS"
    )

  })

  test_that("getAltTSS works with nuancedTSS = TRUE", {

    testPeakDT <- MOCHA:::testPeaks
    set.seed(2023)
    testPeakDT$Log2FC_C <- sample(seq(-4, 4, .1), nrow(testPeakDT), replace=TRUE)

    tpeaksDT <- MOCHA::getAltTSS(
      completeDAPs = testPeakDT,
      returnAllTSS = FALSE,
      nuancedTSS = TRUE,
      nuancedTSSGap = 150,
      threshold = 0.2,
      TxDb = TxDb.Hsapiens.UCSC.hg38.refGene,
      Org = org.Hs.eg.db
    )

    expect_snapshot_output(
      tpeaksDT,
      variant = "nuancedTSS"
    )

  })
}

