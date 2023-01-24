# This test should only be used for local testing
# with TxDb.Hsapiens.UCSC.hg38.refGene and org.Hs.eg.db installed/
if (
  require("TxDb.Hsapiens.UCSC.hg38.refGene", quietly = TRUE) &&
  require("org.Hs.eg.db", quietly = TRUE) &&
  require("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
) {
  if (require("ArchR", quietly = TRUE) & dir.exists("PBMCSmall")) {
    test_that("We can call peaks by sample from an ArchR project", {
      capture.output(
        testProj <- ArchR::loadArchRProject("PBMCSmall"),
        type = "message"
      )

      TxDb <- TxDb.Hsapiens.UCSC.hg38.refGene
      Org <- org.Hs.eg.db
      capture.output(
        tiles <- MOCHA::callOpenTiles(
          ATACFragments = testProj,
          TxDb = TxDb,
          Org = Org,
          cellPopLabel = "Clusters",
          cellPopulations = c("C2", "C5"),
          numCores = 1,
          outDir = tempdir()
        ),
        type = "message"
      )

      expect_snapshot(
        tiles,
        variant = "ArchR"
      )
    })
  }

  test_that("We can call peaks independent of ArchR", {

    capture.output(
      tiles <- MOCHA::callOpenTiles(
        ATACFragments = MOCHA::exampleFragments,
        cellColData = MOCHA::exampleCellColData,
        blackList = MOCHA::exampleBlackList,
        genome = genome,
        TxDb = TxDb,
        Org = Org,
        outDir = tempdir(),
        cellPopLabel = "Clusters",
        cellPopulations = c("C2", "C5"),
        numCores = 1
      ),
      type = "message"
    )

    expect_snapshot(
      tiles,
      variant = "list"
    )
  })
  
  test_that("We throw a warning when a sample has less than 5 cells", {
    
      sample1frags <- GenomicRanges::GRanges(
        seqnames = Rle(c("chr1"), c(1)),
        ranges = IRanges(c(760101:760110), end = c(760111:760120), names = head(letters, 10)),
        strand = "*",
        RG = c("c1","c2","c2","c3","c3","c4","c4","c4","c5","c5")
      )
      sample2frags <- GenomicRanges::GRanges(
        seqnames = Rle(c("chr1"), c(1)),
        ranges = IRanges(c(760101:760110), end = c(760111:760120), names = head(letters, 10)),
        strand = "*",
        RG = c("c6","c7","c7","c8","c8","c8","c9","c9","c9","c9")
      )
      tiny_fragments <- GenomicRanges::GRangesList(sample1frags, sample2frags)
      names(tiny_fragments) <- c("t_cd8_temra#sample1", "t_cd8_temra#sample2")
      
      tiny_cellColData <- data.frame(
        Sample = c(rep("sample1", 5), rep("sample2", 4)),
        cellPop = rep("t_cd8_temra", 9)
      )
      rownames(tiny_cellColData) <- c("c1","c2","c3","c4","c5","c6","c7","c8","c9")
    
      expect_warning(tiles <- MOCHA::callOpenTiles(
          ATACFragments = tiny_fragments,
          cellColData = tiny_cellColData,
          blackList = MOCHA::exampleBlackList,
          genome = genome,
          TxDb = TxDb,
          Org = Org,
          outDir = tempdir(),
          cellPopLabel = "cellPop",
          cellPopulations = c("t_cd8_temra"),
          studySignal = 10, # manually provide or have nFrags col in cellColData
          numCores = 1, verbose=TRUE
      ))

    }
  )
}
