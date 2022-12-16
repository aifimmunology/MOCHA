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

    TxDb <- TxDb.Hsapiens.UCSC.hg38.refGene
    Org <- org.Hs.eg.db
    genome <- BSgenome.Hsapiens.UCSC.hg19
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
}