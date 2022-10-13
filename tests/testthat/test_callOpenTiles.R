# This test should only be used for local testing
if (requireNamespace("ArchR", quietly = TRUE) & dir.exists("PBMCSmall")) {
  test_that("We can call peaks by sample from an ArchR project", {
    capture.output(testProj <- ArchR::loadArchRProject("PBMCSmall"), type = "message")
  
    TxDb <- TxDb.Hsapiens.UCSC.hg38.refGene 
    Org <- org.Hs.eg.db
    capture.output(
      tiles <- MOCHA::callOpenTiles(
        ATACFragments = testProj,
        TxDb = TxDb,
        Org = Org,
        cellPopLabel = "Clusters",
        cellPopulations = c("C1", "C2"),
        numCores = 1
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
  capture.output(
    tiles <- MOCHA::callOpenTiles(
      MOCHA:::ATACFragments,
      MOCHA:::cellColData,
      MOCHA:::blackList,
      MOCHA:::genome,
      TxDb = TxDb,
      Org = Org,
      outDir = "./test_out",
      cellPopLabel = "Clusters",
      cellPopulations = c("C1", "C2"),
      numCores = 1
    ),
    type = "message"
  )
  unlink("./test_out", recursive = TRUE)
  
  expect_snapshot(
    tiles,
    variant = "list"
  )
})

