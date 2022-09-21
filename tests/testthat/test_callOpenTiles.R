if(dir.exists("PBMCSmall")){
  
  test_that("We can call peaks by sample", {
    capture.output(testProj <- loadArchRProject("PBMCSmall"), type = "message")
    outdir <- dirname(getOutputDirectory(testProj))
  
    TxDb <- TxDb.Hsapiens.UCSC.hg38.refGene 
    Org <- org.Hs.eg.db
    capture.output(
      tiles <- MOCHA::callOpenTiles(
        testProj,
        TxDb = TxDb,
        Org = Org,
        cellPopLabel = "Clusters",
        cellPopulations = c("C1", "C2"),
        numCores = 1
      ),
      type = "message"
    )
    
    # expect_snapshot_value(
    #   tiles,
    #   style = "json2"
    # )
  })
}
