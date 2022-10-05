if (requireNamespace("ArchR", quietly = TRUE)) {
  if(dir.exists("PBMCSmall")){
    
    test_that("We can call peaks by sample from an ArchR project", {
      capture.output(testProj <- ArchR::loadArchRProject("PBMCSmall"), type = "message")
    
      TxDb <- TxDb.Hsapiens.UCSC.hg38.refGene 
      Org <- org.Hs.eg.db
      capture.output(
        tiles <- MOCHA::callOpenTiles(
          ArchRProj = testProj,
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
  
} else {
  
  test_that("We can call peaks independent of ArchR", {
    
    TxDb <- TxDb.Hsapiens.UCSC.hg38.refGene 
    Org <- org.Hs.eg.db
    capture.output(
      tiles <- MOCHA::callOpenTiles(
        ArchRProj = testProj,
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

