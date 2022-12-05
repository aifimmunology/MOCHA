# Code to prepare test data

library(ArchR)
library(MOCHA)
library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(org.Hs.eg.db)

ArchR::getTestProject()
testProj <- ArchR::loadArchRProject("PBMCSmall")

TxDb <- TxDb.Hsapiens.UCSC.hg38.refGene 
Org <- org.Hs.eg.db

exampleTileResults <- MOCHA::callOpenTiles(
  ATACFragments = testProj,
  TxDb = TxDb,
  Org = Org,
  cellPopLabel = "Clusters",
  cellPopulations = c("C2", "C5"),
  numCores = 1
)

exampleFragments <- MOCHA::getPopFrags(
  ArchRProj = testProj,
  metaColumn = "Clusters",
  cellSubsets = c("C2", "C5"),
  region = NULL,
  numCores = 10,
  sampleSpecific = TRUE,
  NormMethod = "nfrags",
  blackList = NULL,
  overlapList = 50
)

exampleCellColData <- ArchR::getCellColData(testProj)
exampleCellColData <- cellColData[c("Sample", "nFrags", "Clusters")]

exampleBlackList <- ArchR::getBlacklist(testProj)
exampleGenome <- ArchR::validBSgenome(ArchR::getGenome(testProj))

usethis::use_data(
  exampleTileResults,
  exampleFragments,
  exampleCellColData,
  exampleBlackList,
  exampleGenome,
  internal = FALSE,
  overwrite = TRUE
)
