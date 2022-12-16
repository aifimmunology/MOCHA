# Code to prepare test data

library(ArchR)
library(MOCHA)
library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(org.Hs.eg.db)

ArchR::getTestProject()
testProj <- ArchR::loadArchRProject("PBMCSmall")

TxDb <- TxDb.Hsapiens.UCSC.hg38.refGene 
Org <- org.Hs.eg.db

testTileResults <- MOCHA::callOpenTiles(
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

usethis::use_data(
  exampleFragments,
  exampleCellColData,
  exampleBlackList,
  internal = FALSE,
  overwrite = TRUE,
  compress = "xz"
)

usethis::use_data(
  testTileResults,
  internal = TRUE,
  overwrite = TRUE,
  compress = "xz"
)
