## code to prepare `DATASET` dataset goes here

library(ArchR)
library(MOCHA)
library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(org.Hs.eg.db)

testProj <- ArchR::loadArchRProject("PBMCSmall")
TxDb <- TxDb.Hsapiens.UCSC.hg38.refGene 
Org <- org.Hs.eg.db

tiles <- MOCHA::callOpenTiles(
  ATACFragments = testProj,
  TxDb = TxDb,
  Org = Org,
  cellPopLabel = "Clusters",
  cellPopulations = c("C1", "C2"),
  numCores = 1
)

# Add dataset objects to new environment
ATACFragments <- MOCHA::getPopFrags(
  ArchRProj = testProj,
  metaColumn = "Clusters",
  cellSubsets = c("C1", "C2"),
  region = NULL,
  numCores = 10,
  sampleSpecific = TRUE,
  NormMethod = "nfrags",
  blackList = NULL,
  overlapList = 50
)
cellColData <- ArchR::getCellColData(testProj)
blackList <- ArchR::getBlacklist(testProj)
genome <- ArchR::validBSgenome(ArchR::getGenome(testProj))

usethis::use_data(
  ATACFragments,
  cellColData,
  blackList,
  genome,
  internal = TRUE,
  overwrite = TRUE
)
