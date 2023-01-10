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
  outDir = NULL,
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


########## Internal data ##########
set.seed(1)
addArchRThreads(threads = 2)
inputFiles <- ArchR::getTutorialData("Hematopoiesis")
addArchRGenome("hg19")

ArrowFiles <- ArchR::createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 4, # Dont set this too high because you can always increase later
  filterFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

doubScores <- ArchR::addDoubletScores(
  input = ArrowFiles,
  k = 10, # Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", # Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)

proj <- ArchR::ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "HemeTutorial",
  copyArrows = TRUE # This is recommened so that you maintain an unaltered copy for later usage.
)
proj <- ArchR::filterDoublets(ArchRProj = proj)

proj <- ArchR::addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")

proj <- ArchR::addClusters(input = proj, reducedDims = "IterativeLSI")

proj <- ArchR::addImputeWeights(proj)

proj <- ArchR::saveArchRProject(ArchRProj = proj)

testTileResultsMultisample <- MOCHA::callOpenTiles(
  ATACFragments = proj,
  outDir = NULL,
  TxDb = TxDb,
  Org = Org,
  cellPopLabel = "Clusters",
  cellPopulations = c("C2", "C3"),
  numCores = 4
)
experiments <- MultiAssayExperiment::experiments(testTileResultsMultisample)

testSampleTileMatrix <- MOCHA::getSampleTileMatrix(
  testTileResultsMultisample,
  cellPopulations = cellPopulations,
  threshold = 0
)

for (x in 1:length(experiments)) {
  print(x)
  exp <- testTileResultsMultisample[[x]]
  randIdx <- sort(sample(1:length(exp), 50000, replace = FALSE))

  testTileResultsMultisample[[x]] <- exp[randIdx]
}

usethis::use_data(
  testTileResults,
  testTileResultsMultisample,
  internal = TRUE,
  overwrite = TRUE,
  compress = "xz"
)
