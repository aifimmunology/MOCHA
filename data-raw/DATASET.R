# Code to prepare test data

library(ArchR)
library(MOCHA)
library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(org.Hs.eg.db)

ArchR::getTestProject()
testProj <- ArchR::loadArchRProject("PBMCSmall")

TxDb <- "TxDb.Hsapiens.UCSC.hg38.refGene"
Org <- "org.Hs.eg.db"

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

# Uncomment if HemeFragments is already downloaded
# inputFiles <- c(scATAC_BMMC_R1="HemeFragments/scATAC_BMMC_R1.fragments.tsv.gz", 
#                 scATAC_PBMC_R1="HemeFragments/scATAC_PBMC_R1.fragments.tsv.gz", 
#                 scATAC_CD34_BMMC_R1="HemeFragments/scATAC_CD34_BMMC_R1.fragments.tsv.gz")

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
  copyArrows = TRUE # This is recommended so that you maintain an unaltered copy for later usage.
)
proj <- ArchR::filterDoublets(ArchRProj = proj)

proj <- ArchR::addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")

proj <- ArchR::addClusters(input = proj, reducedDims = "IterativeLSI")

proj <- ArchR::addImputeWeights(proj)

proj <- ArchR::saveArchRProject(ArchRProj = proj)

# Uncomment for reproducing if HemeTutorial already exists
# proj <- ArchR::loadArchRProject("HemeTutorial")

testTileResultsMultisample <- MOCHA::callOpenTiles(
  ATACFragments = proj,
  outDir = "../HemeTutorial/MOCHA/",
  TxDb = TxDb,
  Org = Org,
  cellPopLabel = "Clusters",
  cellPopulations = c("C2", "C3"),
  numCores = 4
)
experiments <- MultiAssayExperiment::experiments(testTileResultsMultisample)

for (x in 1:length(experiments)) {
  print(x)
  exp <- testTileResultsMultisample[[x]]
  randIdx <- sort(sample(1:length(exp), 50000, replace = FALSE))

  testTileResultsMultisample[[x]] <- exp[randIdx]
}

# Motif list for MotifEnrichment
requireNamespace("chromVAR", quietly = TRUE)
requireNamespace("chromVARmotifs", quietly = TRUE)
requireNamespace("motifmatchr", quietly = TRUE)
requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
library(BSgenome.Hsapiens.UCSC.hg19)

testSampleTileMatrix <- MOCHA::getSampleTileMatrix(
  testTileResultsMultisample,
  cellPopulations = c("C2", "C3"),
  threshold = 0
)

STM <- MOCHA::addMotifSet(
  testSampleTileMatrix,
  chromVARmotifs::human_pwms_v2,
  genome = BSgenome.Hsapiens.UCSC.hg19,
  motifSetName = 'CISBP'
)
posList <- STM@metadata$CISBP
smallPosList <- posList[1:10]

# Create a small peakset based on plasmablast peaks from Samir
samplePeaks <- readRDS(
  "/Users/imran.mcgrath/Downloads/1137_cells_sample_1_peaks.rds"
)[[1]]
samplePeaks <- samplePeaks[samplePeaks$seqnames %in% c("chr1", "chr2")]
samplePeaks <- samplePeaks[samplePeaks$peak == TRUE]
miniSamplePeaks <- samplePeaks[
  sample(nrow(samplePeaks), size=floor(nrow(samplePeaks)/3))
]
# Spoof FDR column, making 10% of FDR <0.1
miniSamplePeaks$FDR <- round(runif(nrow(miniSamplePeaks)),2)


testPeaks <- miniSamplePeaks
testPosList <- smallPosList

# Uncomment to preserve current versions of test data
# while adding new metadata
newmeta1 <- testTileResults@metadata
testTileResults <- MOCHA:::testTileResults
metadata(testTileResults) <- newmeta1

newmeta2 <- testTileResultsMultisample@metadata
testTileResultsMultisample <- MOCHA:::testTileResultsMultisample
metadata(testTileResultsMultisample) <- newmeta2

testPeaks <- MOCHA:::testPeaks
testPosList <- MOCHA:::testPosList

# Save internal data
usethis::use_data(
  testTileResults,
  testTileResultsMultisample,
  testPeaks,
  testPosList,
  internal = TRUE,
  overwrite = TRUE,
  compress = "xz"
)
