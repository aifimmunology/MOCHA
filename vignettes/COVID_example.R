# Load libraries
library(ArchR)
library(devtools)
install("/home/jupyter/scMACS")
library(scMACS)

####################################################
# 1. Setting Parameters
####################################################

# Load the ArchR Project 
covidArchR <- ArchR::loadArchRProject("/home/jupyter/FullCovid")

# For example: filter ArchR Project to three samples from each Covid Status
samplesToKeep <- c(
  "B011-AP0C1W3", 
  "B011-AP0C1W8", 
  "B011-AP0C2W1",
  "B025_FSQAAZ0BZZS-01",
  "B025_FSQAAZ0C0YJ-01",
  "B025_FSQAAZ0C00P-07"
)
idxSample <- BiocGenerics::which(covidArchR$Sample %in% samplesToKeep)
cellsSample <- covidArchR$cellNames[idxSample]
covidArchR <- covidArchR[cellsSample, ]

# Parameters for calling open tiles
cellPopLabel <- "CellSubsets"
cellPopulations <- c("CD16 Mono", "DC", "MAIT")
numCores <- 20

# Parameters for downstream analysis
threshold <- 0.2
groupColumn <- "COVID_status"
join <- "union"

####################################################
# 2. Call open tiles (main peak calling step)
#    and get sample-tile matrices
#    for all specified cell populations
####################################################

tileResults <- scMACS::callOpenTiles( 
    covidArchR,
    cellPopLabel = cellPopLabel,
    cellPopulations = cellPopulations,
    numCores = numCores
)


####################################################
# 3. Get consensus sample-tile matrices
#    for all cell populations.
#    These matrices are organized by cell population
#    MultiAssayExperiment object and are the primary 
#    input to downstream analyses.
####################################################

SampleTileMatrices <- scMACS::getSampleTileMatrix( 
    tileResults,
    cellPopulations = 'ALL',
    groupColumn = groupColumn,
    threshold = threshold,
    join = join,
    NAtoZero = FALSE,
    log2Intensity = TRUE
)