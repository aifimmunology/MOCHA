# Load libraries
library(devtools)
install("/home/jupyter/scMACS")
library(scMACS)
library(ArchR)

####################################################
# 1. Setting Parameters
####################################################

# Load the ArchR Project 
covidArchR <- ArchR::loadArchRProject("/home/jupyter/FullCovid")

# For example: filter ArchR Project to three samples
samplesToKeep <- c("B011-AP0C1W3", "B011-AP0C1W8", "B011-AP0C2W1")
idxSample <- BiocGenerics::which(covidArchR$Sample %in% samplesToKeep)
cellsSample <- covidArchR$cellNames[idxSample]
covidArchR <- covidArchR[cellsSample, ]

# Parameters for calling open tiles
sampleLabel <- "Sample"
cellPopulations <- c("DC", "MAIT")
cellPopLabel <- "CellSubsets"
numCores <- 20

# Parameters for downstream analysis
cellPopulation <- "DC"
reproducibility <- 0

####################################################
# 2. Call open tiles (main peak calling step)
#    Done once for all specified cell populations
####################################################

tileResults <- scMACS::callOpenTiles(
    covidArchR,
    cellPopLabel = cellPopLabel,
    cellPopulations = cellPopulations,
    sampleLabel = sampleLabel,
    numCores = numCores
)

####################################################
# 3. Get reproducible sample-peak matrix
#    Done for each cell population individually
####################################################

reproducibleTiles <- scMACS::getReproducibleTiles(
    tileResults,
    cellPopulation = cellPopulation,
    reproducibility = reproducibility
)

sampleTileMatrix <- scMACS::getSampleTileMatrix(
    tileResults,
    cellPopulation = cellPopulation,
    reproducibleTiles,
    NAtoZero = TRUE
)