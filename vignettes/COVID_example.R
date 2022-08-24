# Load libraries
library(ArchR)
library(devtools)
setwd("scMACS")
install()

library(scMACS)

# Load the ArchR Project
# You should substitute this with your own ArchR project.
# You must have completed cell labeling with your ArchR project.
myArchRProj <- ArchR::loadArchRProject("/home/jupyter/FullCovid")

# For our example: filter ArchR Project to three samples from each Covid Status
# This is not strictly necessary for your ArchR Project 
samplesToKeep <- c(
  "B011-AP0C1W3", 
  "B011-AP0C1W8", 
  "B011-AP0C2W1",
  "B025_FSQAAZ0BZZS-01",
  "B025_FSQAAZ0C0YJ-01",
  "B025_FSQAAZ0C00P-07"
)
idxSample <- BiocGenerics::which(myArchRProj$Sample %in% samplesToKeep)
cellsSample <- myArchRProj$cellNames[idxSample]
myArchRProj <- myArchRProj[cellsSample, ]

####################################################
# 1. Setting Parameters
#    These should be set according to your ArchR project
#    and investgative question.
#
#    For more details on each of these parameters, 
#    view the help pages for each function using 
#    ?callOpenTiles and ?getSampleTileMatrix
####################################################

# Parameters for calling open tiles
cellPopLabel <- "CellSubsets" 
cellPopulations <- c("MAIT", "CD16 Mono", "DC")
numCores <- 20

# Parameters for generating the sample-tile matrices
threshold <- 0.2
groupColumn <- "COVID_status"
join <- "union"

####################################################
# 2. Call open tiles (main peak calling step)
#    and get sample-tile matrices
#    for all specified cell populations
####################################################

tileResults <- scMACS::callOpenTiles( 
    myArchRProj,
    cellPopLabel = cellPopLabel,
    cellPopulations = cellPopulations,
    numCores = numCores
)

####################################################
# 3. Get consensus sample-tile matrices
#    for all cell populations.
#    These matrices are organized by cell population
#    RangedSummarizedExperiment object and are the 
#    primary input to downstream analyses.
####################################################
install()

library(scMACS)
SampleTileMatrices <- scMACS::getSampleTileMatrix( 
    tileResults,
    cellPopulations = cellPopulations,
    groupColumn = groupColumn,
    threshold = threshold,
    join = join,
    NAtoZero = TRUE,
    log2Intensity = TRUE
)
traceback()
####################################################
# 3. Get DAPs for specific cell populations.
#    
####################################################

differentials <- getDifferentialAccessibleTiles(
    SampleTileObj = SampleTileMatrices,
    cellPopulation = "MAIT",
    groupColumn = "COVID_status",
    foreground = "Positive",
    background = "Negative",
    fdrToDisplay = 0.2,
    numCores = 20
)
traceback()
print("Done!")