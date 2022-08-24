####################################################
# 0. Load libraries, ArchR project, and annotation
#    databases. Optionally filter the ArchR project
#    to a subset of samples.
####################################################

library(ArchR)
library(scMACS)

# You should substitute this with your own ArchR project.
# You must have completed cell labeling with your ArchR project.
myArchRProj <- ArchR::loadArchRProject("/home/jupyter/FullCovid")

# Define your annotation package for TxDb object(s)
# and genome-wide annotation
# Here our samples are human using hg38 as a reference.
# For more info: https://bioconductor.org/packages/3.15/data/annotation/
library(TxDb.Hsapiens.UCSC.hg38.refGene)
library(org.Hs.eg.db)
TxDb <- TxDb.Hsapiens.UCSC.hg38.refGene
Org <- org.Hs.eg.db

# Optional: Filter your ArchR project by sample.
# For our example we filter ArchR Project to three samples from 
# each Covid Status (3 Positive, 3 Negative).
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
#    These should be set according to your ArchR 
#    project and investgative question.
#
#    For more details on each of these parameters, 
#    view the help pages for each function using 
#    ?callOpenTiles and ?getSampleTileMatrix
####################################################

# Parameters for calling open tiles.
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
    TxDb = TxDb,
    Org = Org,
    numCores = numCores
)

####################################################
# 3. Get consensus sample-tile matrices
#    for all cell populations.
#    These matrices are organized by cell population
#    RangedSummarizedExperiment object and are the 
#    primary input to downstream analyses.
####################################################

SampleTileMatrices <- scMACS::getSampleTileMatrix( 
    tileResults,
    cellPopulations = cellPopulations,
    groupColumn = groupColumn,
    threshold = threshold,
    join = join,
    NAtoZero = TRUE,
    log2Intensity = TRUE
)

####################################################
# 4. Add gene annotations to our SampleTileMatrices.
#    This info will aid further downstream analyses
#    but is not required for differential 
#    accessibility.
#    This function can also take any GRanges object
#    and add annotations to its metadata.
####################################################

SampleTileMatricesAnnotated <- scMACS::annotateTiles( 
  SampleTileMatrices
)

####################################################
# 5. Get differential accessibility for specific 
#    cell populations. Here we are comparing MAIT  
#    cells between samples where our groupColumn 
#    "COVID_status" is Positive (our foreground) 
#    to Negative samples (our background).
####################################################

differentials <- scMACS::getDifferentialAccessibleTiles(
    SampleTileObj = SampleTileMatrices,
    cellPopulation = "MAIT",
    groupColumn = groupColumn,
    foreground = "Positive",
    background = "Negative",
    fdrToDisplay = 0.4,
    numCores = numCores
)
