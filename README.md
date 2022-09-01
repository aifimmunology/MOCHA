README for scMACS
===============

* * *

Table of Contents
-----------------

* [Introduction](#introduction)
* [Installation](#library)
* [Package Vignette on COVID PASC dataset](#vignette)

* [Contact](#contact)
* [License](#license)

* * *

# <a name="introduction"></a> Introduction
scMACS is an R package containing a novel single-cell peak-calling algorithm that leverages single-cell information to determine whether a particular genomic region is open by calculating two measures of intensities, and using these to call peaks via a hierarchical model. 

# <a name="library"></a> Install package and load library

To install in a HISE IDE, run the following lines of code to install directly from GitHub replacing 'your_token' with your Github Personal Access Token. 
    
    Sys.setenv(GITHUB_PAT='your_token') 
    devtools::install_github("aifimmunology/scMACS")

# <a name="vignette"></a> Usage
Example usage can be found in this [vignette](vignettes/COVID_example.R) copied below:
```R
# Load libraries
library(ArchR)
library(devtools)
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
#    MultiAssayExperiment object and are the primary 
#    input to downstream analyses.
####################################################

SampleTileMatrices <- scMACS::getSampleTileMatrix( 
    tileResults,
    cellPopulations = cellPopulations,
    groupColumn = groupColumn,
    threshold = threshold,
    join = join,
    NAtoZero = FALSE,
    log2Intensity = TRUE
)

#Annotate your open chromatin regions as either a promoter, exonic, intronic, or distal region. Gene names are given for all but distal. 
SampleTileMatrices <- annotateTiles(SampleTilesMatrices)

##########################################

## 4. Plot region 

##########################################

countdf <- extractRegionDF(SampleTileObj, cellTypes = 'CD16 Mono',
                           region = 'chr3:38137866-38139912', 
                           groupColumn = 'InfectionStages',
                           numCores = 30,
                           sampleSpecific = FALSE)
pdf('ExamplePlot.pdf')
plotRegion(countSE = countdf, whichGene = 'MYD88')
dev.off()

```

The output of each function is a [MultiAssayExperiment](https://www.bioconductor.org/packages/devel/bioc/vignettes/MultiAssayExperiment/inst/doc/MultiAssayExperiment.html#overview-of-the-multiassayexperiment-class) organizing your results (either open tiles or sample-tile matrices) by cell population.

Sample-level metadata can be accessed using `colData(yourResults)`. 

Individual results for each cell population can be accessed as follows:
```R
# Get the sample-tile matrix for DC cells
DCcellsMatrix <- getCellTypeMatrix(SampleTileMatrices, 'DC')
# Get the sample-tile matrix for MAIT cells
MAITcellsMatrix <- getCellTypeMatrix(SampleTileMatrices, 'MAIT')

```

The results of a specific cell population from `scMACS::callOpenTiles` (as above) are stored as a [RaggedExperiment](https://bioconductor.org/packages/release/bioc/vignettes/RaggedExperiment/inst/doc/RaggedExperiment.html).

This object stores ranged data alongside the sample-level metadata.

```R
# Get open tiles for DC cells
DCcellRaggedExp <- tileResults[["DC"]]

# DCcellTileResults is a RaggedExperiment
# You can get the sample-level metadata (same as above)
sampleData <- colData(DCcellRaggedExp)

# You can access the underlying GRanges object containing ranges for all samples
tilesGRanges <- rowRanges(DCcellRaggedExp)

```


# <a name="contact"></a> Contact

To contact the developers on issues and feature requests, please contact us via the discussions tab for feature requests, or open issues for any bugs. 
    
# <a name="license"></a> License

scMACS follows the Allen Institute Software License - full information about the license can be found on the LICENSE file. 
