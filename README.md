README for scMACS
===============

* * *

Table of Contents
-----------------

* [Introduction](#introduction)
* [Install package and load library](#library)
* [Tutorials](#example-main)
    * [Tutorial-1: Peak-calling on Entire Cell Populations](#example1)
    * [Tutorial-2: Sample-specific Peak-calling](#example2)
    * [Tutorial-3: Differential Accessibility](#example3)
* [Contact](#contact)
* [License](#license)

* * *

# <a name="introduction"></a> Introduction
scMACS is an R package for single-cell peak-calling algorithm that quantifies two measures of reads intensities in scATAC data to call peaks. It leverages single-cell information to determine whether a particular genomic region is open by calculating two measures of intensities, and using these to call peaks. 

# <a name="library"></a> Install package and load library

To install library, simply run the following lines of code to install either from source
   
   
    install.packages("scMACS_0.1.0.tar.gz", repos = NULL, type ="source")
    library("scMACS")

or directly from GitHub

    devtools::install_github("aifimmunology/scMACS")
    library("scMACS")
    
# <a name="example-main"></a> Tutorials
## <a name="example1"></a> Tutorial-1: Peak-calling on Entire Cell Populations

This tutorial demonstrates how to call peaks for a given cell population. 

### Load Library
   
    #Load scMACS library
    library(scMACS)
    library(data.table)
    library(ArchR)
    library(GenomicRanges, lib.loc=libLoc)
    library(plyranges, lib.loc=libLoc)
    
### Load data and assign parameters
    
    ### parameters
    cellSubsets='ALL'
    cellCol_label_name="predictedGroup_Col2.5"
    returnAllPeaks=FALSE
    numCores=10


    ### load ArchR Project
    start <- Sys.time()
    ArchRProj <- ArchR::loadArchRProject('/path/to/ArchRProject')

    peak_list <- scMACS::callPeaks(ArchRProj=ArchRProj,
        cellSubsets=cellSubsets,
        cellCol_label_name=cellCol_label_name,
        returnAllPeaks=returnAllPeaks,
        numCores=numCores
    )