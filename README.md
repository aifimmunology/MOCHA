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
scMACS is an R package containing a novel single-cell peak-calling algorithm that leverages single-cell information to determine whether a particular genomic region is open by calculating two measures of intensities, and using these to call peaks via a hierarchical model. 

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
   
    #Load scMACS and accompanying libraries
    library(scMACS)
    library(data.table)
    library(ArchR)
    library(GenomicRanges)
    library(plyranges)
    
### Load data and assign parameters
    
    ################################################################################
    ### To call peaks using scMACS
    ### the user must input 4 different
    ### parameters: 
    
    ###  @param cellSubsets vector of strings. Cell subsets for which to call peaks. 
    ###         Optional, 'ALL' defaults to all cell populations in metadata file. 
    
    ###  @param cellCol_label_name string indicating which column in the metadata 
    ###         file contains the cell population label

    ###  @param returnAllPeaks boolean. Indicates whether scMACS should return object 
    ###         containing all genomic regions or just the positive (+) called peaks. 
    
    ###  @param numCores integer. Number of cores to parallelize peak-calling 
    ###        across multiple cell populations 
    
    ################################################################################
    
    ### input parameters
    
    cellSubsets='ALL'
    cellCol_label_name="predictedGroup_Col2.5"
    returnAllPeaks=FALSE
    numCores=10


    ### load ArchR Project

    ArchRProj <- ArchR::loadArchRProject('/path/to/ArchRProject')

    peaks_by_sample <-scMACS::callPeaks_by_population(ArchRProj=ArchRProj,
        cellSubsets=cellSubsets,
        cellCol_label_name=cellCol_label_name,
        sampleCol_label_name='Sample',
        returnAllPeaks=returnAllPeaks,
        numCores=numCores)


## <a name="example2"></a> Tutorial-1: Sample-specific Peak-calling

To make peak calls on specific samples rather than by pooling cells across samples the user must call a different function in scMACS (as shown below), which breaks down the fragment files in a Cell X Sample matrix to call sample-specific peaks across all cell populations. 

### Load Library
   
    #Load scMACS and accompanying libraries
    library(scMACS)
    library(data.table)
    library(ArchR)
    library(GenomicRanges)
    library(plyranges)
    
### Load data and assign parameters
    
    ################################################################################
    ### To call peaks using scMACS
    ### the user must input 4 different
    ### parameters: 
    
    ###  @param cellSubsets vector of strings. Cell subsets for which to call peaks. 
    ###         Optional, 'ALL' defaults to all cell populations in metadata file. 
    
    ###  @param cellCol_label_name string indicating which column in the metadata 
    ###         file contains the cell population label

    ###  @param returnAllPeaks boolean. Indicates whether scMACS should return object 
    ###         containing all genomic regions or just the positive (+) called peaks. 
    
    ###  @param numCores integer. Number of cores to parallelize peak-calling 
    ###        across multiple cell populations 
    
    ################################################################################
    
    ### input parameters
    
    cellSubsets='ALL'
    cellCol_label_name="predictedGroup_Col2.5"
    returnAllPeaks=FALSE
    numCores=10


    ### load ArchR Project

    ArchRProj <- ArchR::loadArchRProject('/path/to/ArchRProject')

    ## call sample specific peaks
    peaks_by_sample <- scMACS::callPeaks_by_sample(ArchRProj=ArchRProj,
        cellSubsets=cellSubsets,
        cellCol_label_name=cellCol_label_name,
        sampleCol_label_name='Sample',
        returnAllPeaks=returnAllPeaks,
        numCores=numCores)

The result will be a nested list of lists where the hierarchy is organized as follows (using CD4 Naive & CD14 Monocytes as examples): 
   - CD4 Naive 
       - sample 1 peaks
       - sample 2 peaks 
        ...
       - sample N peaks 
        
   - CD14 Monocytes
       - sample 1 peaks 
       - sample 2 peaks
       ...
       - sample N peaks
 
# <a name="contact"></a> Contact

To contact the developers on issues and feature requests, please contact us via the discussions tab for feature requests, or open issues for any bugs. 
    
# <a name="license"></a> License

scMACS follows the Allen Institute Software License - full information about the license can be found on the LICENSE file. 
