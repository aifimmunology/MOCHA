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
    * [Tutorial-4: Peak Widths](#example4)

* [Contact](#contact)
* [License](#license)

* * *

# <a name="introduction"></a> Introduction
scMACS is an R package containing a novel single-cell peak-calling algorithm that leverages single-cell information to determine whether a particular genomic region is open by calculating two measures of intensities, and using these to call peaks via a hierarchical model. 

<br> ![img](vignettes/scMACS_overview.png) <br>
scATAC processing and overview of the **scMACS** algorithm. 

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

    ### Call peaks across cell populations
    peaks_by_population <-scMACS::callPeaks_by_population(ArchRProj=ArchRProj,
        cellSubsets=cellSubsets,
        cellCol_label_name=cellCol_label_name,
        returnAllPeaks=returnAllPeaks,
        numCores=numCores)
    
    ### Extract peaks for CD14 Cells
    > peaks_by_population['CD14 Mono']

    `CD14 Mono`
    GRanges object with 701223 ranges and 8 metadata columns:
           seqnames            ranges strand |    lambda1     lambda2
              <Rle>         <IRanges>  <Rle> |  <numeric>   <numeric>
       [1]     chr1     816000-816499      * | 0.00151324 0.000346078
       [2]     chr1     817000-817499      * | 0.03011040 0.000840476
       [3]     chr1     819500-819999      * | 0.00111177 0.000840476
       [4]     chr1     820000-820499      * | 0.00213089 0.000840476
       [5]     chr1     820500-820999      * | 0.00172942 0.000840476
       ...      ...               ...    ... .        ...         ...
    [701219]     chrY 21647500-21647999      * | 0.00132795 0.000346078
    [701220]     chrY 21648000-21648499      * | 0.00108089 0.000346078
    [701221]     chrY 21650500-21650999      * | 0.00506472 0.000593277
    [701222]     chrY 21651000-21651499      * | 0.00101912 0.000346078
    [701223]     chrY 26315000-26315499      * | 0.00271766 0.000593277
               lambda3 maxIntensity  numCells Prediction PredictionStrength
             <numeric>    <integer> <integer>  <numeric>          <numeric>
       [1]  0.00432257            2     16009   0.916841         0.00151324
       [2]  0.00469423            4     16009   1.000000         0.03011040
       [3]  0.00515960            4     16009   0.640429         0.00111177
       [4]  0.00523143            4     16009   0.994498         0.00213089
       [5]  0.00522831            4     16009   0.966985         0.00172942
       ...         ...          ...       ...        ...                ...
    [701219] 0.001061903            2     16009   0.826380         0.00132795
    [701220] 0.001061903            2     16009   0.608304         0.00108089
    [701221] 0.001077519            3     16009   1.000000         0.00506472
    [701222] 0.001071272            2     16009   0.539962         0.00101912
    [701223] 0.000480979            3     16009   0.999614         0.00271766
                Peak
           <logical>
       [1]      TRUE
       [2]      TRUE
       [3]      TRUE
       [4]      TRUE
       [5]      TRUE
       ...       ...
    [701219]      TRUE
    [701220]      TRUE
    [701221]      TRUE
    [701222]      TRUE
    [701223]      TRUE
    -------
    seqinfo: 24 sequences from an unspecified genome; no seqlengths
  
  

## <a name="example2"></a> Tutorial-2: Sample-specific Peak-calling

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
        

The trimmed output (below) is a nested list of lists (such as the output above) where the hierarchy is organized as follows (CD14 Monocyte peak calls as examples): 
        
   - CD14 Monocytes
       - sample 1 peaks 
       - sample 2 peaks
       - ...
       - sample N peaks

    > peaks_by_sample['CD14 Mono']
    
    `CD14 Mono`|`X001_PB5206W4-01`
    GRanges object with 284164 ranges and 8 metadata columns:
               seqnames              ranges strand |    lambda1     lambda2
                  <Rle>           <IRanges>  <Rle> |  <numeric>   <numeric>
           [1]     chr1       817000-817499      * | 0.02124364 0.001235993
           [2]     chr1       827000-827499      * | 0.02220926 0.001928150
           [3]     chr1       827500-827999      * | 0.04377477 0.001928150
           [4]     chr1       869500-869999      * | 0.03122171 0.001235993
           [5]     chr1       870000-870499      * | 0.00579372 0.000593277
           ...      ...                 ...    ... .        ...         ...
      [284160]     chrX 155891000-155891499      * | 0.00193124 0.001235993
      [284161]     chrX 155891500-155891999      * | 0.00257499 0.000593277
      [284162]     chrX 155892000-155892499      * | 0.00160937 0.000593277
      [284163]     chrX 155893500-155893999      * | 0.00160937 0.000593277
      [284164]     chrX 155898500-155898999      * | 0.00193124 0.001235993
      -------
      seqinfo: 24 sequences from an unspecified genome; no seqlengths
      
    `CD14 Mono`|`X001_PB5206W7-01`
    GRanges object with 274422 ranges and 8 metadata columns:
               seqnames              ranges strand |    lambda1     lambda2
                  <Rle>           <IRanges>  <Rle> |  <numeric>   <numeric>
           [1]     chr1       817000-817499      * | 0.01837908 0.001582072
           [2]     chr1       820000-820499      * | 0.00183791 0.000741596
           [3]     chr1       821500-821999      * | 0.00183791 0.000741596
           [4]     chr1       822000-822499      * | 0.00275686 0.001582072
           [5]     chr1       822500-822999      * | 0.00229739 0.001582072
           ...      ...                 ...    ... .        ...         ...
      [274418]     chrX 155885500-155885999      * | 0.00183791 0.000741596
      [274419]     chrX 155886000-155886499      * | 0.00275686 0.001582072
      [274420]     chrX 155891500-155891999      * | 0.00229739 0.000741596
      [274421]     chrX 155893500-155893999      * | 0.00229739 0.000741596
      [274422]     chrX 155896500-155896999      * | 0.00275686 0.000741596

      -------
      seqinfo: 24 sequences from an unspecified genome; no seqlengths

    `CD14 Mono`|`X001_PB7626W2-01`
    GRanges object with 211044 ranges and 8 metadata columns:
               seqnames            ranges strand |    lambda1     lambda2
                  <Rle>         <IRanges>  <Rle> |  <numeric>   <numeric>
           [1]     chr1     817000-817499      * | 0.01121082 0.000889915
           [2]     chr1     826500-826999      * | 0.00336325 0.000889915
           [3]     chr1     827000-827499      * | 0.02130057 0.001829270
           [4]     chr1     827500-827999      * | 0.03755626 0.001829270
           [5]     chr1     869500-869999      * | 0.02242165 0.000889915
           ...      ...               ...    ... .        ...         ...
      [211040]     chrY 21422000-21422499      * | 0.01121082 0.001829270
      [211041]     chrY 21440500-21440999      * | 0.00336325 0.000889915
      [211042]     chrY 21650500-21650999      * | 0.00280271 0.000889915
      [211043]     chrY 26315000-26315499      * | 0.00224216 0.000889915
      [211044]     chrY 26434500-26434999      * | 0.00224216 0.000889915

      -------
      seqinfo: 24 sequences from an unspecified genome; no seqlengths

    `CD14 Mono`|`X001_PB7626W4-01`
    GRanges object with 179328 ranges and 8 metadata columns:
               seqnames            ranges strand |    lambda1     lambda2
                  <Rle>         <IRanges>  <Rle> |  <numeric>   <numeric>
           [1]     chr1     817000-817499      * | 0.01622039  0.00202703
           [2]     chr1     827000-827499      * | 0.02205973  0.00202703
           [3]     chr1     827500-827999      * | 0.02725025  0.00202703
           [4]     chr1     869500-869999      * | 0.03568485  0.00202703
           [5]     chr1     870000-870499      * | 0.00583934  0.00202703
           ...      ...               ...    ... .        ...         ...
      [179324]     chrY 21260500-21260999      * | 0.00583934 0.000988795
      [179325]     chrY 21289000-21289499      * | 0.00454171 0.000988795
      [179326]     chrY 21315500-21315999      * | 0.00259526 0.000988795
      [179327]     chrY 21422000-21422499      * | 0.00778579 0.000988795
      [179328]     chrY 21650500-21650999      * | 0.00454171 0.000988795

      -------
      seqinfo: 24 sequences from an unspecified genome; no seqlengths


## <a name="example3"></a> Tutorial-3: Differential Accesibility

This part is currently in development. TBD 
 
# <a name="contact"></a> Contact

To contact the developers on issues and feature requests, please contact us via the discussions tab for feature requests, or open issues for any bugs. 
    
# <a name="license"></a> License

scMACS follows the Allen Institute Software License - full information about the license can be found on the LICENSE file. 
