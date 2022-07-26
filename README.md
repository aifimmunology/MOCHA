README for scMACS
===============

* * *

Table of Contents
-----------------

* [Introduction](#introduction)
* [Install package and load library](#library)
* [Package Vignette on COVID PASC dataset](#vignette)
* [Tutorials](#example-main)
    * [Tutorial-1: Peak-calling on Entire Cell Populations](#example1)
    * [Tutorial-2: Sample-specific Peak-calling](#example2)
    * [Tutorial-3: Differential Accessibility](#example3)
    * [Tutorial-4: Sample X Peak Matrix](#example4)


* [Contact](#contact)
* [License](#license)

* * *

# <a name="introduction"></a> Introduction
scMACS is an R package containing a novel single-cell peak-calling algorithm that leverages single-cell information to determine whether a particular genomic region is open by calculating two measures of intensities, and using these to call peaks via a hierarchical model. 

<br> ![img](vignettes/scMACS_overview.png) <br>
scATAC processing and overview of the **scMACS** algorithm. 

# <a name="library"></a> Install package and load library

To install library, run the following lines of code to install directly from GitHub

    devtools::install_github("aifimmunology/scMACS", ref="package-details")
    library("scMACS")

# <a name="vignette"></a> Usage
Example usage can be found in this [vignette](vignettes/COVID_example.R).


# <a name="example-main"></a> Tutorials - deprecated
## <a name="example1"></a> Tutorial-1: Peak-calling on Entire Cell Populations

This tutorial demonstrates how to call peaks for a given cell population. This process
assumes pooling cells across samples/subjects to call peaks for all cells in a given 
cell population.

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
                        PeakID seqnames    start      end strand     lambda1
     1:     chr1:816000-816499     chr1   816000   816499      * 0.001513241
     2:     chr1:817000-817499     chr1   817000   817499      * 0.030110403
     3:     chr1:819500-819999     chr1   819500   819999      * 0.001111769
     4:     chr1:820000-820499     chr1   820000   820499      * 0.002130890
     5:     chr1:820500-820999     chr1   820500   820999      * 0.001729418
    ---                                                                     
    701219: chrY:21647500-21647999     chrY 21647500 21647999      * 0.001327946
    701220: chrY:21648000-21648499     chrY 21648000 21648499      * 0.001080886
    701221: chrY:21650500-21650999     chrY 21650500 21650999      * 0.005064724
    701222: chrY:21651000-21651499     chrY 21651000 21651499      * 0.001019121
    701223: chrY:26315000-26315499     chrY 26315000 26315499      * 0.002717657


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
        

The trimmed output (below) is a nested list of lists (such as the output above) where the hierarchy is organized as follows (CD14 Monocyte peak calls as examples), with the last entry containing the union of all sample-specific peaks. 
        
   - CD14 Monocytes
       - sample 1 peaks 
       - sample 2 peaks
       - ...
       - sample N peaks

    > peaks_by_sample['CD14 Mono']
    
    `CD14 Mono`|`X002_PB2216W7-01`
                            PeakID seqnames    start      end strand     lambda1
         1:     chr1:817000-817499     chr1   817000   817499      * 0.034305124
         2:     chr1:821000-821499     chr1   821000   821499      * 0.003026923
         3:     chr1:821500-821999     chr1   821500   821999      * 0.003026923
         4:     chr1:822000-822499     chr1   822000   822499      * 0.004035897
         5:     chr1:826500-826999     chr1   826500   826999      * 0.006053845
        ---                                                                     
    386614: chrY:21650500-21650999     chrY 21650500 21650999      * 0.006053845
    386615: chrY:21681500-21681999     chrY 21681500 21681999      * 0.003026923
    386616: chrY:26315000-26315499     chrY 26315000 26315499      * 0.004035897
    386617: chrY:26400500-26400999     chrY 26400500 26400999      * 0.004035897
    386618: chrY:26441000-26441499     chrY 26441000 26441499      * 0.003026923
      
    `CD14 Mono`|`X001_PB7626W7-01`
                        PeakID seqnames    start      end strand     lambda1
     1:     chr1:817000-817499     chr1   817000   817499      * 0.017246420
     2:     chr1:823000-823499     chr1   823000   823499      * 0.002555025
     3:     chr1:827000-827499     chr1   827000   827499      * 0.024911496
     4:     chr1:827500-827999     chr1   827500   827999      * 0.037686622
     5:     chr1:869500-869999     chr1   869500   869999      * 0.030660302
    ---                                                                     
    213064: chrY:21260500-21260999     chrY 21260500 21260999      * 0.003193782
    213065: chrY:21316000-21316499     chrY 21316000 21316499      * 0.003193782
    213066: chrY:21422000-21422499     chrY 21422000 21422499      * 0.016607664
    213067: chrY:21440500-21440999     chrY 21440500 21440999      * 0.003193782
    213068: chrY:21650500-21650999     chrY 21650500 21650999      * 0.004471294
     
    `CD14 Mono`|Union
                     PeakID seqnames    start      end strand
    1:     chr1:817000-817499     chr1   817000   817499      *
    2:     chr1:821000-821499     chr1   821000   821499      *
    3:     chr1:826500-826999     chr1   826500   826999      *
    4:     chr1:827000-827499     chr1   827000   827499      *
    5:     chr1:827500-827999     chr1   827500   827999      *
    ---                                                         
    1856185: chrY:14067000-14067499     chrY 14067000 14067499      *
    1856186: chrY:14370500-14370999     chrY 14370500 14370999      *
    1856187: chrY:16644000-16644499     chrY 16644000 16644499      *
    1856188: chrY:19824000-19824499     chrY 19824000 19824499      *
    1856189: chrY:20655500-20655999     chrY 20655500 20655999      *

## <a name="example3"></a> Tutorial-3: Differential Accesibility

To calculate differential accessibility, the code below can be invoked to 
determine whether a particular region is open or not.

    ## two different cell pop'ns with 
    ## enough samples to determine
    ## differential accessibility

    ## define groups for differential 
    ## analyses 
    groupA = peaks_by_sample[['CD8 TCM']]
    groupB = peaks_by_sample[['B memory']]

    ## define a sample of peaks to look at 
    ## for differential accessibility 
    candidatePeaks <- sample(allPeaks,1000)

    ## run function call 
    tmp <- differential_accessibility(groupA, groupB, candidatePeaks=candidatePeaks, doPDRAnalysis=FALSE)
    
    ## examine output
    > head(tmp, 10) 
  
                            Peak ES_wilc   Wilcoxon     ES_Chisq Chisquare
    W  chr14:100179000-100179499    59.5 0.59298010 7.149237e-32 1.0000000
    W1    chr1:28838000-28838499    59.5 0.59298010 7.149237e-32 1.0000000
    W2   chr22:21102500-21102999    45.0 0.23180711 4.901620e-01 0.4838550
    W3    chr4:95244000-95244499    59.5 0.59298010 7.149237e-32 1.0000000
    W4  chr1:206725500-206725999    38.5 0.06666009 1.137894e+00 0.2860978
    W5   chr15:99754500-99754999    59.5 0.59298010 7.149237e-32 1.0000000
    W6   chr15:60893000-60893499    59.5 0.59298010 7.149237e-32 1.0000000
    W7      chr1:3567000-3567499    56.0 0.40326191 9.642857e-04 0.9752273
    W8  chr1:175029000-175029499    52.5 0.28379398 1.244444e-01 0.7242632
    W9  chr4:113760500-113760999    24.5 0.01133921 2.244801e+00 0.1340642
          MinPval L1A_avg L1B_avg
    W  0.59298010    0.00  0.0007
    W1 0.59298010    0.00  0.0007
    W2 0.23180711    0.02  0.0058
    W3 0.59298010    0.00  0.0007
    W4 0.06666009    0.00  0.0066
    W5 0.59298010    0.00  0.0007
    W6 0.59298010    0.00  0.0007
    W7 0.40326191    0.00  0.0012
    W8 0.28379398    0.00  0.0015
    W9 0.01133921    0.00  0.0106

where the output contains the following pieces of information quantifying the differential accessibility:

- Peak= Peak ID
- ES_wilc = effect size of the wilcoxon rank sum test
- Wilcoxon = the statistical significance of the Rank sum test
- ES_Chisq = effect size of the chi square test
- Chisquare = tthe statistical significance of the chi-square test
- MinPval = minimum p-value across both tests
- L1A_Avg= Avg lambda1 value across samples in group A
- L1B_Avg= Avg lambda1 value across samples in group B

to qualify whether that region is differential accessible or not. 

In addition, an optional parameter (not shown above) allows the user to calcualte PDR (probability of detection rate) he probability of detecting a read (PDR) is a new metric we developed to qualify technical noise across samples. The PDR ( probability detection rate) is a technical measure defined to capture dropout rate in single-cell ATAC data to better qualify whether differential accessibility is a product of technical noise or true biological signal. It calculates this metric, for each group, by comparing the most open sample (highest lambda1) to the smallest sample (smallest N), by randomly downsampling cells of the most open sample (smallest N) and seeing how often was biological signal detected in that subsample. 

## <a name="example4"></a> Tutorial-4: Sample X Peak Matrix 

In order to summarize peak-calls for a given population, the code 
below can be used to create a sample X peak matrix that includes 
the union of all peaks (across samples), and includes the lambda1 measurement
for each peak in which reads were observed across all patients. 

    # Here we generate the 
    # peak-sample matrix 
    # getting metrics for lambda1
    # for each sample 

    cd14_samples <- peaks_by_sample[['CD14 Mono']]

    sample_specific_matrix <- create_peak_sampleMatrix(cd14_samples)
    
    ## trimmed output 
    > sample_specific_matrix
                                PeakID X001_PB5206W2-01 X001_PB5206W4-01
          1:   chr10:10000000-10000499      0.000000000                0
          2: chr10:100000000-100000499      0.000000000                0
          3: chr10:100001000-100001499      0.000000000                0
          4: chr10:100001500-100001999      0.000000000                0
          5: chr10:100002000-100002499      0.003679851                0
         ---                                                            
    1474125:      chrY:9546500-9546999      0.000000000                0
    1474126:      chrY:9811000-9811499      0.000000000                0
    1474127:      chrY:9956000-9956499      0.000000000                0
    1474128:      chrY:9974500-9974999      0.000000000                0
    1474129:      chrY:9993000-9993499      0.000000000                0
    
    ...
    
                X002_PB2216W4-01 X002_PB2216W5-01 X002_PB2216W6-01 X002_PB2216W7-01
          1:                0                0                0      0.002577198
          2:                0                0                0      0.000000000
          3:                0                0                0      0.000000000
          4:                0                0                0      0.000000000
          5:                0                0                0      0.000000000
         ---                                                                    
    1474125:                0                0                0      0.000000000
    1474126:                0                0                0      0.000000000
    1474127:                0                0                0      0.002577198
    1474128:                0                0                0      0.000000000
    1474129:                0                0                0      0.000000000


# <a name="contact"></a> Contact

To contact the developers on issues and feature requests, please contact us via the discussions tab for feature requests, or open issues for any bugs. 
    
# <a name="license"></a> License

scMACS follows the Allen Institute Software License - full information about the license can be found on the LICENSE file. 
