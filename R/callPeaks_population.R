#' @title \code{callPeaks}
#'
#' @description \code{callPeaks} is the main peak-calling function in scMACS
#'              that serves as a wrapper function to call peaks provided a set 
#'              of fragment files and cell metadata
#'
#'
#' @param cellColData The cell-level metadata from the original ArchR Project
#' @param blackList A GRanges object containing a blacklist of regions to exclude
#' @param cellPopulation String/character vector. Cell subset for which to call peaks.
#' @param cellCol_label_name string indicating which column in the meta data file contains 
#'        the cell population label
#' 
#' @param returnAllPeaks boolean. Indicates whether scMACS should return object containing all genomic regions or just the positive (+) called peaks. Default to the latter, only positive peaks. 
#' 
#' @param numCores integer. Number of cores to parallelize peak-calling across
#'                 multiple cell populations 
#'
#' @param totalFrags # of fragments in that sample for that cell population
#' @param StudypreFactor ratio of average signal between new study and training set
#' @return scMACs_PeakList an list containing peak calls for each cell population passed on in the 
#'         cell subsets argument. Each peak call is returned as as Genomic Ranges object.
#' 
#' @details The technical details of the algorithm are found in XX.
#'
#' @references XX
#'
#' @export

callPeaks_by_population <- function(
    cellColData,
    blackList,
    cellPopulation,
    cellCol_label_name,
    returnAllPeaks=FALSE,
    numCores=10,
    totalFrags,
    fragsList,
    StudypreFactor
){

    ### User-input/Parameter Checks
    
      if(class(cellPopulation)!='character'){
        stop('cellPopulation must be a string indicating cell population')
      }
      if(class(cellCol_label_name)!='character'){
        stop('cellCol_label_name must be a string indicating the column in cellColData containing cell population labels')      
      }   
      if(is.null(fragsList)){
          stop('Load fragments prior to running scMACS')
          
      }

    ## coefficients trained on ~ 3600 frags per cell 
    ## and future datasets need to be calibrated to
    ## these coefficients 
    finalModelObject = scMACS::finalModelObject
    thresholdModel = scMACS::thresholdModel
    
    # Rename parameters for downstream use TODO remove in later final cleanup 
    cellPopulation <- cellPopulation
    
    ### get barcnodes by cell pop for 
    ### peak-calling by different 
    ### cell population 
    
    cellNames <- row.names(cellColData)[which(cellColData[,cellCol_label_name]==cellPopulation)]

    cellsPerPop <- length(cellNames)
    
    df <- cbind(cellPopulation, cellsPerPop)
        
    cat('\n\nRunning peak calls for the following cell population:\n')
    print(df)
    
    #########################
    #########################
    ######################### 

    FinalBins <-  determine_dynamic_range(
        AllFragmentsList = fragsList,
        blackList = blackList,
        binSize = 500, 
        doBin = FALSE
    )

    countsMatrix <- calculate_intensities(
        fragMat = fragsList,
        candidatePeaks = FinalBins,
        totalFrags = totalFrags
    )

    countsMatrix$TotalIntensity <- countsMatrix$TotalIntensity * StudypreFactor
    countsMatrix$maxIntensity <- countsMatrix$maxIntensity * StudypreFactor    

    countsMatrix = countsMatrix[countsMatrix$TotalIntensity > 0,]
    scMACS_peaks <- make_prediction(
        X = countsMatrix,
        finalModelObject = finalModelObject
    )

    if(!returnAllPeaks){
        scMACS_peaks <- scMACS_peaks[scMACS_peaks$peak==T]
    }
    
    return(scMACS_peaks)

}
