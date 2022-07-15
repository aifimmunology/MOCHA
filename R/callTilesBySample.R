#' @title \code{callPeaks}
#'
#' @description \code{callPeaks} is the main peak-calling function in scMACS
#'              that serves as a wrapper function to call peaks provided a set 
#'              of fragment files and cell metadata
#'
#' @param blackList A GRanges object containing a blacklist of regions to exclude
#' @param returnAllPeaks boolean. Indicates whether scMACS should return object containing all genomic regions or just the positive (+) called peaks. Default to the latter, only positive peaks. 
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
#' @noRd
#' 

callTilesBySample <- function(
    blackList,
    returnAllTiles=FALSE,
    numCores=10,
    totalFrags,
    fragsList,
    StudypreFactor
){

    # Coefficients trained on ~ 3600 frags per cell 
    # Future datasets need to be calibrated to
    # these coefficients 
    finalModelObject = scMACS::finalModelObject
    thresholdModel = scMACS::thresholdModel

    FinalBins <-  scMACS:::determine_dynamic_range(
        AllFragmentsList = fragsList,
        blackList = blackList,
        binSize = 500, 
        doBin = FALSE
    )

    countsMatrix <- scMACS:::calculate_intensities(
        fragMat = fragsList,
        candidatePeaks = FinalBins,
        totalFrags = totalFrags
    )

    countsMatrix$TotalIntensity <- countsMatrix$TotalIntensity * StudypreFactor
    countsMatrix$maxIntensity <- countsMatrix$maxIntensity * StudypreFactor    

    countsMatrix = countsMatrix[countsMatrix$TotalIntensity > 0,]
    scMACS_tiles <- scMACS:::make_prediction(
        X = countsMatrix,
        finalModelObject = finalModelObject
    )

    if(!returnAllTiles){
        scMACS_tiles <- scMACS_tiles[scMACS_tiles$peak==T]
    }
    
    return(scMACS_tiles)

}
