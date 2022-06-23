#' @title \code{callPeaks_by_sample}
#'
#' @description \code{callPeaks_by_sample} is the main peak-calling function in scMACS
#'              that serves as a wrapper function to call peaks provided a set 
#'              of fragment files and an ArchR Project for meta-data purposes
#'
#'
#' @param ArchRProj an ArchR Project 
#' @param cellPopulations vector of strings. Cell subsets for which to call peaks. Optional, if cellPopulations='ALL', then peak calling is done on all cell populations in the ArchR project metadata
#' @param cellPopLabel string indicating which column in the meta data file contains 
#'        the cell type label
#' @param sampleLabel string indicating which column in the meta data file contains 
#'        the samples
#' @param returnAllPeaks boolean. Indicates whether scMACS should return object containing all genomic regions or just the positive (+) called peaks. Default to the latter, only positive peaks. 
#'
#' @param numCores integer. Number of cores to parallelize peak-calling across
#'                 multiple cell populations
#' @param returnFrags boolean. Determines if fragment files should be returned
#'
#' @return scMACs_PeakList an list containing peak calls for each cell population passed on in the 
#'         cell subsets argument. Each peak call is returned as as Genomic Ranges object.
#'
#' @details The technical details of the algorithm are found in XX.
#'
#' @references XX
#'
#' @export

callPeaks_by_sample <- function(ArchRProj,
                                cellPopulations = "ALL",
                                cellPopLabel, 
                                sampleLabel,
                                returnAllPeaks=FALSE,
                                numCores=30
                     ){

    # Get fragments for our given cell subsets.
    cellPopulations <- unique(cellPopulations)
    frags <- getPopFrags(ArchRProj, 
                         metaColumn = cellPopLabel, 
                         cellSubsets = cellPopulations, 
                         numCores = numCores)
    
    cellPopulations <- names(frags)
    
    # Why not:
    # Get fragments specific to samples+celltype
    # callPeaks_by_population() for each sample+celltype
    ## Would avoid subsetting entire ArchR projects
    

    ## identify cell barcodes 
    normalization_factors <- as.numeric(sapply(frags, length))
    
    
    
    ## add preFactor multiplier across datasets
    curr_frags_median <- median(ArchRProj@cellColData$nFrags)
    study_prefactor = 3668/curr_frags_median # training median
    
    normalization_factors = normalization_factors 
    # ###########################################################
    # ###########################################################
    # ## call scMACS peaks 

    # ## subset archr projects using 
    # ## the barcodes to create 
    # ## different projects per population 

    barcodes <- lapply(frags, function (x) unique(x$RG)
               )

    subset_ArchR_projects <- parallel::mclapply(barcodes, function(x) 
           subsetCells(ArchRProj, x),
                                      mc.cores=numCores
           )

    #pryr::mem_used()

    ### call peaks using 
    ### scMACS functionalities
    peak_lists <- parallel::mclapply(1:length(subset_ArchR_projects),
                           function(x)
                           callPeaks_by_population(ArchRProj=subset_ArchR_projects[[x]],
                                            cellPopulations=cellPopulations[x],
                                            cellPopLabel=cellPopLabel,
                                            returnAllPeaks=returnAllPeaks,
                                            numCores=numCores,
                                            totalFrags=normalization_factors[x],
                                            fragsList =frags[[x]]

                           ),
                          mc.cores=numCores

                        )
    
    rm(subset_ArchR_projects, frags)
    
    return(peak_lists)

}
