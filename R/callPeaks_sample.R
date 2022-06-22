#' @title \code{callPeaks_by_sample}
#'
#' @description \code{callPeaks_by_sample} is the main peak-calling function in scMACS
#'              that serves as a wrapper function to call peaks provided a set 
#'              of fragment files and an ArchR Project for meta-data purposes
#'
#'
#' @param ArchRProj an ArchR Project 
#' @param metadf meta data file
#' @param cellType_Samples vector of strings. CellType samples to call peaks on.
#' @param metaColumn string indicating which column in the meta data file contains 
#'        the cell population label
#'
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
                                metadf,
                                cellType_Samples,
                                metaColumn,
                                returnAllPeaks=TRUE,
                                numCores=30,
                                returnFrags=F
                     
                     ){
    
    ## identify cell population to export 
    ## across samples by cell types
    
    cellType_Samples <- unique(cellType_Samples)
    frags <- getPopFrags(ArchRProj, 
                         metaColumn = metaColumn, 
                         cellSubsets = cellType_Samples, 
                         numCores= numCores)
    
    pryr::mem_used()
    sample_names <- names(frags)
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
    # ## different projects per downsample 

    barcodes <- lapply(frags, function (x) unique(x$RG)
               )
    
    subset_ArchR_projects <- lapply(barcodes, function(x) 
           subsetCells(ArchRProj, x)
           )

    pryr::mem_used()
                               
    ### call peaks using 
    ### scMACS functionalities
    peak_lists <- mclapply(1:length(subset_ArchR_projects),
                           function(x)

                    callPeaks_by_population(ArchRProj=subset_ArchR_projects[[x]],
                                            cellSubsets=sample_names[x],
                                            cellCol_label_name=metaColumn,
                                            returnAllPeaks=TRUE,
                                            numCores=1,
                                            totalFrags=normalization_factors[x],
                                            fragsList =frags[[x]],
                                            StudypreFactor = study_prefactor
                           ),
                          mc.cores=numCores

                        )
    
    rm(subset_ArchR_projects, frags)
    pryr::mem_used()
    
    return(peak_lists)

}
