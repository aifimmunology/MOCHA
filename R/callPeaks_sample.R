#' @title \code{callPeaks_by_sample}
#'
#' @description \code{callPeaks_by_sample} is the main peak-calling function in scMACS
#'              that serves as a wrapper function to call peaks provided a set 
#'              of fragment files and an ArchR Project for meta-data purposes
#'
#'
#' @param ArchRProj an ArchR Project 
#' @param metadf TODO 
#' @param cellSubsets vector of strings. Cell subsets for which to call peaks. Optional, if cellSubsets='ALL', then peak calling is done on all cell populations in the ArchR project metadata
#' @param cellCol_label_name string indicating which column in the meta data file contains 
#'        the cell population label
#' @param cellType_sample_label_name string indicating which column in the meta data file contains 
#'        the samples
#'
#' @param returnAllPeaks boolean. Indicates whether scMACS should return object containing all genomic regions or just the positive (+) called peaks. Default to the latter, only positive peaks. 
#'
#' @param numCores integer. Number of cores to parallelize peak-calling across
#'                 multiple cell populations
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
                                cellSubsets,
                                cellCol_label_name, 
                                cellType_sample_label_name,
                                returnAllPeaks=TRUE,
                                numCores=30
                     ){
    
    ## identify cell population to export 
    ## across samples by cell types
    
    cellTypesToExport <- unique(cellSubsets)
    frags <- getPopFrags(ArchRProj, 
                         metaColumn = cellCol_label_name, 
                         cellSubsets = cellSubsets, 
                         numCores= numCores)
    
    # Why not:
    # Get fragments specific to samples+celltype
    # callPeaks_by_population() for each sample+celltype
    ## Would avoid subsetting entire ArchR projects
    
    names(frags) <- cellSubsets

    ## identify cell barcodes 
    normalization_factors <- sapply(frags, length)

    # ###########################################################
    # ###########################################################
    # ## call scMACS peaks 

    # ## subset archr projects using 
    # ## the barcodes to create 
    # ## different projects per downsample 

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
                                            cellSubsets=cellSubsets[x],
                                            cellCol_label_name=cellType_sample_label_name,
                                            returnAllPeaks=TRUE,
                                            numCores=1,
                                            totalFrags=normalization_factors[x],
                                            fragsList =frags[[x]]
                           ),
                          mc.cores=numCores

                        )
    
    finalObject = list(Sample_peaks=peak_lists,
                       Sample_fragments=frags)
    
    return(finalObject)
    

}
