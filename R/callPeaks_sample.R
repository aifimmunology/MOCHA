#' @title \code{callPeaks_by_sample}
#'
#' @description \code{callPeaks_by_sample} is the main peak-calling function in scMACS
#'              that serves as a wrapper function to call peaks provided a set 
#'              of fragment files and an ArchR Project for meta-data purposes
#'
#'
#' @param ArchRProj an ArchR Project 
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
                                cellType_Samples,
                                metaColumn,
                                returnAllPeaks=TRUE,
                                numCores=30,
                                returnFrags=F
                     
                     ){
    
    ## identify cell population to export 
    ## across samples by cell types
    
    cellType_Samples <- unique(cellType_Samples)
    
    # Get frags where the cellSubsets are pre-grouped cellType_Samples
    frags <- getPopFrags(ArchRProj, 
                         metaColumn = metaColumn, 
                         cellSubsets = cellType_Samples, 
                         numCores= numCores)
    # Check for and remove celltype_samples for which there are no fragments.
    frags_no_null <- frags[lengths(frags) != 0]
    
    
    sample_names <- names(frags_no_null)
    
    # calculate normalization factors as the number of fragments for each celltype_samples
    normalization_factors <- as.integer(sapply(frags_no_null, length))
    
    
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

    barcodes <- lapply(
        frags_no_null,
        function(x){
            unlist(unique(x$RG))
        }
    )
    
    subset_ArchR_projects <- lapply(barcodes, function(x) 
           ArchR::subsetCells(ArchRProj, x)
    )

                           
    ### call peaks using 
    ### scMACS functionalities
    peak_lists <- lapply(
        1:length(subset_ArchR_projects),
        function(x){
            callPeaks_by_population(
                ArchRProj=subset_ArchR_projects[[x]],
                cellSubsets=sample_names[x],
                cellCol_label_name=metaColumn,
                returnAllPeaks=TRUE,
                numCores=1,
                totalFrags=normalization_factors[x],
                fragsList=frags_no_null[[x]],
                StudypreFactor = study_prefactor
            )
        }#,
        #mc.cores=numCores
    )
    
    rm(subset_ArchR_projects, frags)
    
    
    return(peak_lists)

}
