#' @title \code{callPeaks_by_sample}
#'
#' @description \code{callPeaks_by_sample} is the main peak-calling function in scMACS
#'              that serves as a wrapper function to call peaks provided a set 
#'              of fragment files and an ArchR Project for meta-data purposes
#'
#'
#' @param ArchRProj an ArchR Project 
#' @param cellPopulations vector of strings. Cell subsets for which to call peaks. This list of group names must be identical to names that appear in the ArchRProject metadata.  Optional, if cellPopulations='ALL', then peak calling is done on all cell populations in the ArchR project metadata. Default is 'ALL'.
#' @param cellPopLabel string indicating which column in the ArchRProject metadata contains 
#'        the cell population label.
#' @param returnAllPeaks boolean. Indicates whether scMACS should return object containing all genomic regions or just the positive (+) called peaks. Default to the latter, only positive peaks. 
#' @param numCores integer. Number of cores to parallelize peak-calling across
#'                 multiple cell populations
#' @param returnFrags boolean. Determines if fragment files should be returned
#'
#' @return 
#'
#' @details The technical details of the algorithm are found in XX.
#'
#' @references XX
#'
#' @export

callPeaks_by_sample <- function(ArchRProj,
                                cellPopLabel,
                                cellPopulations = "ALL",
                                returnAllPeaks = TRUE,
                                numCores = 30,
                                returnFrags = F
                     ){
    
    # Get cell metadata and blacklisted regions from ArchR Project
    cellColData <- ArchR::getCellColData(ArchRProj)
    blackList <- ArchR::getBlacklist(ArchRProj)

    # Get frags grouped by cell population and sample
    # This will also validate the input cellPopulations
    frags <- getPopFrags(
        ArchRProj = ArchRProj,
        metaColumn = cellPopLabel,
        cellSubsets = cellPopulations,
        region = NULL,
        numCores = numCores,
        sampleSpecific = TRUE,
        NormMethod = "nfrags", 
        blackList = NULL,
        overlapList = 50
    )
    
    # Check for and remove celltype_samples for which there are no fragments.
    fragsNoNull <- frags[lengths(frags) != 0]
    
    # Rename frags by cell population
    renamedFrags <- lapply(
        1:length(fragsNoNull),
        function(y){
            # Split out celltype and sample from the name
            x <- fragsNoNull[y]
            print(names(x))
            celltype_sample <- names(x)
            splits <- unlist(str_split(celltype_sample, "#"))
            celltype <- splits[1]
            sample <- unlist(str_split(splits[2], "__"))[1]
            # Rename the fragments with just the sample
            names(x) <- sample
            # Return as a list named for celltype
            output <- list(x)
            names(output) <- celltype
            output
        }
    )
    # Group them by cell population
    renamedFrags <- unlist(renamedFrags, recursive=FALSE)
    splitFrags <- split(renamedFrags, f=names(renamedFrags))
    
    # Main loop over all cell populations
    experimentList <- list()
    for (cellPop in names(splitFrags)){
        
        print(str_interp("Calling peaks for cell population ${cellPop}"))
        
        # Get our fragments for this cellPop
        popFrags <- unlist(splitFrags[[cellPop]])
        
        # Simplify sample names to remove everything before the first "."
        sampleNames <-  gsub("^([^.]+).", "", names(popFrags))
        names(popFrags) <- sampleNames
        
        # Calculate normalization factors as the number of fragments for each celltype_samples
        normalization_factors <- as.integer(sapply(popFrags, length))

        # Add prefactor multiplier across datasets
        curr_frags_median <- median(cellColData$nFrags)
        study_prefactor = 3668/curr_frags_median # Training median

        # ###########################################################
        # ###########################################################
        # ## call scMACS peaks 

        # ## subset archr projects using 
        # ## the barcodes to create 
        # ## different projects per downsample 

    #     barcodes <- lapply(
    #         popFrags,
    #         function(x){
    #             unlist(unique(x$RG))
    #         }
    #     )

    #     subset_ArchR_projects <- lapply(barcodes, function(x) 
    #            ArchR::subsetCells(ArchRProj, x)
    #     )


        # This mclapply would parallelize over each sample within a celltype.
        # Each arrow is a sample so this is allowed
        # (Arrow files are locked - one access at a time)
        peaksGRangesList <- parallel::mclapply(
            1:length(popFrags),
            function(x){
                callPeaks_by_population(
                    # cellColData = cellColData,
                    blackList = blackList,
                    # cellPopulation = sampleNames[x],
                    # cellPopLabel = cellPopLabel,
                    returnAllPeaks = TRUE,
                    numCores = numCores,
                    totalFrags = normalization_factors[x],
                    fragsList = popFrags[[x]],
                    StudypreFactor = study_prefactor
                )
            },
            mc.cores = numCores
        )

        names(peaksGRangesList) <- sampleNames

        # Package rangeList into a RaggedExperiment
        ragExp <- RaggedExperiment::RaggedExperiment(
            peaksGRangesList
        )
        
        # And add it to the experimentList for this cell population
        experimentList <- append(experimentList, ragExp)
    }
    
    # TODO: Add experimentList to MultiAssayExperiment
    names(experimentList) <- names(splitFrags)
    browser()
    results <- MultiAssayExperiment::MultiAssayExperiment(
        experiments = experimentList
    )
    return(results)
}
