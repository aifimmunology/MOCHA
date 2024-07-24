#' @title Get coverage for a given region
#'
#' @description \code{localFootprint} will extract the coverage files created by
#'   callOpenTiles and return a specific region's coverage
#'
#' @param SampleTileObj The SummarizedExperiment object output from
#'   getSampleTileMatrix
#' @param type Boolean. Default is true, and exports Coverage. If set to FALSE,
#'   exports Insertions.
#' @param region a GRanges object or vector or strings containing the regions of
#'   interest. Strings must be in the format "chr:start-end", e.g.
#'   "chr4:1300-2222".
#' @param cellPopulations vector of strings. Cell subsets for which to call
#'   peaks. This list of group names must be identical to names that appear in
#'   the SampleTileObj.  Optional, if cellPopulations='ALL', then peak calling
#'   is done on all cell populations. Default is 'ALL'.
#' @param groupColumn Optional, the column containing sample group labels for
#'   returning coverage within sample groups. Default is NULL, all samples will
#'   be used.
#' @param subGroups a list of subgroup(s) within the groupColumn from the
#'   metadata. Optional, default is NULL, all labels within groupColumn will be
#'   used.
#' @param sampleSpecific If TRUE, get a sample-specific count dataframe out.
#'   Default is FALSE, average across samples and get a dataframe out.
#' @param approxLimit Optional limit to region size, where if region is larger
#'   than approxLimit basepairs, binning will be used. Default is 100000.
#' @param binSize Optional numeric, size of bins in basepairs when binning is
#'   used. Default is 250.
#' @param sliding Optional numeric. Default is NULL. This number is the size of
#'   the sliding window for generating average intensities.
#' @param numCores integer. Number of cores to parallelize peak-calling across
#'   multiple cell populations
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#'
#' @return countSE a SummarizedExperiment containing coverage for the given
#'   input cell populations.
#'
#' @examples
#' \dontrun{
#' countSE <- MOCHA::localFootprint(
#'   SampleTileObj = SampleTileMatrices,
#'   cellPopulations = "ALL",
#'   region = "chr1:18137866-38139912",
#'   numCores = 30,
#'   sampleSpecific = FALSE
#' )
#' }
#'
#' @export
#' @keywords utils

localFootprint <- function(SampleTileObj,
                          type = TRUE,
                          region,
                          cellPopulations = "ALL",
                          groupColumn = NULL,
                          subGroups = NULL,
                          sampleSpecific = FALSE,
                          normTn5 = TRUE,
                          sumWidth = 10, 
                          medianWidth = 11,
                          numCores = 1,
                          verbose = FALSE) {
  . <- idx <- score <- NULL

  cellNames <- names(SummarizedExperiment::assays(SampleTileObj))
  metaFile <- SummarizedExperiment::colData(SampleTileObj)
  outDir <- SampleTileObj@metadata$Directory
    
  #Extract genome database
  genome_db = SampleTileObj@metadata$Genome
  outDir = SampleTileObj@metadata$Directory

  genome <- getAnnotationDbFromInstalledPkgname(dbName = genome_db, type = 'BSgenome')

  ## Check wthether the sample tile boject has an accurate directory with files, and insertion bias pre-calculated. 
  if (is.na(outDir)) {
    stop("Missing coverage file directory. SampleTileObj$metadata must contain 'Directory'.")
  }

  if (!file.exists(outDir)) {
    stop("Directory given by SampleTileObj@metadata$Directory does not exist. No insertions files can be loaded.")
  }
    
  if (normTn5 & !any(grepl('InsertionBias', names(SampleTileObj@metadata)))) {
    stop("InsertionBias has not been calculated. Please run addInsertionBias() on the SampleTileObj and try again.")
  }else if(normTn5){
      
      insertBias = SampleTileObj@metadata[['InsertionBias']][,'Norm']
      
  }else{
      
      insertionBias = NULL
      
  }

    #### Now process the region into a GRanges object.
  if (is.character(region)) {
    # Convert to GRanges
    regionGRanges <- MOCHA::StringsToGRanges(region)
  } else if (class(region)[1] == "GRanges") {
    regionGRanges <- region
  } else {
    stop("Wrong region input type. Input must either be a string, or a GRanges location.")
  }
    
  #Generate filter window, which is expanded for the rolling sum. 
   filterWindow <- plyranges::stretch(plyranges::anchor_center(regionGRanges), 
                            extend = 2*max(sumWidth, medianWidth)) 
    
   # If the region is too large, throw an error. 
   size1 = GenomicRanges::end(regionGRanges)- GenomicRanges::start(regionGRanges)
   if(size1 > 3000){
    
       stop('Please provide a region smaller than 3000 bps. Local footprint cannot realistically be identified at this scale')
       
   }
    
    
  if (all(toupper(cellNames) == "COUNTS")) {
    stop(
      "The only assay in the SummarizedExperiment is Counts. The names of assays must reflect cell types,",
      " such as those in the Summarized Experiment output of getSampleTileMatrix."
    )
  }

  if (all(toupper(cellPopulations) == "ALL")) {
    cellPopulations <- cellNames
  }
  if (!all(cellPopulations %in% cellNames)) {
    stop("Some or all cell populations provided are not found.")
  }

  # Pull out a list of samples by group.

  if (!is.null(subGroups) & !is.null(groupColumn)) {
    # If the user defined a list of subgroup(s) within the groupColumn from the metadata, then it subsets to just those samples
    subSamples <- lapply(subGroups, function(x) metaFile[metaFile[, groupColumn] %in% x, "Sample"])
    names(subSamples) <- subGroups
  } else if (!is.null(groupColumn)) {

    # If no subGroup defined, then it'll form a list of samples across all labels within the groupColumn
    subGroups <- unique(metaFile[, groupColumn])

    subSamples <- lapply(subGroups, function(x) metaFile[metaFile[, groupColumn] %in% x, "Sample"])
    names(subSamples) <- subGroups
                         
  } else {

    # If neither groupColumn nor subGroup is defined, then it forms one list of all sample names
    subGroups <- "All"
    subSamples <- list("All" = metaFile[, "Sample"])
      
  }

  # Determine if binning is needed to simplify things


  cl <- parallel::makeCluster(numCores)
  # Pull up the cell types of interest, and filter for samples and subset down to region of interest
  cellPopulation_Files <- lapply(cellPopulations, function(x) {
    if (verbose) {
      message(stringr::str_interp("Extracting insertions from cell population '${x}'"))
    }
      
    # Pull up insertions files
    originalCovGRanges <- readRDS(paste(outDir, "/", x, "_CoverageFiles.RDS", sep = ""))
    if ('Insertions' %in% names(originalCovGRanges)){
      originalInsertions <- originalCovGRanges[['Insertions']]
      #originalCoverage <- originalCovGRanges[['Insertions']]
    }else{
      stop('Error around reading insertions files. Check that coverage/insertions files are not corrupted.')
    }
      
    # Edge case: One or more samples are missing coverage for this cell population,
    # e.g. if a cell population only exists in one sample.
    for(y in seq_along(subSamples)) {
      if (!all(subSamples[[y]] %in% names(originalInsertions))) {
        missingSamples <- paste(subSamples[[y]][!subSamples[[y]] %in% 
                                        names(originalInsertions)], collapse = ", ")
          
        if(!force){
            
            stop(stringr::str_interp(c(
              "There is no insertion coverage for cell population '${x}' in the ",
              "following samples in sample grouping '${names(subSamples)[y]}': ",
              "${missingSamples}"
            )))
            
        }else{
            
            subSamples[[y]] = subSamples[[y]][subSamples[[y]] %in% 
                                              names(originalInsertions)]
            
        }
      }
    }
      
    ##Subset coverage down to only the relevant samples.
    originalInsertions =  originalInsertions[unlist(subSamples)]
    
    ## filter down to only insertions within the key region. 
    originalInsertions = lapply(originalInsertions, function(XX){
                          plyranges::join_overlap_intersect(XX, filterWindow)
        })
    
    names(originalInsertions) = unlist(subSamples)
   
    iterList <- lapply(originalInsertions, function(y) {
        list(y, normTn5, insertBias, genome_db, sumWidth, medianWidth)
        })

    normInsert <- pbapply::pblapply(cl = cl, X = iterList, normLocalInserts)
    
      
    cellPopFootprint <- do.call('rbind',lapply(names(normInsert), function(ZZ){
                                 tmpMat = dplyr::mutate(as.data.frame(normInsert[[ZZ]]), Samples = ZZ)[c(1,2,4,5)]
                                 colnames(tmpMat) = c('chr', 'Locus', 'Counts', 'Groups')
                                 tmpMat
                            }))

    # If not sample specific, take the average coverage across samples.
    # if it is sample specific, just subset down the coverage to the region of interest.
    if (!sampleSpecific) {       
        for(ZZ in names(subSamples)){
            cellPopFootprint[cellPopFootprint$Groups %in% subSamples[[ZZ]],'Groups'] = ZZ
        }
        cellPopFootprint <- dplyr::summarize(dplyr::group_by(cellPopFootprint, chr, Locus, Groups), Counts = mean(Counts))
        cellPopFootprint = cellPopFootprint[,c('chr', 'Locus', 'Counts', 'Groups')]
    } 
    gc()
    cellPopFootprint
    
  })

  names(cellPopulation_Files) <- cellPopulations

  parallel::stopCluster(cl)

  newMetadata <- SampleTileObj@metadata
  newMetadata$History <- append(newMetadata$History, paste("localFootprint", utils::packageVersion("MOCHA")))
 
  newMetadata$Type <- 'InsertionFootprint'

  countSE <- SummarizedExperiment::SummarizedExperiment(cellPopulation_Files,
    metadata = newMetadata
  )

  return(countSE)
}


## helper functions.
                         
# Generates rolling sum of insertions
                         
### Helper function should normalize Tn5 (or not), 
### then generate rolling sum of insretions
### Can generate overage after?
                    
normLocalInserts <- function(iterList) {
    
    insertList = iterList[[1]]
    normTn5 = iterList[[2]]
    insertBias = iterList[[3]]
    genome_db = iterList[[4]]
    sumWidth = iterList[[5]]
    medianWidth = iterList[[6]]
      
    ## Transform it into a GPos for a rolling sum
    gpos <- GenomicRanges::GPos(insertList)
    singleInsert <- rep(GenomicRanges::score(insertList),
             GenomicRanges::width(insertList))
    gpos$score <- singleInsert
    
    if(normTn5){
        
        
        genome <- getAnnotationDbFromInstalledPkgname(dbName = genome_db,
                                                  type = 'BSgenome')
        ### Normalize Insertion bias by
        ## 1. Identify pattern of
        ## 2. Pulling up the normalization score. 
        seqWindow = plyranges::reduce_ranges(insertList)
        seqWindow1 = seqWindow
        GenomicRanges::start(seqWindow1) = GenomicRanges::start(seqWindow1) - 3
        GenomicRanges::end(seqWindow1) = GenomicRanges::end(seqWindow1) + 2

        sequences <- Biostrings::getSeq(x = genome, seqWindow1)
        ## Now split the sequences into rolling 6-bp strings for each location, and pull out the normalized bias for each bp position.
        sequences2 <- unlist(lapply(as.vector(sequences), function(XX){
                        substring(XX, seq(1, GenomicRanges::width(seqWindow)[1], 1), 
                                  seq(1, GenomicRanges::width(seqWindow)[1], 1) + 5)
            }))

        singleInsert = singleInsert/insertBias[sequences2]
    }
    
    ######## Transform with smoothing filter ######## 
    ### Thank you to Imran McGrath for his work on this code.
    newx <- zoo::rollsum(singleInsert, k=sumWidth, align = 'center', fill = 0)
    newx <- suppressWarnings(zoo::rollmedian(newx, k=medianWidth, 
                                    align = 'center', fill = 0))
    gpos$score <- newx

    return(gpos)
}
                                     
     
# This function is not used.
## Efficiently subsets and bins coverage for a given region across samples
# subsetSlidingBinCoverage <- function(iterList) {
#   idx <- NULL
#   binnedData <- iterList[[1]]
#   regionGRanges <- plyranges::reduce_ranges(binnedData)
# 
#   tmpCounts <- lapply(subList, function(z) {
#     tmpGR <- plyranges::join_overlap_intersect(z, regionGRanges) %>%
#       tmpGR() <- plyranges::join_overlap_intersect(tmpGR, binnedData)
#     tmpGR <- plyranges::group_by(tmpGR, idx) %>%
#       tmpGR() <- plyranges::reduce_ranges(tmpGR, score = mean(score))
#     tmpGR <- dplyr::ungroup(tmpGR)
#     tmpGR
#   })
# 
#   return(tmpCounts)
# }
