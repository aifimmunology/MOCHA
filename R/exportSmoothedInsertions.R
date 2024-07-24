#' @title Export insertion counts to per-sample BigWig files after applying
#'   a rolling sum and rolling median smoothing filter.
#'
#' @description \code{exportLocalFootprints} Takes a SampleTileMatrix with
#'   linked insertion files, corrects Tn5 insertion bias, and applies a smoothing filter (a rolling sum then
#'   rolling median) to the insertions. This generates a local footprint track for visually identifying evidence of potential binding sites. 
#'   These files are then written to bigwig format.
#'
#' @param SampleTileObj A MultiAssayExperiment or RangedSummarizedExperiment
#'   from MOCHA
#' @param cellPopulation A string denoting the cell population of interest
#' @param outDir Directory to write output bigwig files. Default is NULL, where
#'   the directory in `SampleTileObj@metadata$Directory` will be used.
#' @param windowSize Window size for rolling sum & median over basepairs. Default is 10.
#' @param normTn5 A boolean for whether to normalize by Tn5 insertion bias. 
#' @param force Set TRUE to overwrite existing files. Default is FALSE.
#' @param slow Set TRUE to bypass optimisations and compute smoothing filter
#'   directly on the whole genome. May run slower and consume more RAM. Default
#'   is FALSE.
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#'
#' @return outPaths List of paths of exported insertion files
#'
#' @examples
#' \dontrun{
#' # Depends on and manipulates files on filesystem
#' outPath <- MOCHA::exportSmoothedInsertions(
#'   SampleTileObj,
#'   cellPopulation = "CD4 Naive", sumWidth = 10, medianWidth = 11, verbose = FALSE
#' )
#' }
#'
#' @export
#' @keywords exporting
exportLocalFootprints <- function(SampleTileObj,
                                     cellPopulation,
                                     outDir = NULL,
                                     windowSize = 10,
                                     groupColumn = NULL,
                                     subGroups = NULL,
                                     sampleSpecific = FALSE,
                                     normTn5 = TRUE,
                                     force=FALSE,
                                     slow=FALSE,
                                     verbose=FALSE,
                                     numCores = 5
                                    ){
    
    score <- start <- seqnames <- NULL
    if (!any(names(SummarizedExperiment::assays(SampleTileObj)) %in% cellPopulation)) {
        stop("cellPopulation was not found within SampleTileObj. Check available cell populations with `colData(SampleTileObj)`.")
    }

    sourcedir <- SampleTileObj@metadata$Directory
    metaFile <- SummarizedExperiment::colData(SampleTileObj)
    
    
    if(numCores > 6){
    
        warning("You have set the number of cores (numCores) to be greater than 5. Be careful because more parallelization can lead to much higher RAM usage. We recommend less than 5 cores. ")
        
    }
    
    ## Check and set directory for output. 
    if (is.null(outDir)){
        outdir <- SampleTileObj@metadata$Directory
    } else if (!dir.exists(outDir)){
        dir.create(outDir)
    }
    
    ### Pull out Tn5 insert bias, if necessary.
    if(normTn5 & any(grepl('InsertionBias', names(SampleTileObj@metadata)))){
          ## Pull in genome database
          genome_db = SampleTileObj@metadata$Genome
          genome = getAnnotationDbFromInstalledPkgname(dbName = genome_db, type = 'BSgenome')
          insertBias = SampleTileObj@metadata$InsertionBias

      }else if(normTn5 & !any(grepl('InsertionBias', names(SampleTileObj@metadata)))){

          stop('Attempting to normalize by Tn5 insertion bias, but no bias calculated. Please run addInsertionBias.')

      }else{
        
         insertBias = NULL
        
      }
    
    covFile <- file.path(sourcedir, paste0(cellPopulation, "_CoverageFiles.RDS"))
    if (!file.exists(covFile)){
        stop("Coverage file ", covFile, " for could not be found. ",
            "Please ensure that the directory given by `SampleTileObj@metadata$Directory` exists and contains coverage files.", 
            " To preserved linked files when sharing MOCHA objects, ", 
            " use `packMOCHA()` and `unpackMOCHA()`.")
    }
    insertionsList <- readRDS(covFile)
    if (!"Insertions" %in% names(insertionsList)) {
        stop("No insertions found in the file `", covFile, "`.",
             " Please run peak calling with the latest MOCHA version,", 
             " 1.0.1 or greater.")
    }
    
    insertionsGRangesList <- insertionsList$Insertions
    rm(insertionsList)
    
    # Pull out a list of samples by group.
    if (!is.null(subGroups) & !is.null(groupColumn)) {
        # If the user defined a list of subgroup(s) within the groupColumn from the metadata, 
        # then it subsets to just those samples
        
        subSamples <- lapply(subGroups, function(x) metaFile[metaFile[, groupColumn] %in% x, "Sample"])
        names(subSamples) <- subGroups
                             
    } else if (!is.null(groupColumn)) {

        # If no subGroup defined, then it'll form a list of samples across all labels within the groupColumn
        subGroups <- unique(metaFile[, groupColumn])
        subSamples <- lapply(subGroups, function(x) metaFile[metaFile[, groupColumn] %in% x, "Sample"])
        names(subSamples) = subGroups
                             
    } else {

        # If neither groupColumn nor subGroup is defined, then it forms one list of all sample names
        subGroups <- "All"
        subSamples <- list("All" = metaFile[, "Sample"])
        
    }
    
    ## Filter to specific samples of interest from the SampleTileObject
    for(y in seq_along(subSamples)) {
      if (!all(subSamples[[y]] %in% names(insertionsGRangesList))) {
        missingSamples <- paste(subSamples[[y]][!subSamples[[y]] %in% 
                                        names(insertionsGRangesList)], collapse = ", ")
          
        if(!force){
            stop(stringr::str_interp(c(
              "There is no insertion coverage for cell population '${x}' in the ",
              "following samples in sample grouping '${names(subSamples)[y]}': ",
              "${missingSamples}"
            )))
        }else{
            
            subSamples[[y]] = subSamples[[y]][subSamples[[y]] %in% 
                                              names(insertionsGRangesList)]
        }
      }
    }
      
    ##Subset insertions down to only the relevant samples.
    insertionsGRangesList =  insertionsGRangesList[unlist(subSamples)]
    
    outPaths <- list()
    
    if(!sampleSpecific){
    
        sampleList = list()
        
    }
    
    for (sampleName in names(insertionsGRangesList)) {
        
        ##Print out the name of the insertion file. 
        type <- paste0("window", windowSize)
        outfile <- file.path(outdir, paste0(
            cellPopulation, "__", sampleName, "__", type, ".bw"
        ))
        
        if (file.exists(outfile) && !force) {
            if (verbose) { message("Output file already exists: ", outfile) }
            next # Skip this one!
        }
        if (verbose) { message("Smoothing insertions for ", outfile, " ... ") }

        # Create iteration list        
        iterList <- lapply(levels(GenomicRanges::seqnames(insertionsGRangesList[[sampleName]])), function(XX){
                tmpGR = insertionsGRangesList[[sampleName]][GenomicRanges::seqnames(insertionsGRangesList[[sampleName]]) == XX]
                list(tmpGR, insertBias, windowSize, genome_db)
            })
        
        cl <- parallel::makeCluster(numCores)
        
        #Now process the data. 
        chrList <- pbapply::pblapply(cl = cl, X = iterList, processFile)
        
        allChrGR <-do.call('c', chrList)
        
        parallel::stopCluster(cl)
        
        rm(chrList)
        rm(iterList)
        gc()

        if(sampleSpecific){
            result <- tryCatch({
                plyranges::write_bigwig(allChrGR, outfile)
            }, warning = function(w) {
                message(w)
            }, error = function(e) {
                message("Error occurred writing the bigwig for sample `", sampleName, "` :")
                message(e)
                file.remove(outfile)
                return(allChrGR)
            }, finally = {
                outPaths <- append(outPaths, outfile)
            })
        }else{
        
            sampleList = append(sampleList, list(allChrGR))
            names(sampleList)[length(sampleList)] = sampleName
            
        }
    }
                

    ## Generate group averages
    outPaths = list()
    
  if(!sampleSpecific){
      
      message("Generating group-averages for local footprints and writing to disk.")
      
      cellPopSubsampleCov <- pbapply::pblapply(cl = NULL, X = subSamples, function(y){
          averageCoverage(sampleList[y])
          })
      
      names(cellPopSubsampleCov) <- subGroups
      
      for(one_group in subGroups){
          
        outfile <- file.path(outdir, paste0(
            cellPopulation, "__", groupColumn, "_", one_group, "__", type, ".bw"
        ))
        outfile <- gsub(" |\\.","_", outfile)
        
        if (file.exists(outfile) && !force) {
            if (verbose) { message("Output file already exists: ", outfile) }
            next # Skip this one!
        }
          
          
        result <- tryCatch({
                    plyranges::write_bigwig(cellPopSubsampleCov[[one_group]], outfile)
                }, warning = function(w) {
                    message(w)
                }, error = function(e) {
                    message("Error occurred writing the bigwig for sample `", sampleName, "` :")
                    message(e)
                    file.remove(outfile)
                    return(allChrGR)
                }, finally = {
                    outPaths <- append(outPaths, outfile)
                })
          
      }

  }
    
  return(outPaths)
}

           
#' @title smoothChromosome, an internal function for smoothening insertions for a full chromosome
#'
#' @description Takes a GRanges of single basepair insertions and smoothens it to deal with sparsity. 
#'
#' @param objectList A list of GenomicRanges object with insertions, insertionBias matrix, and windowsize for smoothening. 
#' 
#' @return A GenomicRanges object with corrected and smoothened insertion counts. 
#'
#' @noRd
           
processFile <- function(objectList){
    
    chrIns <- objectList[[1]]
    insertBias <- objectList[[2]]
    windowSize <- objectList[[3]]
    genome_db <- objectList[[4]]

    ##get database 
    genome <- getAnnotationDbFromInstalledPkgname(dbName = genome_db, type = 'BSgenome')
    
    onlyIns <- plyranges::filter(chrIns, score !=0)
    ## Find windows to pay attention to:
    windowsGR <-plyranges::reduce_ranges(plyranges::stretch(onlyIns, extend = 2*windowSize))
    ##Tile those windows. 
    tilesGR <- plyranges::tile_ranges(windowsGR, width = 1)
    #Find insertions at each bp within those windows
    bpInsert = plyranges::join_overlap_left(tilesGR, chrIns)
    bpInsert$score[is.na(bpInsert$score)] = 0
    
    ## Now normalize by tn5 if necessary
    if(!is.null(insertBias)){
        seqLocations = windowsGR
        GenomicRanges::start(seqLocations) = GenomicRanges::start(seqLocations) - 3
        GenomicRanges::end(seqLocations) = GenomicRanges::end(seqLocations) + 2

        sequences2 <- Biostrings::getSeq(x = genome, seqLocations)
        ## Now split the sequences into rolling 6-bp strings for each location, and pull out the normalized bias for each bp position.
        sequences3 <- unlist(lapply(as.vector(sequences2), function(XX){
                        lenSeq = nchar(XX)-5
                        substring(XX, seq(1, lenSeq, 1), 
                                  seq(1, lenSeq, 1) + 5)
            }))
        rm(sequences2)
        ## Unknown nucleotides are presented as an N, which we don't have information on. 
        ## When a give 6 bp window has an N, then set the norm value to NA. 
        normValues = rep(NA, times = length(tilesGR))
        normValues[sequences3 %in% rownames(insertBias)] = 
            insertBias[sequences3[sequences3 %in% rownames(insertBias)], 'Norm']
        ## Then normalize the insertion counts
        bpInsert$score = bpInsert$score/normValues
        bpInsert$score[is.na(bpInsert$score)] = 0
        rm(normValues)
        rm(sequences3)
        
    }
    rm(windowsGR)

    ## Now conduct smoothening
    tilesGR$partition = c(1:length(tilesGR))
    tilesGR$position = paste(GenomicRanges::seqnames(tilesGR), ":",GenomicRanges::start(tilesGR), sep ='')
    tilesGR2 = plyranges::stretch(
            plyranges::anchor_center(tilesGR), windowSize)
    rm(tilesGR)
    #Convert to Data.table for processing speed. 
    bpInsert = plyranges::select(bpInsert, !partition)
    insertDT = as.data.table(plyranges::join_overlap_left(bpInsert, tilesGR2))

    #Use data.table to run the average
    partitionDT = insertDT[,sum(score, na.rm = TRUE),by=list(position, partition)]
    # Transfer the rolling sum score back to the GR object
    positionList = paste(GenomicRanges::seqnames(bpInsert), ":",GenomicRanges::start(bpInsert), sep ='')
    bpInsert$score = partitionDT[,V1][match(positionList, partitionDT[,position])]

    ## Repeat process for median over the windows
    insertDT = as.data.table(plyranges::join_overlap_left(bpInsert, tilesGR2))
    partitionDT = insertDT[,median(score, na.rm = TRUE),by=list(position, partition)]
     # Transfer the rolling median score back to the GR object
    bpInsert$score = partitionDT[,V1][match(positionList, partitionDT[,position])]
    bpInsert <- plyranges::compute_coverage(bpInsert, weight= bpInsert$score)
    GenomeInfoDb::seqinfo(bpInsert) = GenomeInfoDb::seqinfo(chrIns)
    rm(partitionDT)
    rm(tilesGR2)
    return(bpInsert)
}
