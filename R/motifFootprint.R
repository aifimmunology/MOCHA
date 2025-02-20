#' @title Generate motif footprints
#'
#' @description Generate a plot of average normalized insertions around a motif center
#'
#' @param SampleTileObj A Sample-tile object from MOCHA's getSampleTileMatrix()
#' @param motifName The name of metadata entry with Motif location, added via addMotifSet()
#' @param specMotif An optional string specifying which motif to analyze. If blank, it will analyze all motifs. 
#' @param regions An optional GRanges object or list of strings in the format chr1:100-200, specifiying which specific regions to look at when conducting motif footprinting. 
#' @param cellPopulations A list of cell populations to conduct motif footprinting on. 
#' @param windowSize A number, representing the window to analyze around each motif location. Default is 500 bp. 
#' @param groupColumn A string, corresponding to a metadata column within the SampleTileObj, that describes the groups by which you want to summarize motif footprints. If sampleSpecific = FALSE, then motifFootprint will average insertions across samples within each group. 
#' @param subGroups A list of subgroups, if you want to only look at specific groups within the groupColumn. 
#' @param sampleSpecific A boolean for whether to generate average motif footprints within each group, or to return a data.frame for all samples. 
#' @param normTn5 A boolean for whether to normalize by Tn5 insertion bias. 
#' @param smoothTn5 The window size for smoothen Tn5 insertions at each location. Ideal when looking at motif footprints over a smaller number of regions, rarer cell types, or sparse regions, where local noise can make it harder to see the overall pattern. Can be set to 0.
#' @param numCores Number of cores to parallelize over
#' @param force Boolean. If FALSE, it will through an error if there's an empty sample or not overlap with regions and motifs. If TRUE will ignore these issues and continue (or return NULL)
#' @param verbose Boolean. Default is FALSE. Will print more messages if TRUE. 
#' 
#' @return A SummarizedExperiment containing motif footprinting data 
#'
#' @export
#' @keywords downstream


motifFootprint <- function(SampleTileObj,
                           motifName = 'Motifs',
                           specMotif = NULL,
                          regions = NULL,
                          cellPopulations = "ALL",
                          windowSize = 500,
                          normTn5 = TRUE,
                          smoothTn5 = 10,
                          groupColumn = NULL,
                          subGroups = NULL,
                          sampleSpecific = FALSE,
                          numCores = 1,
                           force = FALSE,
                          verbose = FALSE) {
  . <- idx <- score <- width <- Group <- Sample <- Location <- Index <- NULL
  cellNames <- SummarizedExperiment::assayNames(SampleTileObj)
  metaFile <- SummarizedExperiment::colData(SampleTileObj)
  outDir <- SampleTileObj@metadata$Directory

  if (is.na(outDir)) {
    stop("Missing coverage file directory. SampleTileObj$metadata must contain 'Directory'.")
  }
    
  #Make sure the window size is evenly divisible by 2. 
  if(windowSize/2 != round(windowSize/2)){
      windowSize = round(windowSize/2)*2
  }

  if (!file.exists(outDir)) {
    stop("Directory given by SampleTileObj@metadata$Directory does not exist.")
  }
    
  if (all(toupper(cellPopulations) == "ALL")) {
    cellPopulations <- cellNames
  }
  if (!all(cellPopulations %in% cellNames)) {
    stop("Some or all cell populations provided are not found.")
  }
    
  if(!is.null(smoothTn5)){
    if(smoothTn5 == 0){
        smoothTn5 = NULL
    }
  }
    
  ## Pull out cell type-specific motif positions, and subset by specific regions and/or to specific motifs.
  motifs <- getCellTypeMotifs(SampleTileObj, 
                              cellPopulation = cellPopulations,
                              specMotif = specMotif,
                              MotifSetName = motifName, asGRangesList = FALSE)
    
  if(!is.null(regions)){

      if (is.character(regions)) {
          
        # Convert to GRanges
        regionGRanges <- MOCHA::StringsToGRanges(regions)
      } else if (class(regions)[1] == "GRanges") {
        regionGRanges <- regions
      } else {
        stop("Wrong regions input type. regions must either be a string, or a GRanges location.")
      }
      #Combine nearby tiles
      regionGRanges2 = plyranges::reduce_ranges(regionGRanges)
      motifs2 <- lapply(motifs, function(XX){
                          plyranges::filter_by_overlaps(XX, regionGRanges2)
                })
      names(motifs2) = names(motifs)
      motifs = motifs2
      rm(motifs2)
      
      if(any(lengths(motifs) == 0) & !force){
          stop('No motif overlaps with open tiles for cell types and regions provided')
      }else if(all(lengths(motifs) == 0) & force){
      
          return(NULL)
          
      }else if(any(lengths(motifs) == 0) & force){      
          motifs = motifs[lengths(motifs) > 0]
      }
  }

  if (all(toupper(cellNames) == "COUNTS")) {
    stop(
      "The only assay in the SummarizedExperiment is Counts. The names of assays must reflect cell types,",
      " such as those in the Summarized Experiment output of getSampleTileMatrix."
    )
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
  } else {

    # If neither groupColumn nor subGroup is defined, then it forms one list of all sample names
    subGroups <- "All"
    subSamples <- list("All" = metaFile[, "Sample"])
  }
                    
  ### Verify that they've calculate insertion bias if they want to normalize by Insertion bias
  ### Calculate bias. 
  if(normTn5 & any(grepl('InsertionBias', names(SampleTileObj@metadata)))){
      ## Pull in genome database
      genome_db = SampleTileObj@metadata$Genome
      
      genome <- getAnnotationDbFromInstalledPkgname(dbName = genome_db, type = 'BSgenome')
      insertBias = SampleTileObj@metadata$InsertionBias
      #Remove any NAs that might be there
      insertBias = insertBias[!is.na(insertBias[,'Norm']),]
      
  }else if(normTn5 & !any(grepl('InsertionBias', names(SampleTileObj@metadata)))){
      
      stop('Attempting to normalize by Tn5 insertion bias, but no bias calculated. Please run addInsertionBias.')
      
  }
                      
  ## Create experimentList - this will be alist of matrices with insertions for each location/sample across positions.    
  experimentList1 = list()
  colData1_tmp = list()
                         
  #Pull up the cell types of interest, and filter for samples and subset down to region of interest  
  for(x in cellPopulations) {

    if (verbose) {
      message(stringr::str_interp("Extracting insertions from cell population '${x}'"))
    }
      
    # Pull up insertions files
    originalCovGRanges <- readRDS(paste(outDir, "/", x, "_CoverageFiles.RDS", sep = ""))
    if ('Insertions' %in% names(originalCovGRanges)){
      originalInsertions <- originalCovGRanges[['Insertions']]
      rm(originalCovGRanges)
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

    ### iterate over each motif, and subsample down to those regions. 
    for(YY in names(motifs)){
        
            ## Some motifs can be directly overlapping for hundreds of basepairs. 
            ## These motif matches are unlikely to yield useful results because a 
            ## transcription factor can't possible be bound at all of these sites at once, so then
            ## a motif footprint at these sites looses values, and becomes irrelevant (i.e. no footprint around motif)
            ## So we'll reduce out motifs to find these overlapping regions, and remove anything with more than
            ## windowSize/10 bps. 
        
            subMotifs = plyranges::filter(plyranges::reduce_ranges(motifs[[YY]]), width <= windowSize/2)
        
            if (verbose) {
              message(stringr::str_interp("Processing motif footprint for ${YY}."))
            }
        
            cl <- parallel::makeCluster(numCores)  
            stretchMotifs = plyranges::stretch(plyranges::anchor_center(
                                plyranges::reduce_ranges(subMotifs)), extend = windowSize*1.10)
            
            if(normTn5){
                
                newList =  lapply(names(originalInsertions), function(XX){ 
                    list(plyranges::filter_by_overlaps(originalInsertions[[XX]],  stretchMotifs), 
                         subMotifs, windowSize,                                               
                         insertBias[,'Norm', drop =FALSE], genome_db, smoothTn5, XX)})
                rm(stretchMotifs)
                rm(subMotifs)
                 allNorms <- pbapply::pblapply(cl = cl, X = newList, normMotifs2)
                
            }else{
                
               
                newList = lapply(names(originalInsertions), function(XX){ 
                    list(plyranges::filter_by_overlaps(originalInsertions[[XX]],  stretchMotifs),
                         subMotifs, windowSize, smoothTn5, XX)})
                rm(stretchMotifs)
                rm(subMotifs)
                allNorms <- pbapply::pblapply(cl = cl, X = newList, normMotifs)
                
                
            }
        
            parallel::stopCluster(cl)
        
            names(allNorms) = names(originalInsertions)
            #Clean up
            rm(newList)
            gc()
            allNorms = data.table::rbindlist(allNorms)
        
            if (!sampleSpecific) {
                names(subSamples) = subGroups
                allNorms[, Group := ifelse(Sample %in% subSamples[[1]],  names(subSamples)[1], NA)]
                for(subSample1 in names(subSamples)[2:length(subSamples)]){
                    allNorms[, Group := ifelse(Sample %in% subSamples[[subSample1]], subSample1, Group)]
                }
                
                ##Now generate the group-level location-specific average
                allNorms = allNorms[,.(score = mean(score, na.rm = TRUE)), by = .(Group, Position, Location)]
                setnames(allNorms, 'Group', 'Sample')            
            }
            colData2 = unique(allNorms[, c('Sample', 'Position'), with = FALSE])[
                            , Index := paste(Sample, Position, sep='__')]
            
            matrix1 = data.table::dcast(allNorms[, c('Sample', 'Position', 'score', 'Location'), with = FALSE],
                                              Location ~ Sample + Position, fun.aggregate = sum, 
                                                value.var = 'score', sep = "__")
            rm(allNorms)
            matrix2 = methods::as(matrix1[,!'Location'], 'matrix')
            rownames(matrix2) = matrix1[,Location]
            rm(matrix1)
            #CellType_Motif for index name
            list_index_name = paste(x, YY, sep ='__')
            # Add motif & cell type info
            #typeCast matrix2 as a sparse matrix
            experimentList1 = append(experimentList1, list(methods::as(matrix2, 'sparseMatrix')))
            names(experimentList1)[length(experimentList1)] = list_index_name
            
            if(!all(colData2$Index %in% colData1_tmp$Index)){
            
                subColData2 = colData2[!colData2$Index %in% colData1_tmp$Index,]
                colData1_tmp = cbind(colData1_tmp, subColData2)
                
            }
            
            rm(matrix2)

    }    
    rm(originalInsertions)
    gc()
    
    
  }
  newMetadata <- SampleTileObj@metadata
  newMetadata$History <- append(newMetadata$History, paste("motifFootprint", utils::packageVersion("MOCHA")))
                         
  ##Filter metadata down to essentials: summarizedData, metadata, Genome, TxDb, OrgDb, Directory, History, and TYpe
                         
  sampleMetadata = SummarizedExperiment::colData(SampleTileObj)

  newMetadata = newMetadata[c(1:6)]
                         
  newMetadata = append(newMetadata, list('Type' = 'Footprints', 'metadata' = sampleMetadata))
                         
  #Convert the data for the MultiAssayExperiment
  colData1_tmp= as.data.frame(colData1_tmp)
  rownames(colData1_tmp) = colData1_tmp$Index
                         
  footprintSE <- MultiAssayExperiment::MultiAssayExperiment(experiments = experimentList1,
                            colData = colData1_tmp,
                            metadata = newMetadata
                  )

  return(footprintSE)
}
                         
#' @title Add InsertionBias to MOCHA object
#'
#' @description Generated and Add InsertionBias to MOCHA object
#'
#' @param SampleTileObj A Sample-tile object from MOCHA's getSampleTileMatrix()
#' @param numCores Number of cores to parallelize oer
#' @param verbose A boolean to determin 
#'
#' @return A Sample-Tile object with an insertion bias matrix within the metadata list
#'
#' @export
#' @keywords downstream
        
                   
addInsertionBias <- function(SampleTileObj, numCores = 1, verbose = TRUE){


    if (!requireNamespace("Biostrings", quietly = TRUE)) {
          stop(
          "Package 'Biostrings' is required for calculating insertion bias.",
          "Please install Biostrings from Bioconductor to continue: BiocManager::install('Biostrings')"
          )
    }
    
    . <- Count <- score <- NULL
    ## Pull in Genome database
    genome_db = SampleTileObj@metadata$Genome
    outDir = SampleTileObj@metadata$Directory

    genome <- getAnnotationDbFromInstalledPkgname(dbName = genome_db, type = 'BSgenome')
    
    cellPopulations = SummarizedExperiment::assayNames(SampleTileObj)

    cellPopulation_Files <- lapply(cellPopulations, function(x) {
        
        if (verbose) {
          message(stringr::str_interp("Extracting insertions from cell population '${x}'"))
        }

        # Pull up insertions files
        originalCovGRanges <- readRDS(paste(outDir, "/", x, "_CoverageFiles.RDS", sep = ""))
        if ('Insertions' %in% names(originalCovGRanges)){
          originalInsertions <- originalCovGRanges[['Insertions']]
          rm(originalCovGRanges)
        }else{
         stop('Error around reading insertions files. Check that coverage/insertions files are not corrupted.')
        }
        
        insertList = lapply(originalInsertions, function(XX){
                    if(!is.null(XX)){
                        tmpGR = plyranges::filter(XX, score !=0)
                        if(length(tmpGR) > 10){
                            list(plyranges::filter(XX, score !=0), genome_db)
                        }else {NULL}
                    }else{NULL}
            })
         rm(originalInsertions)

         cl = parallel::makeCluster(numCores)
         sequencesTmp = pbapply::pblapply(cl = cl, X = insertList, getBias)
         parallel::stopCluster(cl)
        
         gc()
         rm(insertList)
         full_seq = do.call('rbind', sequencesTmp)
         full_seq = dplyr::summarize(dplyr::group_by(full_seq, seq),
                                     Count = sum(Count))
         full_seq$Norm = full_seq$Count/sum(full_seq$Count)
         rm(sequencesTmp)
         return(full_seq)
        
        })
    
    populationBias = dplyr::group_by(do.call('rbind', cellPopulation_Files), seq)
    populationBias = dplyr::summarize(populationBias, Count = sum(Count))
    populationBias$Norm = populationBias$Count/sum(populationBias$Count)
    biasMat = do.call('cbind', lapply(cellPopulation_Files, function(ZZ){
                        bindTmp = ZZ
                        if(any(!populationBias$seq %in% ZZ$seq)){
                             tmpMat = data.frame(seq = populationBias$seq[!populationBias$seq %in% ZZ$seq])
                             tmpMat$Count = 0
                             tmpMat$Norm = 0
                             bindTmp = rbind(bindTmp, tmpMat)
                        }
                        rownames(bindTmp) = bindTmp$seq
                        bindTmp[populationBias$seq, 'Norm']
        }))
    colnames(biasMat) = cellPopulations
    rownames(biasMat) = populationBias$seq
    biasMat = cbind(populationBias[,'Norm'], biasMat)
    
    genome_freq <- Biostrings::oligonucleotideFrequency(
        x = Biostrings::getSeq(x = genome, GenomicRanges::seqnames(genome)),
        width = 6
    )
    
    if (dim(genome_freq)[2] > 1) {
        genome_freq <- colSums(genome_freq)
    }
    genome_freq = genome_freq[!grepl('N',names(genome_freq))]
    genome_freq = genome_freq/sum(genome_freq)
                      
    if (any(!names(genome_freq) %in% populationBias$seq)) {
        warning("Some possible hexamers do not appear within observed insertions.")
        missingMers = names(genome_freq)[! names(genome_freq) %in% populationBias$seq]
        
        tmpMat = do.call('cbind', lapply(colnames(biasMat), function(XX){ rep(0, length(missingMers)) }))
        rownames(tmpMat) = missingMers
        biasMat = rbind(biasMat, tmpMat)
    }

   if (verbose) {
      message(stringr::str_interp("Normalizing hexamer bias across all celltypes'"))
    }
    biasMat2 = pbapply::pbapply(biasMat, 2, function(XX) 
        XX/genome_freq[rownames(biasMat)])
                      
    SampleTileObj@metadata = append(SampleTileObj@metadata,
                            list('InsertionBias' = biasMat2))
                                
    return(SampleTileObj)
}
    

     
#' @title Internal function for parallelizing motif footprinting over samples
#'
#' @description Generates normalized insertions for a given sample at all motif positions for a given motif, without Tn5 correction.
#' @param list1 A list, including A GRanges object of normalized insertions for one sample, motif positions, and window size.
#'
#' @return A tibble containing motif footprinting data 
#'
#' @noRd
      
normMotifs <- function(list1){

    position <- V1 <- partition <- score <- start <- end <- Mid <- seqnames <- NULL
    insertList = list1[[1]]
    motifPos = list1[[2]]
    windowSize = list1[[3]]
    smoothWindow = list1[[4]]
    sampleName = list1[[5]]
    
    #Add the middle position
    motifPos$Mid = round(c(GenomicRanges::end(motifPos) - GenomicRanges::start(motifPos))/2) +
                            GenomicRanges::start(motifPos)
    GenomicRanges::start(motifPos) = motifPos$Mid
    GenomicRanges::end(motifPos) = motifPos$Mid
    
    ##Generate the 500 bp window around the center of the motif
    windows1 = plyranges::anchor_center(motifPos)
    windows1 = plyranges::stretch(windows1, windowSize)
    rm(motifPos)
    
    ##Filter the insertions
    subInsert = plyranges::filter_by_overlaps(insertList, windows1)
    #Tile the windows into single basepairs
    windows2 = plyranges::tile_ranges(windows1, 1)
     # Using the single base-pair positions to get insertin counts, and motif mid points, 
    subInsert1 = plyranges::join_overlap_left(windows2, subInsert)
    if(!is.null(smoothWindow)){
        ## If we want to smoothen things: 
        ## Define smoothening function
         ## Define smoothening function
        smoothRegions <- function(insertGR, windowsGR, windowSize){
                data.table::setDTthreads(threads = 1)
                partition <- score <- NULL
                ### Let's run a for loop over each chromosome to slow this down and minimize memory usage
                windowsGR$partition = c(1:length(windowsGR))
                windowsGR$position = paste(GenomicRanges::seqnames(windowsGR), ":",GenomicRanges::start(windowsGR), sep ='')
                insertList = GenomicRanges::GRanges()
                for(chr in unique(GenomicRanges::seqnames(windowsGR))){
                    windowsGR2 = windowsGR[GenomicRanges::seqnames(windowsGR) == chr]
                    insertGR2 = insertGR[GenomicRanges::seqnames(insertGR) == chr]
                    windowsGR2 = plyranges::stretch(
                        plyranges::anchor_center(windowsGR2), windowSize)
                    ## Join expanded windows and inserts
                    insertDT = as.data.table(plyranges::join_overlap_left(insertGR2, windowsGR2))
                    #Use data.table to run the average
                    meanRow = insertDT[,sum(score, na.rm = TRUE),by=list(position, partition)]
                    rm(insertDT)
                    # Transfer the rolling sum score back to the GR object
                    positionList = paste(GenomicRanges::seqnames(insertGR2), ":",GenomicRanges::start(insertGR2), sep ='')
                    insertGR2$score = meanRow[,V1][match(positionList, meanRow[,position])]
                    rm(meanRow)
                    ## Repeat process for median over the windows
                    insertDT = as.data.table(plyranges::join_overlap_left(insertGR2, windowsGR2))
                    meanRow = insertDT[,stats::median(score, na.rm = TRUE),by=list(position, partition)]
                    rm(insertDT)
                    ## overwrite previous score with the new one.
                    insertGR2$score = meanRow[,V1][match(positionList, meanRow[,position])]
                    rm(meanRow)
                    insertList = c(insertList, insertGR2)
                    rm(insertGR2)
                    
                }
                return(insertList)
        }
        
        subInsert1 = plyranges::select(subInsert1, !partition)
        subInsert1 = smoothRegions(subInsert1, windowsGR = windows2, windowSize =smoothWindow)
    }
    rm(windows2)
    ##Remove score column (if present) from motif locations. 
    if(any(colnames(GenomicRanges::mcols(windows1)) == 'score')){
    
        windows1 <- plyranges::select(windows1, !score)
    }
    
    subInsert1 = plyranges::join_overlap_left(subInsert1,  windows1)
    rm(windows1)
    
    # Turn it into a data.frame, group it by mid point position, and 
    ## find the average over all positions
    subInsert1 = dplyr::mutate(as.data.frame(subInsert1), 
                               score = ifelse(is.na(score), 0, score), 
                  Position = start - Mid,
                  Location = paste(seqnames, ":", Mid, sep = ''), 
                  Sample = sampleName)
    subInsert1 = subInsert1[abs(subInsert1$Position) <= windowSize/2, c('Sample', 'Position', 'Location', 'score')]
    return(as.data.table(subInsert1))
}

                                
#' @title Second internal function for parallelizing motif footprinting over samples
#'
#' @description Generates normalized insertions for a given sample at all motif positions for a given motif, with Tn5 correction.
#' @param list1 A list, including A GRanges object of normalized insertions for one sample, motif positions, and window size.
#' 
#' @return A tibble containing motif footprinting data. 
#' 
#' @noRd
      
normMotifs2 <- function(list1){

    position <- V1 <- partition <- score <- start <- end <- Mid <- seqnames <- NULL
    
    insertList = list1[[1]]
    motifPos = list1[[2]]
    windowSize = list1[[3]]
    biasMat = list1[[4]]
    genome_db = list1[[5]]
    smoothWindow = list1[[6]]
    sampleName = list1[[7]]
    
    ##get database 
    genome <- getAnnotationDbFromInstalledPkgname(dbName = genome_db, type = 'BSgenome')

    #Add the middle position
    motifPos$Mid = round(c(GenomicRanges::end(motifPos) - GenomicRanges::start(motifPos))/2) +
                            GenomicRanges::start(motifPos)
    GenomicRanges::start(motifPos) = motifPos$Mid
    GenomicRanges::end(motifPos) = motifPos$Mid
    
    
    ##Generate the 500 bp window around the center of the motif
    windows1 = plyranges::anchor_center(motifPos)
    windows1 = plyranges::stretch(windows1, windowSize)
    rm(motifPos)
    ##Filter the insertions
    subInsert = plyranges::filter_by_overlaps(insertList, windows1)
    
    #Tile the windows into single basepairs
    windows2 = plyranges::tile_ranges(windows1, 1)
    
    seqLocations = windows1
    GenomicRanges::start(seqLocations) = GenomicRanges::start(seqLocations) - 3
    GenomicRanges::end(seqLocations) = GenomicRanges::end(seqLocations) + 2
    
    sequences2 <- Biostrings::getSeq(x = genome, seqLocations)
    ## Now split the sequences into rolling 6-bp strings for each location, and pull out the normalized bias for each bp position.
    sequences3 <- unlist(lapply(as.vector(sequences2), function(XX){
                    substring(XX, seq(1, GenomicRanges::width(windows1)[1], 1), 
                              seq(1, GenomicRanges::width(windows1)[1], 1) + 5)
        }))
    ## Unknown nucleotides are presented as an N, which we don't have information on. 
    ## When a give 6 bp window has an N, then set the norm value to NA. 
    normValues = rep(NA, times = length(windows2))
    
    normValues[sequences3 %in% rownames(biasMat)] = 
        biasMat[sequences3[sequences3 %in% rownames(biasMat)], 'Norm']
    
    ## Then add it to the windows2  GenomicRanges. 
    windows2$Norm = normValues   
    
    # Using the single base-pair positions to get insertion counts, and motif mid points, 
    subInsert1 = plyranges::join_overlap_left(windows2, subInsert)
    subInsert1$score = subInsert1$score/subInsert1$Norm
    
    if(!is.null(smoothWindow)){
        ## If we want to smoothen things: 
        
         ## Define smoothening function
        smoothRegions <- function(insertGR, windowsGR, windowSize){
                data.table::setDTthreads(threads = 1)
                partition <- score <- NULL
                ### Let's run a for loop over each chromosome to slow this down and minimize memory usage
                windowsGR$partition = c(1:length(windowsGR))
                windowsGR$position = paste(GenomicRanges::seqnames(windowsGR), ":",GenomicRanges::start(windowsGR), sep ='')
                insertList = GenomicRanges::GRanges()
                for(chr in unique(GenomicRanges::seqnames(windowsGR))){
                    windowsGR2 = windowsGR[GenomicRanges::seqnames(windowsGR) == chr]
                    insertGR2 = insertGR[GenomicRanges::seqnames(insertGR) == chr]
                    windowsGR2 = plyranges::stretch(
                        plyranges::anchor_center(windowsGR2), windowSize)
                    ## Join expanded windows and inserts
                    insertDT = as.data.table(plyranges::join_overlap_left(insertGR2, windowsGR2))
                    #Use data.table to run the average
                    meanRow = insertDT[,sum(score, na.rm = TRUE),by=list(position, partition)]
                    rm(insertDT)
                    # Transfer the rolling sum score back to the GR object
                    positionList = paste(GenomicRanges::seqnames(insertGR2), ":",GenomicRanges::start(insertGR2), sep ='')
                    insertGR2$score = meanRow[,V1][match(positionList, meanRow[,position])]
                    rm(meanRow)
                    ## Repeat process for median over the windows
                    insertDT = as.data.table(plyranges::join_overlap_left(insertGR2, windowsGR2))
                    meanRow = insertDT[,stats::median(score, na.rm = TRUE),by=list(position, partition)]
                    rm(insertDT)
                    ## overwrite previous score with the new one.
                    insertGR2$score = meanRow[,V1][match(positionList, meanRow[,position])]
                    rm(meanRow)
                    insertList = c(insertList, insertGR2)
                    rm(insertGR2)
                    
                }
                return(insertList)
        }
        
        subInsert1 = plyranges::select(subInsert1, !partition)
        subInsert1 = smoothRegions(subInsert1, windowsGR = windows2, windowSize =smoothWindow)

    }
    
    ##Remove score column (if present) from motif locations. 
    if(any(colnames(GenomicRanges::mcols(windows1)) == 'score')){
    
        windows1 <- plyranges::select(windows1, !score)
    }
    
    subInsert1 = plyranges::join_overlap_left(subInsert1, windows1)

    # Turn it into a data.frame, group it by mid point position, and 
    ## find the average over all positions
    subInsert1 = dplyr::mutate(as.data.frame(subInsert1), 
                  Position = start - Mid,
                  Location = paste(seqnames, ":", Mid, sep = ''), 
                  Sample = sampleName)
    subInsert1 = subInsert1[abs(subInsert1$Position) <= windowSize/2, c('Sample', 'Position', 'Location', 'score')]
    return(as.data.table(subInsert1))
}


                         
        
#' @title Internal function for parallelizing insertion bias within one sample.
#'
#' @description Function for getting hexamer occurrence frequency around insetions sites within one sample, meant for parallelization. Generates a data.frame of relative occurence of 6-bp sequences
#' @param insertList1 A list of A GRanges object of normalized insertions for one sample, and the genome database.
#'
#' @return A tibble containing relative occurence of sequences around insertion sites. 
#'
#' @noRd

getBias <- function(insertList1){
    value <- NULL
    insertions = insertList1[[1]]
    genome_db = insertList1[[2]]
    genome <- getAnnotationDbFromInstalledPkgname(dbName = genome_db, type = "BSgenome")
    ## We want to tile insertion coverage to single basepair locations, so that we can grab the sequences around that. 
    fInsert  = plyranges::tile_ranges(plyranges::anchor_center(insertions), width = 1)
    fInsert2 = plyranges::join_overlap_inner(fInsert, insertions)
    ##The minimum count should 1. By scaling by 1/minimum, you can undo the insertion normalization to get raw counts at each location.
    fInsert2$score = fInsert2$score*(1/min(fInsert2$score))
    insertions2 <- GenomicRanges::trim(plyranges::stretch(fInsert2, 5))
    
    filtGR = GenomicRanges::GRanges(seqnames = 
                        names(BSgenome::seqinfo(x = genome)), 
                     ranges = IRanges::IRanges(start = 1, 
                                end = BSgenome::seqinfo(genome)@seqlengths), 
                     strand = '*')
    
    insertions2 = plyranges::join_overlap_intersect(insertions2, filtGR)

    sequences <- as.vector(x = Biostrings::getSeq(x = genome, insertions2))
    #Weight the insertions by the coverage score present.
    seq_table = dplyr::summarize(dplyr::group_by(data.frame(seq = sequences, 
                                value = insertions2$score), seq),
                                 Count = sum(value))
    # remove sequences containing N or other wierd characters
    #seq_table = dplyr::filter(seq_table, !grepl(pattern = "N|[^ATCG]", seq))
    return(seq_table)
}
                                
#' @title \code{findFootprints}
#'
#' @description Test for footprints at motif positions within a given set of regions. 
#'              This function will basically take a window around a motif site and test if it's
#'              significantly different from a window around that site across all motif sites within
#'              a given set of sites across a given set of samples. 
#'
#' @param SampleTileObj A Sample-tile object from MOCHA's getSampleTileMatrix()
#' @param motifName The name of metadata entry with Motif location, added via addMotifSet()
#' @param specMotif An optional string specifying which motif to analyze. If blank, it will analyze all motifs within a given set of regions. 
#' @param regions An optional GRanges object or list of strings in the format chr1:100-200, specifiying which specific regions to look at when conducting motif footprinting. 
#' @param cellPopulations A list of cell populations to conduct motif footprinting on. 
#' @param subGroups A list of subgroups, if you want to only look at specific groups within the groupColumn. 
#' @param sampleSpecific A boolean for whether to generate average motif footprints within each group, or to return a data.frame for all samples. 
#' @param normTn5 A boolean for whether to normalize by Tn5 insertion bias. 
#' @param smoothTn5 The window size for smoothen Tn5 insertions at each location. Ideal when looking at motif footprints over a smaller number of regions, rarer cell types, or sparse regions, where local noise can make it harder to see the overall pattern. Can be set to 0.
#' @param footprintSize An estimate of the footprint size that you want to test for. 
#' @param flankPoint The estimated distance between the motif center and the footprint flank (right and left). A window of equal size to the footprint size will be used to estimate background insertions at a given site. 
#' @param numCores Number of cores to parallelize over
#' @param force Boolean. If FALSE, it will through an error if there's an empty sample or not overlap with regions and motifs. If TRUE will ignore these issues and continue (or return NULL)
#'
#' @return A SummarizedExperiment containing motif footprinting data 
#'
#' @noRd
#' @keywords downstream
findFootprints <- function(SampleTileObj,
                           motifName = 'Motifs',
                          regions = NULL,
                          cellPopulations = "ALL",
                          specMotif = NULL,
                          footprintSize = 20,
                          windowSize = 100,
                          normTn5 = TRUE,
                          smoothTn5 = 10,
                          groupColumn = NULL,
                          subGroups = NULL,
                          sampleSpecific = FALSE,
                          numCores = 1,
                           force = FALSE,
                          verbose = FALSE) {
  . <- idx <- score <- Group <- Sample <- Location <- pval <- pval_adj <- NULL
  cellNames <- SummarizedExperiment::assayNames(SampleTileObj)
  metaFile <- SummarizedExperiment::colData(SampleTileObj)
  outDir <- SampleTileObj@metadata$Directory

  if (is.na(outDir)) {
    stop("Missing coverage file directory. SampleTileObj$metadata must contain 'Directory'.")
  }
    
  if (!file.exists(outDir)) {
    stop("Directory given by SampleTileObj@metadata$Directory does not exist.")
  }
    
  if (all(toupper(cellPopulations) == "ALL")) {
    cellPopulations <- cellNames
  }
  if (!all(cellPopulations %in% cellNames)) {
    stop("Some or all cell populations provided are not found.")
  }
    
  if(!is.null(smoothTn5)){
    if(smoothTn5 == 0){
        smoothTn5 = NULL
    }
  }
    
  if (all(toupper(cellNames) == "COUNTS")) {
    stop(
      "The only assay in the SummarizedExperiment is Counts. The names of assays must reflect cell types,",
      " such as those in the Summarized Experiment output of getSampleTileMatrix."
    )
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
  } else {

    # If neither groupColumn nor subGroup is defined, then it forms one list of all sample names
    subGroups <- "All"
    subSamples <- list("All" = metaFile[, "Sample"])
  }
                    
  ### Verify that they've calculate insertion bias if they want to normalize by Insertion bias
  ### Calculate bias. 
  if(normTn5 & any(grepl('InsertionBias', names(SampleTileObj@metadata)))){
      ## Pull in genome database
      genome_db = SampleTileObj@metadata$Genome
      
      genome <- getAnnotationDbFromInstalledPkgname(dbName = genome_db, type = 'BSgenome')
      insertBias = SampleTileObj@metadata$InsertionBias
      #Remove any NAs that might be there
      insertBias = insertBias[!is.na(insertBias[,'Norm']),]
      
  }else if(normTn5 & !any(grepl('InsertionBias', names(SampleTileObj@metadata)))){
      
      stop('Attempting to normalize by Tn5 insertion bias, but no bias calculated. Please run addInsertionBias.')
      
  }
                      
  ## Create experimentList - this will be alist of matrices with insertions for each location/sample across positions.    
  experimentList1 = list()

  #Pull out sample metadata
  sampleMetadata = SummarizedExperiment::colData(SampleTileObj)
                         
  #Pull up the cell types of interest, and filter for samples and subset down to region of interest  
  for(x in cellPopulations) {


     ## Pull out cell type-specific motif positions, and subset by specific regions and/or to specific motifs.
     motifs <- getCellTypeMotifs(SampleTileObj, 
                                  cellPopulation = x,
                                  specMotif = specMotif,
                                  MotifSetName = motifName, asGRangesList = FALSE)
     ## Now filter motifs by regions provided.   
     if(!is.null(regions)){

      if (is.character(regions)) {
          
        # Convert to GRanges
        regionGRanges <- MOCHA::StringsToGRanges(regions)
          
      } else if (class(regions)[1] == "GRanges") {
          
        regionGRanges <- regions
          
      } else {
          
        stop("Wrong regions input type. regions must either be a string, or a GRanges location.")
          
      }
      #Combine nearby tiles
      regionGRanges2 = plyranges::reduce_ranges(regionGRanges)
      motifs2 <- lapply(motifs, function(XX){
                          plyranges::filter_by_overlaps(XX, regionGRanges2)
                })
      names(motifs2) = names(motifs)
      motifs = motifs2
      rm(motifs2)
      
      if(any(lengths(motifs) == 0) & !force){
          stop('No motif overlaps with open tiles for cell types and regions provided')
      }else if(all(lengths(motifs) == 0) & force){
      
          return(NULL)
          
      }else if(any(lengths(motifs) == 0) & force){      
          motifs = motifs[lengths(motifs) > 0]
      }
    }
      
    if (verbose) {
      message(stringr::str_interp("Extracting insertions from cell population '${x}'"))
    }
      
    # Pull up insertions files
    originalCovGRanges <- readRDS(paste(outDir, "/", x, "_CoverageFiles.RDS", sep = ""))
    if ('Insertions' %in% names(originalCovGRanges)){
      originalInsertions <- originalCovGRanges[['Insertions']]
      rm(originalCovGRanges)
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

    ### iterate over each motif, and subsample down to those regions. 
    for(YY in names(motifs)){
        
            ## Some motifs can be directly overlapping for hundreds of basepairs. 
            ## These motif matches are unlikely to yield useful results because a 
            ## transcription factor can't possible be bound at all of these sites at once, so then
            ## a motif footprint at these sites looses values, and becomes irrelevant (i.e. no footprint around motif)
            ## So we'll reduce out motifs to find these overlapping regions, and remove anything with more than
            ## windowSize/10 bps. 
        
            subMotifs = plyranges::reduce_ranges(motifs[[YY]])
            
            if (verbose) {
              message(stringr::str_interp("Processing motif footprint for ${YY}."))
            }
        
            cl <- parallel::makeCluster(numCores)  
            stretchMotifs = plyranges::stretch(plyranges::anchor_center(
                                plyranges::reduce_ranges(subMotifs)), extend = windowSize*1.10)
            
            if(normTn5){
                
                newList =  lapply(names(originalInsertions), function(XX){ 
                    list(plyranges::filter_by_overlaps(originalInsertions[[XX]],  stretchMotifs), 
                         subMotifs, windowSize,                                               
                         insertBias[,'Norm', drop =FALSE], genome_db, smoothTn5, XX)})
                rm(stretchMotifs)
                rm(subMotifs)
                 allNorms <- pbapply::pblapply(cl = cl, X = newList, normMotifs2)
                
            }else{
                
               
                newList = lapply(names(originalInsertions), function(XX){ 
                    list(plyranges::filter_by_overlaps(originalInsertions[[XX]],  stretchMotifs),
                         subMotifs, windowSize, smoothTn5, XX)})
                rm(stretchMotifs)
                rm(subMotifs)
                allNorms <- pbapply::pblapply(cl = cl, X = newList, normMotifs)
                
                
            }
        
            parallel::stopCluster(cl)
        
            names(allNorms) = names(originalInsertions)
            #Clean up
            rm(newList)
            gc()
            allNorms = data.table::rbindlist(allNorms)
        
            if (!sampleSpecific) {
                names(subSamples) = subGroups
                allNorms[, Group := ifelse(Sample %in% subSamples[[1]],  names(subSamples)[1], NA)]
                for(subSample1 in names(subSamples)[2:length(subSamples)]){
                    allNorms[, Group := ifelse(Sample %in% subSamples[[subSample1]], subSample1, Group)]
                }
                
                ##Now generate the group-level location-specific average
                allNorms = allNorms[,.(score = mean(score, na.rm = TRUE)), by = .(Group, Position, Location)]
                setnames(allNorms, 'Group', 'Sample')            
            }
           
            ## Get the wilcoxon test on all centers by sample. 
            locCore = allNorms[abs(Position) <= footprintSize/2,.(score = mean(score), Region = 'Motif'), 
                               by = c('Sample', 'Location')]
            locFlank = allNorms[abs(Position) >= footprintSize/2,.(score = mean(score), Region = 'Flank'), 
                                by = c('Sample', 'Location')]
            fullMerge = data.table::dcast(rbind(locCore, locFlank)[, 
                                    pval := stats::wilcox.test(score ~ Region, .SD)$p.value, Sample],
                                   Sample + pval ~ Region, fun = list(mean,stats::sd), value.var = "score")
            fullMerge =  as.data.frame(fullMerge[,pval_adj := stats::p.adjust(pval, method = 'BH')])
            
            #CellType_Motif for index name
            list_index_name = paste(x, YY, sep ='__')
            # Add motif & cell type info
            rownames(fullMerge) = fullMerge$Sample
            fullMerge$CellType = x
            fullMerge$Motif = YY
            fullMerge = dplyr::inner_join(fullMerge[rownames(sampleMetadata),], as.data.frame(sampleMetadata), by = 'Sample')
            
            experimentList1 = append(experimentList1, list(fullMerge))
            names(experimentList1)[length(experimentList1)] = list_index_name
            rm(fullMerge)

    }    
    rm(originalInsertions)
    gc()
    
    
  }
  if(length(experimentList1) > 1){
      fullMat = do.call('rbind', experimentList1)
    }else{
      fullMat = experimentList1[[1]]
    }

  return(fullMat)
}
         
