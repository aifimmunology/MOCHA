#' @title Run Zero-inflated Generalized Linear Mixed Modeling on pseudobulked
#'   scATAC data
#'
#' @description \code{runZIGLMM} Runs linear mixed-effects modeling for
#'   zero-inflated data using \code{\link[glmmTMB]{glmmTMB}}.
#'
#' @param TSAM_Object A SummarizedExperiment object generated from
#'   getSampleTileMatrix.
#' @param cellPopulation Name of a cell type(s), or 'all'. The function will
#'   combine the cell types mentioned into one matrix before running the model.
#' @param continuousFormula The formula for the continuous data that should be
#'   used within glmmTMB. It should be in the format (exp ~ factors). All
#'   factors must be found in column names of the TSAM_Object metadata, except
#'   for CellType, FragNumber and CellCount, which will be extracted from the
#'   TSAM_Object. modelFormula must start with 'exp' as the response. See
#'   \link[glmmTMB]{glmmTMB}.
#' @param ziformula The formula for the zero-inflated data that should be used
#'   within glmmTMB. It should be in the format ( ~ factors). All factors must
#'   be found in column names of the TSAM_Object colData metadata, except for
#'   CellType, FragNumber and CellCount, which will be extracted from the
#'   TSAM_Object.
#' @param zi_threshold Zero-inflated threshold ( range = 0-1), representing the
#'   fraction of samples with zeros. When the percentage of zeros in the tile is
#'   between 0 and zi_threshold, samples with zeroes are dropped and only the
#'   continous formula is used. Use this parameter at your own risk. Default is
#'   0.
#' @param initialSampling Size of data to use for pilot
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @param numCores integer. Number of cores to parallelize across.
#'
#' @return results a SummarizedExperiment containing LMEM results
#'

simulateMetaData <- function(simulatedFragmentData, sampleNumber = 1){

    subMeta <- as.data.frame(GenomicRanges::mcols(simulatedFragmentData))
    subMeta <- dplyr::distinct(dplyr::select(subMeta, !c(inPeak,PeakID)))
    allCells <- unique(subMeta$CellID)

    if(sampleNumber > 1){
        averageSampleRatio = length(allCells)/sampleNumber
        cellNumbers = rpois(n = sampleNumber, lambda = ceiling(averageSampleRatio))
        normCellNumbers = round(cellNumbers/sum(cellNumbers)*length(allCells))
        while(sum(normCellNumbers) != length(allCells)){
            cellNumbers = rpois(n = sampleNumber, lambda = ceiling(averageSampleRatio))
            normCellNumbers = round(cellNumbers/sum(cellNumbers)*length(allCells))
        }

        sampleNames = unlist(lapply(1:sampleNumber, function(x){ paste(sample(LETTERS, 18, replace= TRUE), collapse = '')}))
        cellMatch <- unlist(lapply(seq_along(sampleNames), function(XX) rep(sampleNames[XX],normCellNumbers[XX])))
        subMeta$Sample = cellMatch
        
    }else{

        subMeta$Sample <- paste(sample(LETTERS, 20, replace= TRUE), collapse = '')
        
    }

    rownames(subMeta) <- subMeta$CellID
    
    return(subMeta)

}


#' @title simulate scATAC-seq fragments
#'
#' @description \code{runZIGLMM} Runs linear mixed-effects modeling for
#'   zero-inflated data using \code{\link[glmmTMB]{glmmTMB}}.
#'
#' @param TSAM_Object A SummarizedExperiment object generated from
#'   getSampleTileMatrix.
#' @param cellPopulation Name of a cell type(s), or 'all'. The function will
#'   combine the cell types mentioned into one matrix before running the model.
#' @param continuousFormula The formula for the continuous data that should be
#'   used within glmmTMB. It should be in the format (exp ~ factors). All
#'   factors must be found in column names of the TSAM_Object metadata, except
#'   for CellType, FragNumber and CellCount, which will be extracted from the
#'   TSAM_Object. modelFormula must start with 'exp' as the response. See
#'   \link[glmmTMB]{glmmTMB}.
#' @param ziformula The formula for the zero-inflated data that should be used
#'   within glmmTMB. It should be in the format ( ~ factors). All factors must
#'   be found in column names of the TSAM_Object colData metadata, except for
#'   CellType, FragNumber and CellCount, which will be extracted from the
#'   TSAM_Object.
#' @param zi_threshold Zero-inflated threshold ( range = 0-1), representing the
#'   fraction of samples with zeros. When the percentage of zeros in the tile is
#'   between 0 and zi_threshold, samples with zeroes are dropped and only the
#'   continous formula is used. Use this parameter at your own risk. Default is
#'   0.
#' @param initialSampling Size of data to use for pilot
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @param numCores integer. Number of cores to parallelize across.
#'
#' @return results a SummarizedExperiment containing LMEM results
#'

simulateFragments <- function(nCells = 500, meanFragsPerCell = 5000, fragThreshold = 1000, 
                                peakNumber = 500, peakCenters = 1,
                                largePeakWindow = 500, Genome = BSgenome.Hsapiens.UCSC.hg38,
                                FRIP = 0.9, meanLengths = c(75, 200), lengthProbability = c(0.9, 0.1)){


    if(length(meanLengths) != length(lengthProbability)){
        stop('meanLengths and lengthProbability must be of equal length.')
    }
    if(FRIP <= 0 | FRIP > 1){
        stop('FRIP score must be greater than zero and less then or equal to 1.')
    }


    #########################################################################
    ## Generate distribution of fragments per cell. 
    #########################################################################

    fragNumber = rpois(n = nCells, lambda = meanFragsPerCell)
    while(any(fragNumber < fragThreshold)){
        fragNumber[fragNumber < fragThreshold] = rpois(n = sum(fragNumber < fragThreshold), lambda = meanFragsPerCell)
    }

    #########################################################################
    ## Generate fragment number and widths for each simulated cell.
    #########################################################################
    
    ## Iterate over each cell, and generate the corresponding number of fragments with 
    ## the fragment length distribution matching the expectations. 

    ## Generate unique CellIDs
    allCellNames <- unlist(lapply(1:nCells, function(x){ paste(sample(LETTERS, 20 , replace= TRUE), collapse = '')}))
    peakWeight <- runif(n = peakNumber, min = 1, max = 8)

    message('Generating fragment lengths across theoretical cells.')

    fragList <- pbapply::pblapply(cl = NULL, X = seq_along(fragNumber), function(XX){

        # If things don't divide exactly, there will be one less fragment than technically generated,
        # which is within tolerances. 
        allLengths <- unlist(lapply(seq_along(meanLengths), function(YY){

            rpois(n = round(fragNumber[[XX]]*lengthProbability[[YY]]), lambda = meanLengths[[YY]])

        }))
        fragNum_tmp = length(allLengths)
        ## Classify the first portion as in peak, according to FRIP score, and the rest outside of the peak.
        InPeak = c(1:fragNum_tmp) <= rpois(1,FRIP*fragNum_tmp)
        ## Classify all the in-peak fragments to a given largePeakWindow
        PeakWindow <- sample(1:peakNumber, size = sum(InPeak), prob = peakWeight, replace = TRUE)

        peakID <- NA
        peakID[InPeak] <- PeakWindow

        data.frame(width = sample(allLengths, size = length(allLengths)), 
                 CellID = rep(allCellNames[XX], length(allLengths)),
                inPeak = InPeak, PeakID = peakID)

    })

    fragDF <- do.call('rbind', fragList)
    #Now let's see how many fragments fell within peaks across samples
    subFragNumber <- unlist(lapply(fragList, function(x) sum(x$inPeak)))
    #What are the peak intensities?
    peakIntensities_tmp <- as.data.frame(table(fragDF$PeakID))
    peakIntensities <- peakIntensities_tmp$Freq
    names(peakIntensities) <- peakIntensities_tmp$Var1

    #########################################################################
    ## Generate peak locations and background locations across the genome
    #########################################################################

    ### Now let's determine where our peaks are.
    ChromLengths = GenomeInfoDb::seqlengths(Genome)
    ChromLengths <- ChromLengths[!grepl('_', names(ChromLengths))]

    ### Generate all possible bins, and mark peak bins at evenly spaced intervals. 
    allBins = sum(ChromLengths)/largePeakWindow

    if(allBins <= peakNumber){
        stop('peakNumber is greater than or equal to the total number of larger open regions possible. Please change peakNumber, largePeakWindow, or both.')
    }

    binDistance = allBins/peakNumber
    peakLocation = seq(1, allBins, by = floor(binDistance))[1:peakNumber]
    startBin = 500

    message('Placing large peak windows across genome.')
    #This arrangement may slightly depass the ends of the chromosome when it comes to bins and trimming may be necessary before export.
    allLocations = do.call('rbind',pbapply::pblapply(cl = NULL, X = seq_along(ChromLengths), function(XX){

        #Start each chromosome at 1/2* largePeakWindow, so that the positions are the center of each bin.
        posList <- seq(startBin,ChromLengths[[XX]], largePeakWindow)
        ## check whether that bin is the index of a true peak, or in-between noise.
        data.frame(chr = rep(names(ChromLengths)[[XX]], length(posList)),
                    start = posList, end = posList + (startBin -1 ), mid = posList+startBin)

    }))
    allLocations$PeakID = NA
    allLocations$isPeak = c(1:length(allLocations$start) %in% peakLocation)
    allLocations$PeakID[allLocations$isPeak] = c(1:peakNumber)
    
    #########################################################################
    ## Generate the center of each peak, and the starting position for fragments within a peak
    #########################################################################

    ##Let's identify the peakNumber that each in-peak fragment will fall into.
    peakStarts <- do.call('rbind',
        pbapply::pblapply(cl = NULL,
                  X = seq_along(peakIntensities), function(XX) {

            if(peakIntensities[XX] != 0){
                data.frame(PeakID = rep(as.numeric(names(peakIntensities)[XX]), peakIntensities[XX]))
            }else{
                NULL
            }

            })
        )


    ## This peak center should also have some noise so that it generates a sharp peak, or multiple potential peaks. 
    ## randomly generate the types of peak centers. 
    message('Generating peak centers')
    ## for each peak, randomly choose a number of peak centers (peakCenters) at positions that are between 25% and 75% of the larger peak window. 
    ## The number of centers should be weighted by the peakIntensity - that is, peaks with more centers should have more fragments. 

    if(peakCenters != 1){

        numberOfCenters = rpois(n = length(peakIntensities), lambda = peakCenters)
        numberOfCenters[numberOfCenters == 0] = 1
        ## Now that we have a good distribution, let's sample from it, according to a weighted probability, where peaks with more fragments are more likely to have more centers. 
        lambda_choices = pbapply::pblapply(cl = NULL, X = seq_along(peakIntensities), function(XX)
                        c(peakIntensities[XX], sample(c(floor(largePeakWindow*0.25):ceiling(largePeakWindow*0.75)), size = numberOfCenters[XX]))
                    ) 
        ## Generate start positions according to multiple peakcenters
        message('Generating fragment start positions & peakIDs')
        allStarts <- unlist(pbapply::pblapply(cl = NULL, X = lambda_choices, simulateFragmentStarts))

    }else{
        allStarts = round(rpois(n = length(peakStarts$PeakID), largePeakWindow/3) + (largePeakWindow/2 - largePeakWindow/3))
    }

    peakStarts$FragStarts = allStarts 

    #########################################################################
    ## Merge the fragment information together, including fragment width, cell ID, peak ID, relative start, and bin ID. 
    #########################################################################

    # Now we need to merge the peakStart data frame (PeakID, Start positions, CellID) with the fragDF (fragment width, CellID, inPeak)
    # Because a given CellID can have multiple fragments fall within a peak (lambda2), we can't merge by CellID. 
    # Instead, we can simply arrange both peakStart and subFrags by CellID and then cbinds. Since both peakStart and subFrags have the same fragment number per cell, this will align,
    # and the peakIDs will be preserved. 
    subFrags <- dplyr::arrange(fragDF[fragDF$inPeak, ], PeakID)
    peakStarts <- dplyr::arrange(peakStarts, PeakID)
    #all(subFrags$PeakID == peakStarts$PeakID) All True
    subFrags$FragStarts = peakStarts$FragStarts
    
    ## Now we merge these fragments with peak locations to get their real locations. 
    peakFrags = dplyr::left_join(subFrags, allLocations[!is.na(allLocations$PeakID),],  by = 'PeakID')
    peakFrags <- dplyr::mutate(peakFrags, end = start + FragStarts + width, start = start + FragStarts)
    peakFrags <- dplyr::select(peakFrags, chr, start, end, CellID, PeakID, inPeak)

    #########################################################################
    ## Generate background fragments across background positions
    #########################################################################

    #### Now we need to place the non-peak fragments randomly across the genome. 
    #### let's get a complete list of all the possible locations by generating a non-peak GRanges and reducing it. 
    if(FRIP < 1){

        message('Placing background fragments into the genome')
        backgroundDF <- allLocations[is.na(allLocations$PeakID),]
        backgroundDF$backBin <- 1:dim(backgroundDF)[1]
        backFrags <- fragDF[!fragDF$inPeak,c(1:3)]
        ## Randomly select a position between 1 and the maximum peak window area. This will be the middle of the background tile. 
        backFrags$FragMid <- sample(1:largePeakWindow-1, size = sum(is.na(fragDF$PeakID)), replace = TRUE) 
        ## Sample with replacement, so that a background fragment could theoretically fall within the same background region. 
        backFrags$backBin <- sample(1:dim(backgroundDF)[1], size = sum(is.na(fragDF$PeakID)), replace = TRUE)

        backFrags <- dplyr::left_join(backFrags, backgroundDF, by = 'backBin')
        backFrags <- dplyr::mutate(backFrags, end = start + FragMid + ceiling(width/2), start = start + FragMid - floor(width/2))
        backFrags <- dplyr::select(backFrags, chr, start, end, CellID, PeakID, inPeak)


        #Combine background and foreground into one GRangesList
        allFrags <- GenomicRanges::makeGRangesFromDataFrame(rbind(peakFrags, backFrags), keep.extra.columns= TRUE)


    }else{
        allFrags <- GenomicRanges::makeGRangesFromDataFrame(peakFrags, keep.extra.columns= TRUE)
    }

    #Add Genome info and trim
    suppressWarnings(GenomicRanges::seqinfo(allFrags) <- GenomicRanges::seqinfo(Genome))
    allFrags <- GenomicRanges::trim(allFrags)

    allFrags$nCells = nCells
    allFrags$meanFragsPerCell = meanFragsPerCell
    allFrags$fragThreshold = fragThreshold
    allFrags$peakNumber = peakNumber
    allFrags$peakCenters = peakCenters
    allFrags$largePeakWindow = largePeakWindow
    allFrags$FRIP = FRIP
    allFrags$meanLengths = paste(meanLengths, collapse =', ')
    allFrags$lengthProbability= paste(lengthProbability, collapse =', ')

    peakSet <- GenomicRanges::makeGRangesFromDataFrame(allLocations[allLocations$isPeak,])

    return(list(peakSet, allFrags))
}


simulateFragmentStarts <- function(LambdaVector){

    centers_tmp =  LambdaVector[-1]
    #Generate weights for each center, so that fragments are not equally distributed across centers.
    if(length(centers_tmp) > 1){
        
        centWeights <- runif(n = length(centers_tmp), min = 0, max = 1)
        #Randomly select centers by the probability
        centerList <- sample(centers_tmp, size = LambdaVector[1], 
                            prob = centWeights/sum(centWeights), replace = TRUE)

        iterList = as.list(table(centerList))

        #Then do poisson sampling by the weight number of fragments for each center. 
        newStarts <- unlist(lapply(seq_along(iterList), function(x){

                rpois(iterList[[x]], lambda = as.numeric(names(iterList))[x])

        }))
    }else{

        newStarts = rpois(n = LambdaVector[1], lambda = LambdaVector[-1])

    }

    return(newStarts)
}


getSimulatedPeakSet <- function(peakNumber = 500, largePeakWindow = 500, Genome = BSgenome.Hsapiens.UCSC.hg38){

    #########################################################################
    ## Generate peak locations and background locations across the genome
    #########################################################################

    ### Now let's determine where our peaks are.
    ChromLengths = GenomeInfoDb::seqlengths(Genome)
    ChromLengths <- ChromLengths[!grepl('_', names(ChromLengths))]

    ### Generate all possible bins, and mark peak bins at evenly spaced intervals. 
    allBins = sum(ChromLengths)/largePeakWindow

    if(allBins <= peakNumber){
        stop('peakNumber is greater than or equal to the total number of larger open regions possible. Please change peakNumber, largePeakWindow, or both.')
    }

    binDistance = allBins/peakNumber
    peakLocation = seq(1, allBins, by = floor(binDistance))[1:peakNumber]
    startBin = 500

    browser()

    message('Placing large peak windows across genome.')
    #This arrangement may slightly depass the ends of the chromosome when it comes to bins and trimming may be necessary before export.
    allLocations = do.call('rbind',pbapply::pblapply(cl = NULL, X = seq_along(ChromLengths), function(XX){

        #Start each chromosome at 1/2* largePeakWindow, so that the positions are the center of each bin.
        posList <- seq(startBin,ChromLengths[[XX]], largePeakWindow)
        ## check whether that bin is the index of a true peak, or in-between noise.
        data.frame(chr = rep(names(ChromLengths)[[XX]], length(posList)),
                    start = posList, end = posList + (startBin -1), mid = posList+startBin)

    }))
    allLocations$PeakID = NA
    allLocations$isPeak = c(1:length(allLocations$start) %in% peakLocation)
    allLocations$PeakID[allLocations$isPeak] = c(1:peakNumber)

    locationGR <- GenomicRanges::makeGRangesFromDataFrame(allLocations[allLocations$isPeak,], keep.extra.columns = TRUE)
    
    return(locationGR)
}
