
## tile1 and tile2 should be the numeric indices of the TSAM to test

testCoAccessibility <- function(STObj, tile1, tile2, numCores = 1, ZI = TRUE, backNumber = 100, highMem = FALSE, verbose = TRUE){

    if(length(tile1) != length(tile2)){

        stop('tile1 and tile2 must be the same length.')

    }

    fullObj <- getBackGroundObj(STObj)

    backPeaks <- chromVAR::getBackgroundPeaks(fullObj)

    if(is.character(tile1) & is.character(tile2)){

        nTile1 = match(tile1, rownames(fullObj))
        nTile2 = match(tile1, rownames(fullObj))
    }else if(is.numeric(tile1) & is.numeric(tile2)){

        nTile1 = tile1
        nTile2 = tile2

        tile1 = rownames(fullObj)[nTile1]
        tile2 = rownames(fullObj)[nTile2]

    }else{stop('tile1 and tile 2 must both be either numbers (indices) or strings')}

    if(backNumber > 2450){
        backNumber = 2450
        warning('backNumber too high. Reset to 2500')

    }else if(backNumber <= 10){
      
      stop('backNumber too low (<=10). We recommend at least 100.')

    }

    cl <- parallel::makeCluster(numCores)
    parallel::clusterExport(cl, varlist = c('nTile1', 'nTile2','backPeaks'), envir = environment())
    
    message('Finding background peak pairs.')

    backgroundCombos <- pbapply::pblapply(seq_along(nTile1), function(x){
        tmpMat <- expand.grid(backPeaks[nTile1[x],], backPeaks[nTile2[x],])
        tmpMat <- tmpMat[tmpMat[,1] != tmpMat[,2],] 
        tmpMat[sample.int(dim(tmpMat)[1], backNumber),] 
    }, cl = cl)

    parallel::stopCluster()

    gc()

    cl <- parallel::makeCluster(numCores)

    accMat <- SummarizedExperiment::assays(fullObj)[[1]]
    
    ## Test original pairs of locations
    combPairs <- data.frame(tile1, tile2)
    subAccMat <- accMat[unique(c(nTile1, nTile2)),]

    message('Identifying foreground.')
    parallel::clusterExport(cl, varlist = c('subAccMat', "combPairs"), envir = environment())
    foreGround <- runCoAccessibility(subAccMat, combPairs, ZI, verbose, cl)
    parallel::stopCluster()

    gc()

    ## Now we need to test the background set

    message('Identifying background correlations.')

    if(highMem){

      allBackCombos <- do.call('rbind',.)

      uniqueBackCombos <- unique(allBackCombos)

      combPairs <- data.frame(Tile1 = rownames(accMat)[uniqueBackCombos[,1]],
                                  Tile2 = rownames(accMat)[uniqueBackCombos[,2]])
      
      subAccMat <- accMat[unique(c(uniqueBackCombos[,1],uniqueBackCombos[,2])),]

      cl <- paralle::makeCluster(numCores)
      message(paste('Generating Background correlations for all', length(tile1), 'pairs', sep = " "))
      parallel::clusterExport(cl, varlist = c('subAccMat', "combPairs"), envir = environment())

      uniquebackGround <- runCoAccessibility(subAccMat, combPairs, ZI, verbose, cl)

      backGround <- dplyr::left_join(allBackCombos, ) %>%  group_by(row_number() %/% backNumber) %>% group_map(~ .x) 

      parallel::stopCluster()

    }else{

      ###low memory implementation

      backGround <- list()

      for (i in 1:length(backgroundCombos)) {

          combPairs <- data.frame(Tile1 = rownames(accMat)[backgroundCombos[[i]][,1]],
                                  Tile2 = rownames(accMat)[backgroundCombos[[i]][,2]])
      
          subAccMat <- accMat[unique(c(backgroundCombos[[i]][,1],backgroundCombos[[i]][,2])),]

          cl <- parallel::makeCluster(numCores)
          message(paste('Generating Background correlations for', i, 'of', length(tile1), 'pairs', sep = " "))
          parallel::clusterExport(cl, varlist = c('subAccMat', "combPairs"), envir = environment())

          tmp_background <- runCoAccessibility(subAccMat, combPairs, ZI, verbose, cl)
          
          parallel::stopCluster(cl)

          gc()

          backGround <- append(backGround, tmp_background)
      }

    }


    cl <- parallel::makeCluster(numCores)
    parallel::clusterExport(cl, varlist = c('foreGround','backGround'), envir = environment())

    pValues <- pbapply::pblapply(1:length(backGround), function(x){
    
        cor1 <- foreGround$Correlation[x]
        
        if(cor1 >= 0){
            
            sum(cor1 > backGround[[x]]$Correlation)/length(backGround[[x]]$Correlation)
            
        }else if(cor1 < 0){
            
            sum(cor1 < backGround[[x]]$Correlation)/length(backGround[[x]]$Correlation)
            
        }
        
    },cl =cl) %>% unlist()
    parallel::stopCluster(cl)

    gc()

    foreGround$pValues <- pValues

    return(foreGround)

}


testCoAccessibility2 <- function(STObj, tile1, tile2, numCores = 1, ZI = TRUE, backNumber = 1000, highMem = FALSE, verbose = TRUE){

    if(length(tile1) != length(tile2)){

        stop('tile1 and tile2 must be the same length.')

    }

    fullObj <- getBackGroundObj(STObj)

    #backPeaks <- chromVAR::getBackgroundPeaks(fullObj)

    if(is.character(tile1) & is.character(tile2)){

        nTile1 = match(tile1, rownames(fullObj))
        nTile2 = match(tile1, rownames(fullObj))
    }else if(is.numeric(tile1) & is.numeric(tile2)){

        nTile1 = tile1
        nTile2 = tile2

        tile1 = rownames(fullObj)[nTile1]
        tile2 = rownames(fullObj)[nTile2]

    }else{stop('tile1 and tile 2 must both be either numbers (indices) or strings')}

    if(backNumber >=  length(rownames(fullObj)) - c(length(tile1) + length(tile2))){
        backNumber = length(rownames(fullObj)) - c(length(tile1) + length(tile2))
        warning('backNumber too high. Reset to all background combinations.')

    }else if(backNumber <= 10){
      
      stop('backNumber too low (<=10). We recommend 1000.')

    }

    cl <- parallel::makeCluster(numCores)

    accMat <- SummarizedExperiment::assays(fullObj)[[1]]
    
    ## Test original pairs of locations
    combPairs <- data.frame(tile1, tile2)
    subAccMat <- accMat[unique(c(nTile1, nTile2)),]

    message('Identifying foreground.')
    parallel::clusterExport(cl, varlist = c('subAccMat', "combPairs"), envir = environment())
    foreGround <- runCoAccessibility(subAccMat, combPairs, ZI, verbose, cl)
    parallel::stopCluster(cl)

    gc()
    
    message('Finding background peak pairs.')

    backGroundTiles <- rownames(fullObj)[!rownames(fullObj) %in% c(tile1, tile2)]
    
    backgroundCombos <- data.frame(Tile1 = sample(backGroundTiles, backNumber),
                          Tile2 = sample(backGroundTiles, backNumber))
    backgroundCombos <- backgroundCombos[backgroundCombos[,1] != backgroundCombos[,2],] 

    ## Now we need to test the background set

    message('Identifying background correlations.')

    cl <- parallel::makeCluster(numCores)
    subAccMatB <- accMat[c(backgroundCombos$Tile1, backgroundCombos$nTile2),]
    parallel::clusterExport(cl, varlist = c('backgroundCombos','subAccMatB'), envir = environment())
    backGround <- runCoAccessibility(subAccMatB, backgroundCombos, ZI, verbose, cl)
    parallel::stopCluster()

    cl <- parallel::makeCluster(numCores)
    parallel::clusterExport(cl, varlist = c('backGround','foreGround'), envir = environment())

    pValues <- pbapply::pblapply(1:length(tile1), function(x){
    
        cor1 <- foreGround$Correlation[x]
        
        if(cor1 >= 0){
            
            sum(cor1 > backGround$Correlation)/length(backGround$Correlation)
            
        }else if(cor1 < 0){
            
            sum(cor1 < backGround$Correlation)/length(backGround$Correlation)
            
        }
        
    },cl = cl) %>% unlist()
    parallel::stopCluster(cl)

    gc()

    foreGround$pValues <- pValues

    return(foreGround)

}


runCoAccessibility <- function(accMat, pairs, ZI = TRUE, verbose = TRUE, numCores = 1){

    zero_inflated_spearman <- unlist(pbapply::pblapply(1:dim(pairs)[1],
      function(x) {
        MOCHA:::weightedZISpearman(
          x = accMat[pairs[x, 1], ],
          y = accMat[pairs[x, 2], ],
          verbose = verbose,
          ZI = ZI
        )
      },
      cl = numCores
    ))

    # Create zero-inflated correlation matrix from correlation values
    zi_spear_mat_tmp <- data.table::data.table(
      Correlation = zero_inflated_spearman,
      Tile1 = pairs[, 1],
      Tile2 = pairs[, 2]
    )

    return(zi_spear_mat_tmp)

}


getBackGroundObj <- function(STObj, NAtoZero = TRUE){

  # Extract all the Sample-Tile Matrices for each cell type
  temp <- SummarizedExperiment::assays(STObj)

  meta1 <- SummarizedExperiment::colData(STObj)
  
  # Let's generate a new assay, that will contain the
  # the intensity for a given cell, as well as the
  # median intensity per sample-tile for all other cell types (i.e. the background)

   	  

  newAssays <- list(do.call('cbind', temp))
  newSamplesNames <- unlist(lapply(names(temp), function(x){
                            paste(x, colnames(STObj), sep = "__") %>% gsub(" ","",.)
            }))

  names(newAssays) <- 'counts'
  colnames(newAssays[[1]]) = newSamplesNames

  if(NAtoZero){  
    newAssays[[1]][is.na(newAssays[[1]])] = 0 
  }


  allSampleData <- do.call('rbind', lapply(names(temp), function(x){

		tmp_meta = meta1
		 tmp_meta$Sample <- paste(x,tmp_meta$Sample, sep = "__") %>% gsub(" ","",.)
                 tmp_meta$CellType = rep(x, dim(tmp_meta)[1])
		 rownames(tmp_meta) <- tmp_meta$Sample
		 tmp_meta
	}
  ))

  cellCounts <- as.data.frame(metadata(STObj)$CellCounts) %>%
    dplyr::mutate(Sample = gsub(" ", "", paste(cellTypeLabelList,Var1, sep = "__"))) %>%
    dplyr::select(Sample, Freq)
  
  allSampleData  <- dplyr::left_join(as.data.frame(allSampleData), cellCounts, by = 'Sample')


  allRanges <- SummarizedExperiment::rowRanges(STObj)
  for (i in names(temp)) {
    mcols(allRanges)[, i] <- rep(TRUE, length(allRanges))
  }

  #return(list(newAssays, allSampleData, allRanges))

  newObj <- SummarizedExperiment(
    assays = newAssays,
    colData =  allSampleData,
    rowRanges = allRanges,
    metadata = metadata(STObj)
  )

  newObj <- chromVAR::addGCBias(newObj, genome = metadata(STObj)$Genome)

   if(any(is.na(SummarizedExperiment::rowData(newObj)$bias))){

                    naList <- is.na(SummarizedExperiment::rowData(newObj)$bias)
                    message(paste(sum(naList), "NaNs found within GC Bias", sep =" "))

                    SummarizedExperiment::rowData(newObj)$bias[which(naList)] = mean(rowData(newObj)$bias, na.rm = TRUE)

                }

  return(newObj)

}