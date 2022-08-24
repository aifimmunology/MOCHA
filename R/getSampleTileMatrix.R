#' @title \code{getSampleTileMatrix}
#'
#' @description \code{getSampleTileMatrix} takes the output of peak calling with 
#'   callOpenTiles and creates sample-tile matrices containing the signal intensity
#'   at each tile. 
#'
#'
#' @param tileResults a MultiAssayExperiment returned by callOpenTiles containing
#'   containing peak calling results.
#' @param cellPopulations vector of strings. Cell subsets in TileResults for which to generate sample-tile matrices. This list of group names must be identical to names that appear in
#'   the ArchRProject metadata.  If cellPopulations='ALL', then peak
#'   calling is done on all cell populations in the ArchR project metadata.
#'   Default is 'ALL'.
#' @param groupColumn Optional, the column containing sample group labels for determining consensus tiles within sample groups. Default is NULL, all samples will be used for determining consensus tiles.
#' @param threshold Threshold for consensus tiles, the minimum % of samples (within a sample group) that a peak must be called in to be retained.
#' @param join The join method to combine consensus tiles across sample groups. Can be "union" (default) or "intersect". 
#' @param NAtoZero Boolean, set to TRUE to convert NA intensities from peak calling (tiles with too low signal) to zeroes. Optional, default is FALSE.
#' @param log2Intensity Boolean, set to TRUE to return the log2 of the sample-tile intensity matrix. Optional, default is FALSE.
#' @param numCores Optional, the number of cores to use with multiprocessing. Default is 1.
#' 
#' @return SampleTileMatrices a MultiAssayExperiment containing a sample-tile intensity matrix 
#'   for each cell population
#'   
#' @details The technical details of the algorithm are found in XX.
#'
#' @references XX
#'
#' @export

getSampleTileMatrix <- function(tileResults,
                                cellPopulations = "ALL",
                                groupColumn = NULL,
                                threshold = 0.2,
                                join = "union",
                                NAtoZero = FALSE,
                                log2Intensity = TRUE,
                                numCores = 1,
                                verbose = FALSE
                               ) {

  if (!(join %in% c("union", "intersect"))){
    stop("`join` must be either 'union' or 'intersect'")
  }
  
  if(class(tileResults)[1] != 'MultiAssayExperiment'){
    stop('tileResults is not a MultiAssayExperiment')  
  }
  
  # Any column can be used to group samples
  # Note that these are case-sensitive
  sampleData <- MultiAssayExperiment::colData(tileResults)
  validGroups <- colnames(sampleData)
  if (!is.null(groupColumn)){
    if (!(groupColumn %in% validGroups)){
      stop("`groupColumn` not found in the column data of tileResults.")
    }
  }
  
  if(length(cellPopulations) == 1 & any(cellPopulations == 'ALL')){
  
      subTileResults <- tileResults
      cellPopulations <- names(tileResults)
  
  }else{
    
    if (all(cellPopulations %in% names(tileResults))){
       subTileResults <- tileResults[names(tileResults) %in% cellPopulations]
    }else{
      stop(paste(
          "All of `cellPopulations` must present in tileResults.",
          "Check `names(tileResults)` for possible cell populations."
      ))
    }
    
  }
  

  if(verbose){
    message(stringr::str_interp("Extracting consensus tile set for each population:  ${paste(cellPopulations, collapse=', ')} "))
  }

  tilesByCellPop <- parallel::mclapply(MultiAssayExperiment::experiments(subTileResults), function(x){
    
    # Get consensus tiles for this cell population for filtering
    scMACS:::singlePopulationConsensusTiles(
      x,
      sampleData,
      threshold,
      groupColumn,
      join
    )

  }, mc.cores = numCores)

  allTiles <- sort(unique(do.call('c',tilesByCellPop)))

  if(verbose){
    message(stringr::str_interp("Generating sample-tile matrix across all populations: ${paste(cellPopulations, collapse=', ')} "))
  }

  # consensusTiles is used to  extract rows (tiles) from this matrix
  sampleTileIntensityMatList <- parallel::mclapply(MultiAssayExperiment::experiments(tileResults), function(x){
    scMACS:::singlePopulationSampleTileMatrix(
      x,
      allTiles,
      NAtoZero,
      log2Intensity
    )
  }, mc.cores = numCores)

  tilePresence <- lapply(tilesByCellPop, function(x) (allTiles %in% x)) %>% do.call('cbind', .) %>% as.data.frame()
  allTilesGR <- scMACS::StringsToGRanges(allTiles)
  mcols(allTilesGR) <- tilePresence

  results <- SummarizedExperiment::SummarizedExperiment(
    sampleTileIntensityMatList,
    rowRanges = allTilesGR,
    colData = sampleData,
    metadata = append(
      list('Log2Intensity' = log2Intensity, 'NAtoZero' = NAtoZero ),
      MultiAssayExperiment::metadata(tileResults)
    )
  )
  return(results)
}


##### function for annotating the tiles for the Sample Tile Matrix.
## Returns a sample tile matrix with the transcript Db annotated and tile annotations. 


annotateTiles <- function(obj,
			  TxDb = NULL,
			   Org = NULL, 
                          promoterRegion = c(2000, 100)){

	if(class(obj)[1] == "RangedSummarizedExperiment" & is.null(TxDb) & is.null(Org)){
	
		if(!all(c('TxDb', 'Org') %in% names(S4Vectors::metadata(obj)))){
			
			stop('Error: Wrong input. Obj must either be an RangedSummarizedExperiment with a TxDb and Org, or a GRanges obj')

		}
		tileGRanges <- SummarizedExperiment::rowRanges(obj)
		TxDb = AnnotationDbi::loadDb(S4Vectors::metadata(obj)$TxDb)
		Org = AnnotationDbi::loadDb(S4Vectors::metadata(obj)$Org)

	}else if(class(tileList)[[1]] == 'GRanges' & !is.null(TxDb) & !is.null(Org)){

		tileGRanges = obj		

	}else{

		stop('Error: Wrong inputs. Verify proper inputs for obj, TxDb, and/or Org')

	}

	txList <- suppressWarnings(GenomicFeatures::transcriptsBy(TxDb, by = ('gene')))
        names(txList) <- suppressWarnings(AnnotationDbi::mapIds(Org, names(txList), "SYMBOL", "ENTREZID"))
	exonList <-  suppressWarnings(ensembldb::exonsBy(TxDb, by = "gene"))
	#Same TxDb, same gene names in same order. 
	names(exonList) <-  names(txList) 

	txs <- stack(txList) %>% GenomicRanges::trim() %>% S4Vectors::unique(.)
	exonSet <- stack(exonList) 

    	promoterSet <- stack(txList) %>% GenomicRanges::trim(.) %>% S4Vectors::unique(.) %>% 
		suppressWarnings(GenomicRanges::promoters(., upstream = promoterRegion[1], downstream = promoterRegion[2]))

	getOverlapNameList <- function(rowTiles, annotGR){
    
    		overlapGroup <- findOverlaps(rowTiles, annotGR) %>% as.data.frame() 
    		overlapGroup$Genes = as.character(annotGR$name[overlapGroup$subjectHits])
    		last <- overlapGroup %>% dplyr::group_by(subjectHits) %>% 
                dplyr::summarize(Genes = paste(unique(Genes), collapse =", "))
    		return(last)

	}

	txs_overlaps <- getOverlapNameList(tileGRanges, txs)
	exon_overlaps <- getOverlapNameList(tileGRanges, exonSet) 
	promo_overlaps <- getOverlapNameList(tileGRanges, promoterSet) 

	tileType <- as.data.frame(mcols(rowTiles)) %>% dplyr::mutate(Index = 1:nrow(.)) %>%   
            dplyr::mutate(Type = dplyr::case_when(Index %in% promo_overlaps$subjectHits ~ 'Promoter',
               Index %in% exon_overlaps$subjectHits ~ 'Exonic',
                Index %in% txs_overlaps$subjectHits ~ 'Intronic',
                TRUE ~ 'Distal')) %>% 
            dplyr::left_join(promo_overlaps, by = c('Index' = 'subjectHits')) %>%
            dplyr::rename('Promo' = Genes)%>% 
            dplyr::left_join(exon_overlaps, by = c('Index' = 'subjectHits')) %>%
            dplyr::rename('Exons' = Genes)%>% 
            dplyr::left_join(txs_overlaps, by = c('Index' = 'subjectHits')) %>%
            dplyr::rename('Txs' = Genes)%>% 
            dplyr::mutate(Genes = ifelse(Type == 'Promoter',Promo, NA)) %>%
            dplyr::mutate(Genes = ifelse(Type == 'Exonic',Exons, Genes)) %>%
            dplyr::mutate(Genes = ifelse(Type == 'Intronic',Txs, Genes)) 

	tileGRanges$tileType = tileType$Type
	tileGRanges$Gene = tileType$Genes

	#If input was as Ranged SE, then edit the rowRanges for the SE and return it. 
	#Else, return the annotated tile GRanges object. 
	if(class(obj)[1] == "RangedSummarizedExperiment"){

		SummarizedExperiment::rowRanges(obj) <- tileGRanges
		return(obj)		

	}else{
		return(tileGRanges)

	}
	
}
