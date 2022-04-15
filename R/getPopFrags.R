#' This function outputs seperate out list of fragments per cell subset. The name for the list is in the format: populationName_cellCount
#' @export

getPopFrags <- function(ArchRProj, metaColumn,cellSubsets = 'ALL' , region = NULL, numCores){
    
    #Extract metadata
    metadf <- getCellColData(ArchRProj)

    #Determine if all or just some cell types need to be extracted.
    if(all(cellSubsets=='ALL')){
       tmp <- metadf[,metaColumn]
       tmp <- tmp[!is.na(tmp)]
       cellCounts = table(tmp)
       cellPopulations=names(cellCounts)
      
    } else{
       tmp <- table(metadf[,metaColumn])
       cellPopulations=cellSubsets[!is.na(cellSubsets)]
       cellCounts = tmp[cellPopulations]
    }

    #Extract fragments in an efficient manner
    allArrows <- getArrowFiles(ArchRProj)

    #Only extract from arrows of interest
    arrowSearch <- paste(unique(metadf$Sample[which(metadf[,metaColumn] %in% names(cellCounts))]), collapse = "|")
    
    arrows <- allArrows[grepl(arrowSearch, allArrows)]

    if(length(arrows) == 0){

	stop("No arrows: Check your cell subset or ArchR Project input.")

     }    

    if(is.null(region)){

	fragsList<-  mclapply(seq_along(arrows), function(x){
                    getFragmentsFromArrow(arrows[x])
            }, mc.cores = numCores)

    }else{
	
	chrom <- gsub(":.*","",region)
	startSite <- gsub(".*:","",region) %>% gsub("-.*", "",.) %>% as.numeric()
	endSite <- gsub(".*-","",region) %>% as.numeric()

	regionGRanges <- GRanges(seqnames = chrom, ranges = IRanges(start = startSite,end = endSite), strand = "*")

	fragsList<-  mclapply(seq_along(arrows), function(x){
			getFragmentsFromArrow(arrows[x],chr = chrom) %>%
							plyranges::join_overlap_intersect(regionGRanges)
            }, mc.cores = numCores)
    }
    
    #From scMACS - sorts cell barcodes by population
    barcodes_by_cell_pop <- lapply(cellPopulations, function(x)
        row.names(metadf)[which(metadf[,metaColumn]==x)]
           )
    names(barcodes_by_cell_pop) = paste(gsub(" |_",".",cellPopulations), cellCounts, sep = "_")

    
    #From scMACS - Function to sort fragments by populations based on cell barcode lists
    subset_Frag <- function(cellNames, tmp){
            tmp_df <- GenomicRanges::as.data.frame(tmp)
            idx <- which(tmp_df$RG %in% cellNames)
            tmp[ idx] 
        }
    #Sort fragments into a list by cell population
    popFrags <- lapply(seq_along(barcodes_by_cell_pop),function(x){
                      tmp <- mclapply(seq_along(fragsList) , function(y) {
                            subset_Frag(barcodes_by_cell_pop[[x]], fragsList[[y]])
                        }, mc.cores = numCores)
		      tmp_ranges <- tmp[unlist(lapply(tmp, function(x) {length(x) >0 }))]
                      if(length(tmp_ranges) > 1){
                      	plyranges::bind_ranges(tmp_ranges)
		      }else{
			tmp_ranges
		      }
        })
    names(popFrags) <- names(barcodes_by_cell_pop)
    popFrags2 <- popFrags[unlist(lapply(popFrags, length)) > 0]
    return(popFrags2)
    
}
