#' @title \code{getPopFrags}
#'
#' @description \code{getPopFrags} is an R helper function that extracts fragments  
#'              for a cell population from an archR Project

#' @param ArchRProj an ArchR Project 
#' @param metaColumn a string indicating which column to use to extract cells
#' @param cellSubsets vector of strings. Cell subsets for which to call peaks. Optional, if cellSubsets='ALL', then peak calling is done on all cell populations in the ArchR project metadata
#' @param region a chromosomal location. If NULL, it pull from the entire genome. 
#' @param numCores integer number of cores to parallelize on
#' @param NormMethod = string indicating which normalization method to use. For scMACS modeling
#' we require the nFrags normalization
#' @param blacklist  a granges object indicating the blacklist region
#' @param overlapList integer indicating acceptable overlap with blacklist regions in bp 
#'
#'
#' @return a list of fragments with a given normalization factor.
#'
#' @export


getPopFrags <- function(ArchRProj, metaColumn,cellSubsets = 'ALL' , region = NULL, numCores = 1, 
                        NormMethod = "NFrags", blackList = NULL, overlapList = 50){
    
    #Extract metadata
    metadf <- getCellColData(ArchRProj)
    
    if(!any(colnames(metadf) %in% metaColumn)){
    
        stop('Error: Metadata column not found')
    }
    # Determine if all or just some cell types need to be extracted.
    # Calculate the number of cells in those groups.
    # Calculate the total nFrags for those groups. 
    if(all(cellSubsets=='ALL')){
      tmp <- metadf[,metaColumn]
      tmp <- tmp[!is.na(tmp)]
      cellCounts = table(tmp)
      cellPopulations=names(cellCounts)
      Norm = as.data.frame(metadf[,c(metaColumn,"nFrags")]) %>% 
                    dplyr::group_by(get(metaColumn))  %>% 
                    dplyr::summarize(nFrags = sum(nFrags))
      
    } else{
      tmp <- table(metadf[,metaColumn])
      cellPopulations=cellSubsets[!is.na(cellSubsets)][order(cellSubsets)]
      cellCounts = tmp[cellPopulations] 
      Norm = as.data.frame(metadf[,c(metaColumn,"nFrags")]) %>% 
                    dplyr::filter(get(metaColumn) %in% cellSubsets) %>% 
                    dplyr::group_by(get(metaColumn)) %>% 
                    dplyr::summarize(nFrags = sum(nFrags)) 
    }

    #Check to make sure cellCounts were calculated correctly.
    if(any(is.na(cellCounts))){
      
      stop('Wrong cell subset input. Please verify that cell subset exists within the correct metadata column.')
      
    }
    #Check to make sure nFrag normalization matches up to cell table. 
    if(!all(Norm[,1] == cellPopulations)){
        print(Norm[,1])
        stop("Issue with nFrag normalization.")
    }
    
    #Extract fragments in an efficient manner
    allArrows <- getArrowFiles(ArchRProj)

    #Only extract from arrows of interest
    arrowSearch <- paste(unique(metadf$Sample[which(metadf[,metaColumn] %in% names(cellCounts))]), collapse = "|")
    arrows <- allArrows[grepl(arrowSearch, allArrows)]

    if(length(arrows) == 0 | length(cellCounts) == 0){

	        stop("No arrows or no cell groups: Check your cell subset or ArchR Project input.")

     }    

    if(is.null(region)){

	        fragsList<-  mclapply(seq_along(arrows), function(x){
                    getFragmentsFromArrow(arrows[x])
            }, mc.cores = numCores)
            

    }else{
      
       if(is.character(region)){
         
           regionGRanges <- StringsToGRanges(region)
           
        }else if(class(region)[1] == "GRanges"){
    
            regionGRanges = region
            
        }else if(!is.null(region)){
          
          stop('Wrong region input type. Input must either be a string, or a GRanges location.')
          
        }

      chrom = regionGRanges %>% as.data.frame() %>% dplyr::select(seqnames)
	  fragsList<-  mclapply(seq_along(arrows), function(x){
			      
                  tmp1 <- getFragmentsFromArrow(arrows[x],chr = as.character(chrom[,1])) %>%
							plyranges::join_overlap_intersect(regionGRanges)
                  if(is.null(blackList)){
                      tmp1
                    }else if(class(blackList)[1] == "GRanges"){
                      
                      tmp1 %>% plyranges::filter_by_non_overlaps(blackList, minoverlap = overlapList)
                      
                    }else{
                      
                          stop('Error: Wrong format for blackList region. Please provide GRanges')
                    }
                      
            }, mc.cores = numCores)
    }
    
    #This calculates total fragments per sample. 
    totalsampleFrags = unlist(fragsList, length)
    
    #From scMACS - sorts cell barcodes by population
    barcodes_by_cell_pop <- lapply(cellPopulations, function(x)
        row.names(metadf)[which(metadf[,metaColumn]==x)]
           )
    
    # Add normalization factor. 
    if(tolower(NormMethod) == "ncells"){
        
        names(barcodes_by_cell_pop) = paste(gsub(" |_",".",cellPopulations), cellCounts/1000, sep = "_")
        
    }else if(tolower(NormMethod) == "nfrags"){
        
        names(barcodes_by_cell_pop) = paste(gsub(" |_",".",cellPopulations), Norm$nFrags/10^6, sep = "_")
        
    }else if(tolower(NormMethod) == "median"){
        
        names(barcodes_by_cell_pop) = paste(gsub(" |_",".",cellPopulations), "Median", sep = "_")
        
    }else if(tolower(NormMethod) == "medmax"){
        
        names(barcodes_by_cell_pop) = paste(gsub(" |_",".",cellPopulations), "MedianMax", sep = "_")
        
    }else if(tolower(NormMethod) == "TotalSampleByNCells"){
        
        names(barcodes_by_cell_pop) = paste(gsub(" |_",".",cellPopulations), totalsampleFrags/10^6*cellCounts/1000, sep = "_")
        
    }else{
    
        stop("Error: Incorrect NormMethod given.")
    }

    
    #From scMACS - Function to sort fragments by populations based on cell barcode lists
    subset_Frag <- function(cellNames, tmp){
            tmp_dt = data.table::as.data.table(tmp)
            idx <- which(tmp_dt$RG %in% cellNames)
            tmp[ idx] 
        }
    ## Identify which subset of arrows the population can be found it. 
    ## Speeds up sample-specific population extraction
        
    fragsListIndex <- lapply(cellPopulations, function(x)
        names(arrows) %in% unique(metadf$Sample[which(metadf[,metaColumn]==x)])
           )
       
    #Sort fragments into a list by cell population
        
    popFrags <- mclapply(seq_along(barcodes_by_cell_pop),function(x){
                    print(paste('Sorting ', names(barcodes_by_cell_pop)[x], sep = ""))
                    if(sum(fragsListIndex[[x]]) > 1){
                       tmp <- lapply(which(fragsListIndex[[x]]) , function(y) {
                            subset_Frag(barcodes_by_cell_pop[[x]], fragsList[[y]])
                        }) 
                      stack(as(tmp,'GRangesList'))
                    }else{
                       subset_Frag(barcodes_by_cell_pop[[x]], fragsList[[which(fragsListIndex[[x]])]])
                    }
               }, mc.cores=numCores)
    names(popFrags) <- names(barcodes_by_cell_pop)
    return(popFrags)
    
}