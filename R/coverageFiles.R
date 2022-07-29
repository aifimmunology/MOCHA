
### getCoverage turns sample/celltype-specific fragments lists into 
### sample-specific coverage files for each sample-celltype. 
# @popFrags - GRangesList of fragments for all sample/celltypes 
# @filterEmpty - True/False flag on whether or not to carry forward regions without coverage. 
# @numCores - number of cores to parallelize over

getCoverage <- function(popFrags, filterEmpty = FALSE, numCores = 1, verbose = FALSE){
  #Extract cell counts for each group   
  groups <- gsub("__.*","",names(popFrags))
  Norm <- gsub(".*__","", names(popFrags))
  names(Norm) = groups
     
	#Summarize the coverage ove the region window at a single basepair resolution
	tmpCounts <- mclapply(seq_along(groups), function(x){
	    if(verbose){print(paste("Counting", names(popFrags)[x], sep=" "))}
	  
	    Num <- as.numeric(Norm[groups[x]])
	    tmp <- plyranges::compute_coverage(popFrags[[x]]) %>% plyranges::mutate(score = score/Num) 
	    if(filterEmpty) {plyranges::filter(tmp, score > 0)}else{tmp}
	}, mc.cores = numCores)
	names(tmpCounts) = groups
    
   return(tmpCounts)
}
    


######## getSpecificCoverage: Function that takes in a GRangesList of coverage per sample/celltype and finds 
## the average coverage intensity for just specific regions. 
## @covFiles: GRangesList of coverage for each sample
## @regions: Regions to count intensities over, must be non-overlapping and non-adjacent ( > 1 bp apart).  
## @numCores: number of cores to parallelize over. 
getSpecificCoverage <- function(covFiles, regions, numCores=1){

	counts <- mclapply(covFiles, function(x){ 
                	x %>% plyranges::mutate(NewScore = score) %>%
                    plyranges::join_overlap_intersect(regions) %>%
                    plyranges::mutate(WeightedScore = NewScore*width(.)) %>%
                             plyranges::reduce_ranges(score = mean(WeightedScore)) 
    
    	}, mc.cores = numCores)


	return(counts)
}
    