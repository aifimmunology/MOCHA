
### getCoverage turns sample/celltype-specific fragments lists into 
### sample-specific coverage files for each sample-celltype. 
# @popFrags - GRangesList of fragments for all sample/celltypes 
# @filterEmpty - True/False flag on whether or not to carry forward regions without coverage. 
# @numCores - number of cores to parallelize over

getCoverage <- function(popFrags, normFactor, filterEmpty = FALSE, numCores = 1, TxDb = TxDb.Hsapiens.UCSC.hg38.refGene, verbose = FALSE){

  if(length(normFactor) == 1){

	normFactor = rep(normFactor, length(popFrags))

  }  

  #Summarize the coverage over the region window at a single basepair resolution
  tmpCounts <- parallel::mclapply(seq_along(popFrags), function(x){

	    if(verbose){print(paste("Counting", names(popFrags)[x], sep=" "))}
	  
	    Num <- normFactor[x]
	    tmp <- plyranges::compute_coverage(popFrags[[x]]) %>% plyranges::mutate(score = score/Num) 

	    GenomeInfoDb::seqinfo(tmp) <- GenomeInfoDb::seqinfo(TxDb)[GenomicRanges::seqnames(GenomeInfoDb::seqinfo(tmp))]

	    if(filterEmpty) {plyranges::filter(tmp, score > 0)}else{tmp}

	}, mc.cores = numCores)

   names(tmpCounts) = names(popFrags)
    
   return(tmpCounts)
}
    


######## getSpecificCoverage: Function that takes in a GRangesList of coverage per sample/celltype and finds 
## the average coverage intensity for just specific regions. 
## @covFiles: GRangesList of coverage for each sample
## @regions: Regions to count intensities over, must be non-overlapping and non-adjacent ( > 1 bp apart).  
## @numCores: number of cores to parallelize over. 
getSpecificCoverage <- function(covFiles, regions, numCores=1){

	counts <- parallel::mclapply(covFiles, function(x){ 
                	x %>% plyranges::mutate(NewScore = score) %>%
                    plyranges::join_overlap_intersect(regions) %>%
                    plyranges::mutate(WeightedScore = NewScore*width(.)) %>%
                             plyranges::reduce_ranges(score = mean(WeightedScore)) 
    
    	}, mc.cores = numCores)


	return(counts)
}
    