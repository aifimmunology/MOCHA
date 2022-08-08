##########
### Functions for motif analysis (enrichment and chromVar)



### Identify motifs within peakset
## @SE_Object - Any SummarizedExperiment-class object for which rowRanges works. 
## @p.cutoff - p.value cutoff to pass to motifmatchr
## @w - 

getMotifSet <- function(SE_Object, pwms = human_pwms_v2, 
                        genome =  BSgenome.Hsapiens.UCSC.hg38, p.cutoff =  5e-05, w = 7){

    TotalPeakSet <- rowRanges(SE_Object)
    motif_ix <- motifmatchr::matchMotifs(pwms = pwms, 
                                         TotalPeakSet, 
                                         genome = genome, 
                                      out = "positions", w = 15)
    names(motif_ix) <- sub("_D_.*|_I_.*","", names(motif_ix)) %>% sub("_I$|_D$","", .) %>%
                    sub(".*_LINE","",.) %>% sub(".*_","",.)


    return(motif_ix)

}


## runChromVar is a wrapper for chromVAR from scMACS
## Obj could be a ragged experiment or a sumarized experiment
## motifGRangesList is a list of all motif positions in a GRangesList format
## genome is the reference genome we use


runChromVar <- function(Obj, motifGRangesList,
                       genome =  BSgenome.Hsapiens.UCSC.hg38){
    
    if(class(Obj) == "RaggedSummarizedExperiment"){
        
        Obj1 <- RaggedExperiment::compactSummarizedExperiment(RepTiles, i = 'TotalIntensity')
        tmp <- SummarizedExperiment::assays(Obj1)
        tmp[[1]][is.na(tmp[[1]])] = 0
        names(tmp) <- "counts"
        assays(Obj1) <- tmp
        
    }else if(class(Obj) == "RangedSummarizedExperiment"){
        
        if( !names(assays(Obj)) %in% 'counts'){
        
            stop("Error: Assay names must include an assay named 'counts'")
            
        }else{
            
            Obj1 = Obj
        }
        
    }else{
    
        stop('Error: Wrong Input Object. Must be either a RaggedExperiment or a RangedSummarizedExperiment')
    }
    
    CisbpAnno <- chromVAR::getAnnotations(motifGRangesList, rowRanges = rowRanges(Obj1))

    Obj1 <- chromVAR::addGCBias(Obj1, genome = genome)
    
    dev <- chromVAR::computeDeviations(object = Obj1, 
                         annotations = CisbpAnno)
    
    return(dev)
}




#################      Category          Not in Category       Total
## Group1      length(Group1Cat)        length(OnlyGroup1)       m 
## Group2      length(Group2Cat)        length(OnlyGroup2)       n
# Total                  k
                 
### Enrichment test for GRanges. Test for enrichment of Category within Group1
### @Group1 - A GRanges object for one set of positions
### @Group2 - The background GRanges object, non-overlapping with Group1.
### @Category - A GRanges object of known locations, such as motifs, that you want to test for enrichment in Group1.
### @type - Default is null. You can use this to pull out or simplify the test to a metadata column within the GRanges 
###          for Group1 and Group2. For example, if you want to test for enrichment of all genes, instead of open regions. 
###          If type = null, then it will just use the number of Ranges instead of the number of unique 
###           entries in column 'type'
EnrichedRanges <- function(Group1, Group2, Category, type = NULL, returnTable = FALSE){
    
    Group1Cat <- filter_by_overlaps(Group1, Category) 
    Group2Cat <- filter_by_overlaps(Group2, Category)
   
    OnlyGroup1 <- filter_by_non_overlaps(Group1, Category)
    OnlyGroup2 <- filter_by_non_overlaps(Group2, Category)
    
    if(returnTable & is.null(type)){

        dt_table <- data.frame(Group1 = c(length(Group1Cat), length(OnlyGroup1)), 
                               Group2 = c(length(Group2Cat), length(OnlyGroup2)), 
                       row.names = c('In Category', 'Not in Category')) 
    
        return(t(dt_table))
        
    }else if(returnTable & 
             sum(c(colnames(mcols(Group1)),colnames(mcols(Group2))) %in% type) == 2 &
             length(type) == 1){
        
       dt_table <- data.frame(Group1 = c(length(unique(mcols(Group1Cat)[,type])), 
                                          length(unique(mcols(OnlyGroup1)[,type]))), 
                               Group2 = c(length(unique(mcols(Group2Cat)[,type])), 
                                         length(unique(mcols(OnlyGroup2)[,type]))), 
                       row.names = c('In Category', 'Not in Category')) 
    
        return(t(dt_table))
       
    }else if(returnTable){
        
        stop('Error: Incorrect method or column name. Please check input')      
    }
    
    if(is.null(type)){
        
       pVal <- phyper(q = length(Group1Cat), 
           m = length(Group1), 
           n = length(Group2),
           k = length(Group1Cat) + length(Group2Cat),
            lower.tail=FALSE)

       enrichment <- (length(Group1Cat)/length(Group1))/(length(Group2Cat)/length(Group2))
        
    }else if(sum(c(colnames(mcols(Group1)),colnames(mcols(Group2))) %in% type) == 2 &
             length(type) == 1){
        
       pVal <- phyper(q = length(unique(mcols(Group1Cat)[,type])), 
           m = length(unique(mcols(Group1)[,type])), 
           n = length(unique(mcols(Group2)[,type])),
           k = length(unique(mcols(Group1Cat)[,type])) + length(unique(mcols(Group2Cat)[,type])),
                     lower.tail=FALSE)
        
       enrichment <- (length(unique(mcols(Group1Cat)[,type]))/length(unique(mcols(Group1)[,type])))/
                    (length(unique(mcols(Group2Cat)[,type]))/length(unique(mcols(Group2)[,type])))
        
    }else{
        
        stop('Error: Incorrect method or column name. Please check input')
        
    }
    
    return(data.frame(p_value = pVal, enrichment = enrichment))
    
}
    
######## Test all motifs for enrichment.
    

MotifEnrichment <- function(Group1, Group2, motifPosList, type = NULL, numCores = 1, returnTable = FALSE){
    
    
    allEnrichmentList <- mclapply(motifPosList, function(x){
        
        tmp_df <- EnrichedRanges(Group1, Group2, Category = x, type = type)
        
    }, mc.cores = numCores)
    df_final <- do.call('rbind',  allEnrichmentList)

    df_final$adjp_val <- p.adjust(df_final$p_value)
    df_final$mlog10Padj <- -log10(adjp_val)
    
    return(df_final)

}
    
############# Pull out all the motifs associated with each gene according to a set of TSS sites

## @TSS_Sites - GRanges objects that are the list of TSS sites of interest. 
##              Must include a column 'name' which has the associated gene name
## @allPeaks - GRanges object of all peaks 
## @TSS_Links - a data.table object that record all the peak-peak links by co-accessibility
##              Must include columns named 'Peak1' and 'Peak2' which contain a string describing
##              each peak in the format 'chr1:100-2000' and must be identical to peaks listed in
##              allPeaks
## @motifPosList - a GRangesList, which each index is a GRanges of all positions 
##                  for a given motif. GRangesList must be named. 
## @numCores - number of cores to multithread over. 
    
    
Gene2Motif <- function(TSS_Sites, allPeaks, TSS_Links, motifPosList, 
                       numCores = 1, verbose = FALSE){
    
    if(verbose){ print('Generating TSS-Peak Network.')}
    
    TSS_Network <- c(TSS_Links$Peak1, TSS_Links$Peak2, 
                           GRangesToString(TSS_Sites)) %>%
                   unique() %>%
                  StringsToGRanges(.) %>% 
                plyranges::filter_by_overlaps(allPeaks, .)
    
    if(verbose){ print('Finding all motifs related to each peak within the TSS-Peak Network.')}
    ##Let's find all the motifs that overlap with each peak within the altTSS Network
    tmpOverlap <- mclapply(seq_along(motifPosList), function(x){
    
        ifelse(count_overlaps(TSS_Network, motifPosList[[x]]) > 0,
           names(motifPosList)[x], NA)
    
    }, mc.cores= numCores)
    
    
    overlap_df <- do.call('cbind', tmpOverlap)
    motifList <- mclapply(c(1:dim(overlap_df)[1]), function(x){
        
        
        ifelse(any(!is.na(overlap_df[x,])),
            list(overlap_df[x,which(!is.na(overlap_df[x,]))]),
            NA)
    
    }, mc.cores= numCores)

    if(verbose){ print('Finding all peaks related to each gene within the TSS-Peak Network.')}
    ##Find all the peaks related to each gene. 
    
    Peak2Gene <- mclapply(unique(TSS_Sites$name), function(x){
    
        geneTSS <- plyranges::filter(TSS_Sites, name == x)  %>% 
                plyranges::filter_by_overlaps(allPeaks, .) %>%
            plyranges::ungroup() %>%
            GRangesToString(.)
    
        tmp <- TSS_Links[Peak1 %in% geneTSS | Peak2 %in% geneTSS,]
    
        if(dim(tmp)[1] > 0){ 
            unique(c(tmp$Peak1, tmp$Peak2, geneTSS))
        }else{
               geneTSS
        }
    
    }, mc.cores = numCores)
    names(Peak2Gene) <- unique(TSS_Sites$name)
    
    

    if(verbose){ print('Linking Motifs to each gene within the TSS-Peak Network')}
    ## Link all the genes to motifs via Peak2Gene and the motifList
    Gene2Motif <- mclapply(Peak2Gene, function(x){
    
        #Find which indices of the AltTSS Network GRanges are linked to that Gene
        tmp <- findOverlaps(StringsToGRanges(x), TSS_Network)
        #Pull up and unlist all the motifs associated with those tiles.
        unlist(motifList[subjectHits(tmp)])
    
    }, mc.cores = numCores)
    
    return(Gene2Motif)
    
}
    
    