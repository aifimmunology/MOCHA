####################################################################################################
###### TF Motif Set Enrichment Analysis (MSEA) 
####################################################################################################

## This analaguous to Gene Set Enrichment Analysis. Instead of testing for enrichment of a geneset with a given gene set in a pathway,
## we are testing the enrichment of a given TF Motif set against a motif set downstream of a given ligand. 
## If there is enrichment, it's a sign that that ligand could drive that set of motifs. 

######## Function for testing TFs interactin with a ligand.
#a. NicheNet - count all TFs regulated by any ligand - n-all
#b. NicheNet, count all TFs regulated by the ligand of interest - n-1
#c. scATAC - count all enriched TFs of interest - m-all
#d. NicheNet-scATAC - count all overlab between B & C - m-1
#E-- phyper(m1, n1, na-n1, ma, lower.tail=F, log.p=F)

# @ ligand_tf_matrix - NicheNet Ligand-tf matrix
# @ MotifEnrichment - Motif enrichment dataframe, unfiltered, from ArchR's peakAnnoEnrich step.
# - make sure that motif enrichment was only run on either up or downregulated peaks, and that you do not filter 
#    after running motif enrichment.
# @ motifColumn is the column within the MotifEnrichmentUnfiltered dataframe that has motif names. Default is 'feature'
# @ specLigand - the name of the specific ligand you want to test
# @ stat_threshold : Significance threshold used to select significant motif set
# @ verbose: Prints the n-all, n-1, m-all, and m-1

PHyperLigandTF <- function(ligand_tf_matrix, MotifEnrichment,
                            specLigand, 
                            motifColumn = "feature", 
                            stat_column = 'log10Padj',
                            stat_threshold = 2,
                           verbose = FALSE){
                           
    MotifEnrichment <- as.data.frame(MotifEnrichment)
    
    allMotifNames <- MotifEnrichment[,motifColumn]
    
    #Test if you need to cleanup the motif names
    if(any(grepl("_",allMotifNames))){
        
       MotifEnrichment[,motifColumn] = gsub("_.*", "", allMotifNames)
    }
 
    allMotifNames <- unlist(MotifEnrichment[,motifColumn])
    
    #Check inputs, to make sure errors won't be generated. 
    if(!any(rownames(ligand_tf_matrix) %in% allMotifNames)){
    
        stop("Error: Transcription Factor names don't match.")
    }
    
    if(!specLigand %in% colnames(ligand_tf_matrix)){
        
        stop("Error: Specific Ligand does not appear in NicheNet matrix.")
        
    }
    
    if(any(colnames(ligand_tf_matrix) %in% colnames(MotifEnrichment))){
        #Motif Enrichment columns and ligand_tf_matrix are named the same. This should not be. 
	#Let's rename the column to TranscriptionFactor
	colnames(MotifEnrichment)[colnames(MotifEnrichment) == motifColumn] = 'TranscriptionFactor'
	motifColumn = "TranscriptionFactor"
    }
    
    #filter ligand_tf_matrix down to those interactions with TFs within MotifEnrichment
    allTFsByAnyLigand =  ligand_tf_matrix[rownames(ligand_tf_matrix) %in% allMotifNames, ]
    allTFsByAnyLigand =  allTFsByAnyLigand[,colSums(allTFsByAnyLigand) > 0]
    
    if(! specLigand %in% colnames( allTFsByAnyLigand)){
     
        print(paste("Warning: the ligand named ", specLigand," does not have any interactions with motif set", sep =""))
        return(NA)
    }
    
    TFMat <- as.data.frame(allTFsByAnyLigand)
    TFMat$TranscriptionFactor <- rownames(allTFsByAnyLigand)
    
    #Check for number of overlaping Transcription Factors between NicheNet and MotifEnrichment set.
    NicheNet_Motif_overlap <- sum(TFMat$TranscriptionFactor %in% allMotifNames)
       
    #Merge data into one. 
    otherMotifs <- "TranscriptionFactor" 
    mergedDF <- inner_join(MotifEnrichment, TFMat, by = structure(names = motifColumn, .Data = otherMotifs)) %>% 
                    pivot_longer(cols = colnames(.)[c((dim(.)[2] - dim(allTFsByAnyLigand)[2] +1):dim(.)[2])], 
                                 names_to = 'Ligand', values_to = 'Score')
    
    nall <- dim(allTFsByAnyLigand)[1]
    n1 <- sum(allTFsByAnyLigand[,specLigand] > 0)

    #Filter to just significantly enriched motif 
    mergedDF_f <-  mergedDF[mergedDF[,stat_column] > stat_threshold,]
    
    #Number of significantly enrichment motifs
    mall <- length(unique(mergedDF_f$TranscriptionFactor)) 
    
    #Number of enriched TFs that could be downstream of a given ligand
    m1 <- sum(unique((mergedDF_f$TranscriptionFactor)) %in% rownames(allTFsByAnyLigand)[allTFsByAnyLigand[,specLigand] > 0])
    
    if(verbose){
    
            print(paste("TF Overlap: ",NicheNet_Motif_overlap, sep = ""))
            print(paste("All TF by Any Ligand: ",nall, sep = ""))
            print(paste("Significantly Enriched Motifs: ",mall, sep = ""))
            print(paste("Enriched TF Downstream of Ligand: ",m1, sep = ""))
    }
    phyper(m1, n1, nall-n1, mall, lower.tail=F, log.p=F)   
    
}
 


############################
### Wrapper for iterating over multiple ligands, and testing each one's motif set enrichment.   
# @ allMotifNames - list of names for all tested TF motifs.
# @ ligand_tf_matrix - NicheNet Ligand-tf matrix
# @ MotifEnrichment - Motif enrichment dataframe, unfiltered, from ArchR's peakAnnoEnrich step.
# - make sure that motif enrichment was only run on either up or downregulated peaks, and that you do not filter 
#    after running motif enrichment.
# @ - motifColumn is the column within the MotifEnrichmentUnfiltered dataframe that has motif names. Default is 'feature'
# @ ligands - vector of ligands to test

# @ stat_threshold : Significance threshold used to select significant motif set
# @ annotationName : If you want to annotate the output motifs with labels, you can provide the annotationName, which will be the column name for the annotation. 
# @ annotation : This is the annotation value that will be added to all rows of the output dataframe. Useful if you wrap this function in a loop or lapply statement over multiple cell types
# @ numcores: Number of cores to parallelize this over. 
# @ verbose: Prints the n-all, n-1, m-all, and m-1

MSEA <- function(ligand_tf_matrix, MotifEnrichment, 
                                       motifColumn = "feature", 
				                        ligands,
                                        stat_column = 'mlog10Padj',
                                       stat_threshold = 2, 
                                       annotationName = 'CellType', annotation = "none", 
                                       numCores = 1, verbose = FALSE){
    
    if(any(!ligands %in% colnames(ligand_tf_matrix)) | any(duplicated(ligands))){
        
        stop("Error: Some ligands does not appear in NicheNet matrix, or are duplicated.")
        
    }

     specificLigands <- mclapply(ligands, function(x){
    
            
        PHyperLigandTF(ligand_tf_matrix,MotifEnrichment, x, motifColumn = motifColumn, 
                                stat_column = stat_column,
                                stat_threshold = stat_threshold, verbose = verbose) 
                
    
            }, mc.cores = numCores)
    names(specificLigands) <- ligands
    #create data frame of ligand, p-val, adjusted p-value, and percentage of ligand-TFs in motif set.
    specDF <- data.frame(ligand = names(specificLigands), p_val = unlist(specificLigands))
    
    if(any(specDF$p_val == 0, na.rm =TRUE)){
        specDF$p_val[specDF$p_val == 0] = rep(1e-323, sum(specDF$p_val == 0))
    }
    specDF$adjp_val <- p.adjust(specDF$p_val)
    
    #Subset ligand matrix down to all TFs related to ligands
    subsetMat <- ligand_tf_matrix[rownames(ligand_tf_matrix) %in% MotifEnrichment[,motifColumn], 
                                  colnames(ligand_tf_matrix) %in% ligands]
	
    #Identify Significant motifs
    sigMotifs <- unlist(unique(MotifEnrichment[MotifEnrichment$mlog10Padj > stat_threshold,motifColumn]))
	
    #Subset ligand matrix down to Significant motifs that are associated with ligands
    sigSubsetMat <- ligand_tf_matrix[rownames(ligand_tf_matrix) %in% sigMotifs, 
                             colnames(ligand_tf_matrix) %in% ligands]
	
    #Calculate percentage of significant motifs against background of all motifs for each tested ligand
    specDF$PercentSigTF <- colSums(sigSubsetMat > 0)/colSums(subsetMat > 0)
    
    #Calculate percentage of significant motifs that interact with a given ligand against the background of all significant motifs
    specDF$PercInNicheNet <- colSums(sigSubsetMat > 0)/
                    length(sigMotifs)
    #Label these values
    specDF[, annotationName] <-  rep(annotation, nrow(specDF))
    specDF
    
}


