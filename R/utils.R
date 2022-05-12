




### This function converst scMACS DAPs data table to a Summarized Experiment 
## for use with ArchR's functions
## This functions expects a DAPs data.table, with a column named 'Peak'
DAPsToSE <- function(daps){
    
    #Generate a GRanges from the character strings found in the peaks column
    tmp1 <- StringsToGRanges(daps$Peak)
    #Generate a data table with the rest of the info.
    tmp2 <- daps[,!'Peak']
    mcols(tmp1) <- tmp2
    tmp3 <- sort(tmp1)
    #Generate the assay list to put into the summarized experiment
    assaysList = lapply(colnames(tmp2), function(x){
                    as.matrix(mcols(tmp3)[,x])
              })
    #Make sure it's named appropriately.
    names(assaysList) = colnames(tmp2)

    rowData = as.data.frame(granges(tmp3))
    rowData$idx = seq(1:length(tmp3))
    #Now load into a summarized experiment format. 
    SummarizedExperiment(
         assays = assaysList,
          #rowRanges = granges(tmp3),
          rowData = rowData,
          metadata = list(Params = list(useMatrix = "PeakMatrix"))
    )
}



### This function converst scMACS sample-specific peak calls, raw fragments by sample, 
## the sample-specific peak matrix, and sample metadata into one Summarized Experiment Object
## It expects a column in the peak matrix called 'tileID', and a column in the metadata called 'Sample'

PeakMatToSE <- function(sample_peak_matrix, metadata, sampleSpecific){
    
    peaks <- StringsToGRanges(sample_peak_matrix$tileID)
    lambdas <- sample_peak_matrix[,!'tileID']
    
    colData1 <- as.data.table(metadata)
    rownames(colData1) <- metadata$Sample
    colData1 <- colData1[,!'Sample']

    SummarizedExperiment(
         assays = list(ReproduciblePeaks = as.matrix(lambdas)),
         rowRanges = peaks,
         colData = colData1,
         metadata = sampleSpecific
    )
}
    


### Quick support functions
##Strings to GRanges

StringsToGRanges <- function(regionString){
  
  chrom <- gsub(":.*","",regionString)
  startSite <- gsub(".*:","",regionString) %>% gsub("-.*", "",.) %>% as.numeric()
  endSite <- gsub(".*-","",regionString) %>% as.numeric()
  
  regionGRanges <- GRanges(seqnames = chrom, ranges = IRanges(start = startSite,end = endSite), strand = "*")
  
  return(regionGRanges)
  
}


