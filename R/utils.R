# Function to get sample-level metadata,
# from an ArchR project's colData
sampleDataFromCellColData <- function(cellColData, sampleLabel){
  
    if (!(sampleLabel %in% colnames(cellColData))){
      stop(paste(
          "`sampleLabel` must present in your ArchR Project's cellColData",
          "Check `names(getCellColData(ArchRProj)` for possible sample columns."
      ))
    }
  
    # Drop columns where all values are NA
    cellColDataNoNA <- BiocGenerics::Filter(function(x){!all(is.na(x))}, cellColData)

    # Convert to data.table
    cellColDT <- data.table::as.data.table(cellColDataNoNA)

    BoolDT <- cellColDT[, lapply(.SD, function(x){length(unique(x)) == 1}), by = c(sampleLabel)]
    trueCols <- apply(BoolDT, 2, all)
    trueCols[[sampleLabel]] <- TRUE
    cellColDF <- as.data.frame(cellColDT)                             

    sampleData <- dplyr::distinct(cellColDF[,names(which(trueCols))])

    # Set sampleIDs as rownames
    rownames(sampleData) <- sampleData[[sampleLabel]]
    return(sampleData)
}


# Function to split the output of getPopFrags into a list 
# of lists of GRanges, one named for each celltype.
# Resulting in a list containing each celltype, and each 
# celltype has a list of GRanges name for each sample.
splitFragsByCellPop <- function(frags) {
    
    # Rename frags by cell population
    renamedFrags <- lapply(
        1:length(frags),
        function(y){
            # Split out celltype and sample from the name
            x <- frags[y]
            celltype_sample <- names(x)
            splits <- unlist(stringr::str_split(celltype_sample, "#"))
            celltype <- splits[1]
            sample <- unlist(stringr::str_split(splits[2], "__"))[1]
            # Rename the fragments with just the sample
            names(x) <- sample
            # Return as a list named for celltype
            output <- list(x)
            names(output) <- celltype
            output
        }
    )
    
    # Group frags by cell population
    renamedFrags <- unlist(renamedFrags, recursive=FALSE)
    splitFrags <- split(renamedFrags, f=names(renamedFrags))
    return(splitFrags)
}


# This function converst scMACS DAPs data table to a Summarized Experiment 
# for use with ArchR's functions
# This functions expects a DAPs data.table, with a column named 'Peak'
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



# This function converst scMACS sample-specific peak calls, raw fragments by sample, 
# the sample-specific peak matrix, and sample metadata into one Summarized Experiment Object
# It expects a column in the peak matrix called 'tileID', and a column in the metadata called 'Sample'

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
    


# Tests if a string is a in the correct format to convert to GRanges 
validRegionString <- function(regionString) {
  if (!is.character(regionString)) {
    return(FALSE)
  }

  pattern <- "([0-9]{1,2}|chr[0-9]{1,2}|chr[X-Y]{1,1}):[0-9]*-[0-9]*"
  matchedPattern <- str_extract(regionString, pattern)

  if (any(is.na(matchedPattern))) {
    return(FALSE)
  } else if (any(!matchedPattern == regionString)) {
    return(FALSE)
  }

  splits <- str_split(regionString, "[:-]")[[1]]
  start <- splits[2]
  end <- splits[3]
  if (any(start > end)) {
    return(FALSE)
  }

  # All conditions satisfied
  return(TRUE)
}

#' @title \code{StringsToGRanges}
#'
#' @description \code{StringsToGRanges} Turns a list of strings in the format chr1:100-200 
#'   into GRanges
#'
#' @export
StringsToGRanges <- function(regionString) {
  # if (length(regionString)>1){
  #   boolList <- lapply(regionString, function(x){validRegionString(x)})
  #   if (!all(boolList)){
  #     stop("Some region strings are invalid. Given regions must all be strings matching format 'seqname:start-end', where start<end e.g. chr1:123000-123500")
  #   }
  # } else if(!validRegionString(regionString)) {
  #   stop("Region must be a string matching format 'seqname:start-end', where start<end e.g. chr1:123000-123500")
  # }

  chrom <- gsub(":.*", "", regionString)
  startSite <- gsub(".*:", "", regionString) %>%
    gsub("-.*", "", .) %>%
    as.numeric()
  endSite <- gsub(".*-", "", regionString) %>% as.numeric()

  if(any(startSite >= endSite)){stop('Error in region string: Make sure the start of the genomic range occurs before the end') }
  regionGRanges <- GRanges(seqnames = chrom, ranges = IRanges(start = startSite, end = endSite), strand = "*")
  return(regionGRanges)
}        

#' @title \code{GRangesToString} Converts a GRanges object to a string in the format 'chr1:100-200'
#'
#' @description \code{GRangesToStrin} Turns a GRanges Object into
#' 	a list of strings in the format chr1:100-200
#'
#'
#' @export            
GRangesToString <- function(GR_obj){
  
  paste(seqnames(GR_obj), ":",start(GR_obj),"-",end(GR_obj),sep="")
  
}

#' @title \code{differentialsToGRanges} Converts a data.frame matrix to a GRanges,
#'   preserving additional columns as GRanges metadata
#'
#' @param differentials a matrix/data.frame with a column tileColumn containing 
#'   region strings in the format "chr:start-end"
#' @param tileColumn name of column containing region strings. Default is "Tile".
#'
#' @return a GRanges containing all original information 
#' @export 
differentialsToGRanges <- function(differentials, tileColumn = "Tile"){
  regions <- scMACS::StringsToGRanges(differentials[[tileColumn]])
  GenomicRanges::mcols(regions) <- differentials
  regions
}

### Function to set seqinfo for GRanges objects.
### This makes them easier to merge and export to BWs by setting the standard
### chromosome lengths and names to a standard, as defined by TxDb.

setSeqInfo <- function(GRangesObj, TxDb = TxDb.Hsapiens.UCSC.hg38.refGene) {
  seqinfo(GRangesObj) <- seqinfo(TxDb)[seqnames(seqinfo(GRangesObj))]
  return(GRangesObj)
}


# This functions takes the output of coaccessibility, filters by a 
# correlation threshold, and then returns a data.frame
#
# tielCorrelations = output of coaccessibility
# nearbyTiles - data.table of peaks within the correlation
# threshold - Absolute correlation threshold
createLinks <- function(nearByTiles, tileCorrelations, threshold = 0.5){
  
  filCorr <- tileCorrelations %>% filter(Peak1 != Peak2) %>%
    filter(abs(Correlation) > threshold)
  region <- makeGRangesFromDataFrame(nearByTiles) %>% 
    plyranges::mutate(idx = c(1:length(.))) %>%
    plyranges::mutate(mid = as.integer((start+end)/2)) 
  filCorr$midPeak1 <- region$mid[match(filCorr$Peak1,region$idx)]
  filCorr$midPeak2 <- region$mid[match(filCorr$Peak2,region$idx)]
  filCorr[!duplicated(t(apply(filCorr[2:3], 1, sort))), ]
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
    
  }else if(sum(c(colnames(mcols(Group1)),colnames(mcols(Group2))) %in% type) == 2 &
           length(type) == 1){
    
    pVal <- phyper(q = length(unique(mcols(Group1Cat)[,type])), 
                   m = length(unique(mcols(Group1)[,type])), 
                   n = length(unique(mcols(Group2)[,type])),
                   k = length(unique(mcols(Group1Cat)[,type])) + length(unique(mcols(Group2Cat)[,type])),
                   lower.tail=FALSE)
    
  }else{
    
    stop('Error: Incorrect method or column name. Please check input')
    
  }
  
  return(pVal)
  
}


