#' @title \code{calculate_intensities}
#'
#' @description \code{calculate_intensities} is an R helper function, part of the single-cell peak calling
#' algorithm scMACS by (XXX et al, 2021) that calculates the design matrix, X,
#' used to make peak calls.
#'
#' @param fragMat a genomic ranges object containing the fragments
#' @param candidatePeaks a genomic ranges object indicating where to call peaks. This could be a pre-selected peak set or a dynamic bins approach to searching the whole genome.
#' @param normalizeBins a boolean indicating whether varying bin widths should be normalized to a fixed width
#'
#' @return a data.table, countsByBin, that returns the two intensity parameters required
#' to calculate the probability of a (+) peak
#'
#' @details The technical details of the algorithm are found in XX.
#'
#' @usage calculate_intensities(fragMat, candidatePeaks, FALSE)
#'
#' @references XX
#' 
#' @export

calculate_intensities <- function(fragMat,
                                  candidatePeaks,
                                  NBDistribution,
                                  normalizeBins=FALSE) {
  print("calculating intensities from ArchR Project!")

  if(class(fragMat) != 'GRanges'){
    stop('fragMat user-input must be a Genomic Ranges (GRanges) object')
  }

  if(class(candidatePeaks) != 'GRanges'){
    stop('candidatePeaks user-input must be a Genomic Ranges (GRanges) object')
  }

  if(class(normalizeBins) != 'logical'){
    stop('normalizeBins user-input must be a TRUE/FALSE boolean indicating whether
         genomic regions should be normalized to fixed-width intervals')
  }

  fragMat_dt <- data.table::as.data.table(fragMat)

  ## transform granges to data.table
  candidatePeaksDF <- data.table::as.data.table(candidatePeaks)
  candidatePeaksDF$bin <- paste(candidatePeaksDF$seqnames, ':',candidatePeaksDF$start,
                               '-',candidatePeaksDF$end,sep='')

  ###
  fragsPerBin <- GenomicRanges::findOverlaps(fragMat,
                              candidatePeaks,
                              minoverlap=0)

  fragsPerBin <- data.table::as.data.table(fragsPerBin)

  numCells = length(unique(fragMat_dt$RG))

  fragsPerBin$cell <- fragMat_dt$RG[fragsPerBin$queryHits]
  fragsPerBin$bin <-  candidatePeaksDF$bin[fragsPerBin$subjectHits]
  fragsPerBin$chr <-  candidatePeaksDF$seqnames[fragsPerBin$subjectHits]
  fragsPerBin$start <- candidatePeaksDF$start[fragsPerBin$subjectHits]
  fragsPerBin$end <- candidatePeaksDF$end[fragsPerBin$subjectHits]
  fragsPerBin$width <- candidatePeaksDF$width[fragsPerBin$subjectHits]

  fragsPerBin = as.data.table(fragsPerBin)
  
  ### get cell count matrix
  cell_counts = fragsPerBin[,list(N=.N),
                            by=list(bin,cell)]
  #### doing bin-level summaries
  setkey(cell_counts, bin)

  #### calculate features
  countsByBin <-  cell_counts[, list(lambda1=sum(N)/numCells,
                                     maxIntensity = max(N)), by=bin
                             ]
  
  ### join to the original dynamic bins 
  ### to retain original variables

  countsByBin <- dplyr::left_join(candidatePeaksDF[, c('seqnames','start','end','bin')],
                                  countsByBin,
                                  by='bin')

  ### if the dynamic bins is calculated
  ### on a cohort, and applied to samples
  ### some entries will be NAs
  ### and this fixes with 0s
  
  countsByBin[is.na(countsByBin)] <- 0
  
  ### order by chromosome
  countsByBin <- countsByBin[order(countsByBin$seqnames, countsByBin$start, decreasing=FALSE),]
  countsByBin$numCells <- numCells

  chromosomeOrder <- c("chr1"  ,"chr2" , "chr3",  "chr4",
                       "chr5" , "chr6" , "chr7"  ,"chr8",
                       "chr9" ,"chr10" ,"chr11", "chr12",
                       "chr13", "chr14", "chr15", "chr16",
                       "chr17", "chr18" ,"chr19", "chr20",
                       "chr21", "chr22", "chrX", 'chrY')


  countsByBin$seqnames <- factor(countsByBin$seqnames, levels=chromosomeOrder, ordered=T)
  countsByBin <- countsByBin[order(countsByBin$seqnames, countsByBin$start, decreasing=FALSE),]
  
  NBDistribution <- calculateNBDistribution(countsByBin, theta=0.01)
    
  countsByBin$lambda2<-NBDistribution[countsByBin$maxIntensity+1]
  
  ### retain only features of interest and in order required
  countsByBin = countsByBin[,c('seqnames','start','end','lambda1','lambda2','numCells'),with=F]
  print(paste('Analysis finished on ', numCells, 'cells'))
    

  return(countsByBin)
}
