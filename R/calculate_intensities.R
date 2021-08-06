#' @title \code{calculate_intensities}
#'
#' @description \code{calculate_intensities} is an R helper function, part of the single-cell peak calling
#' algorithm scMACS by (XXX et al, 2021) that calculates the design matrix, X,
#' used to make peak calls.
#'
#' @param fragMat a genomic ranges object containing the fragments
#' @param PeakRegion a genomic ranges object indicating where to call peaks. This could be a pre-selected peak set or a dynamic bins approach to searching the whole genome.
#' @param normalizeBins a boolean indicating whether varying bin widths should be normalized to a fixed width
#'
#' @return a data.table, countsByBin, that returns the two intensity parameters required
#' to calculate the probability of a (+) peak
#'
#' @details The technical details of the algorithm are found in XX.
#'
#' @references XX
#'


calculate_intensities <- function(fragMat,
                                  positivePeaks,
                                  normalizeBins=TRUE) {
  print("calculating intensities from ArchR Project!")

  fragMat_dt <- as.data.table(fragMat)

  ## transform granges to data.table
  positivePeaksDF <- as.data.table(positivePeaks)
  positivePeaksDF$bin <- paste(positivePeaksDF$seqnames, ':',positivePeaksDF$start,
                               '-',positivePeaksDF$end,sep='')

  ###
  fragsPerBin <- findOverlaps(fragMat,
                              positivePeaks,
                              minoverlap=0)

  fragsPerBin <- as.data.table(fragsPerBin)

  numCells = length(unique(fragMat_dt$RG))

  # 		fragsPerBin <- fragsPerBin[order(fragsPerBin$subjectHits),]
  fragsPerBin$cell <- fragMat_dt$RG[fragsPerBin$queryHits]

  fragsPerBin$bin <-  positivePeaksDF$bin[fragsPerBin$subjectHits]
  fragsPerBin$chr <-  positivePeaksDF$seqnames[fragsPerBin$subjectHits]
  fragsPerBin$start <- positivePeaksDF$start[fragsPerBin$subjectHits]
  fragsPerBin$width <- positivePeaksDF$width[fragsPerBin$subjectHits]


  if(normalizeBins){
    ### get cell count matrix
    cell_counts = fragsPerBin[,list(N=.N,
                                    width=unique(width)),
                              by=list(bin,cell)
    ]
    cell_counts$normed_width <- pmax(1, round(cell_counts$width/500))
    cell_counts$normed_counts <- pmax(1,round(cell_counts$N/cell_counts$normed_width))

    #### doing bin-level summaries
    setkey(cell_counts, bin)
    cat('\n\nCalculating Region-Normalized reads per cell using 500bp wide bins\n\n')

    countsByBin <-  cell_counts[, list(CellsWithNoReads=(numCells-sum(N>=1))/numCells,
                                       ReadsPerCell=sum(N)/numCells,
                                       maxIntensity = max(N)
    ), by=bin
    ]

  } else{
    ### get cell count matrix
    cell_counts = fragsPerBin[,list(N=.N),
                              by=list(bin,cell)
    ]

    #### doing bin-level summaries
    setkey(cell_counts, bin)

    #### calculate features
    countsByBin <-  cell_counts[, list(CellsWithNoReads=(numCells-sum(N>=1))/numCells,
                                       ReadsPerCell=sum(N)/numCells,
                                       maxIntensity = max(N)
    ), by=bin
    ]
  }

  countsByBin <- dplyr::left_join(positivePeaksDF[, c('seqnames','start','bin','width')],
                                  countsByBin,
                                  by='bin')

  countsByBin <- countsByBin[order(countsByBin$seqnames, countsByBin$start, decreasing=FALSE),]

  ## creating counts for bins with no reads
  countsByBin$CellsWithNoReads[is.na(countsByBin$CellsWithNoReads)] <- 1
  countsByBin$ReadsPerCell[is.na(countsByBin$ReadsPerCell)] <- 0
  countsByBin$maxIntensity[is.na(countsByBin$maxIntensity)] <- 0

  countsByBin$threeBin <- rollapply(countsByBin$ReadsPerCell, FUN=mean, width=3, by=1 ,partial=TRUE, alight='left')
  countsByBin$numCells <- numCells


  chromosomeOrder <- c("chr1"  ,"chr2" , "chr3",  "chr4",
                       "chr5" , "chr6" , "chr7"  ,"chr8",
                       "chr9" ,"chr10" ,"chr11", "chr12",
                       "chr13", "chr14", "chr15", "chr16",
                       "chr17", "chr18" ,"chr19", "chr20",
                       "chr21", "chr22", "chrX", 'chrY' )

  countsByBin$chromosome <- factor(countsByBin$seqnames, levels=chromosomeOrder, ordered=T)
  countsByBin <- countsByBin[order(countsByBin$chromosome, countsByBin$start, decreasing=FALSE),]

  countsByBin$threeBin <- rollapply(countsByBin$ReadsPerCell, FUN=mean, width=3, by=1 ,partial=TRUE, alight='left')


  print(paste('Analysis finished on ', numCells, 'cells'))
  return(countsByBin)
}
