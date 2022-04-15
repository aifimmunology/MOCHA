#' @title \code{format_peaks.R}
#'
#' @description \code{format_peaks.R} allows you to format data.tables into
#'              granges by creating the seqnames, start, end fields from a string
#' 
#'
#' @param peaks a peak matrix with a string containing genomic address,
#'         'chr:XXX-XXX', where chromosome is delimited with colon,
#'                        and start/end sites with a dash 
#'
#'
#' @return a granges object 
#'
#'
#' @export

format_peaks <- function(peaks){

    fields <- strsplit(peaks$Tile, '\\:|\\-')
    ## extract chromosome
    peaks$seqnames <- sapply(fields, function(x) unlist(x)[1]
                                 )
    ## extract start site    
    peaks$start <- sapply(fields, function(x)  unlist(x)[2]
                                 )
    ## extract end site    
    peaks$end <- sapply(fields, function(x) unlist(x)[3]
                                 )
    peaks_gr <- GenomicRanges::makeGRangesFromDataFrame(peaks,
                                                       keep.extra.columns=TRUE)
    return(peaks_gr)
    
}