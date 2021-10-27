#' @title \code{create_peak_sampleMatrix}
#'
#' @description \code{create_peak_sampleMatrix} is a function that can 
#'              transform a set of sample-specific peak calls into 
#'              peak X sample matrix containing lambda1 measurements for 
#'              each sample.
#'
#'
#' @param sample_specific_peaks output of calls_peaks_by_sample (peaks only)
#'
#' @return sample_peak_matrix a sample X peak matrix containing observed 
#'         measurements for each sample at each peak. 
#' 
#' @details The technical details of the algorithm are found in XX.
#'
#' @references XX
#'
#' @export

create_peak_sampleMatrix <- function(sample_specific_peaks){
    
    
    tmp_peaks <- sample_specific_peaks[['scMACS_peaks']]

    for(i in 1:length(tmp_peaks)){

        tmp_peaks[[i]]$Sample <- names(tmp_peaks)[i]

    }

    fullMatrix <- do.call(rbind, tmp_peaks[!names(tmp_peaks) %in% 'Union'])
    fullMatrix <-  fullMatrix[,c('lambda1','PeakID','Sample')]


    sample_peak_matrix <- dcast(fullMatrix, formula=PeakID ~ Sample, value.var='lambda1')

    sample_peak_matrix[is.na(sample_peak_matrix)] <- 0
    return(sample_peak_matrix)
    
}