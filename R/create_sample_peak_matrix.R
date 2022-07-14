#' @title \code{create_peak_sampleMatrix}
#'
#' @description \code{create_peak_sampleMatrix} is a function that can 
#'              transform a set of sample-specific peak calls into 
#'              peak X sample matrix containing lambda1 measurements for 
#'              each sample.
#'
#'
#' @param peakCallResults output of calls_peaks_by_sample
#' @param reproduciblePeaks a vector containing the tileIDs to subset the peak matrix
#' @return samplePeakIntensityMat a sample X peak matrix containing observed 
#'         measurements for each sample at each peak. 
#' 
#' @details The technical details of the algorithm are found in XX.
#'
#' @references XX
#'
#' @export

create_peak_sampleMatrix <- function(
    peakCallResults,
    cellPopulation,
    reproduciblePeaks
){
    
    # Get the RaggedExperiment for this cellPopulation
    peaksExperiment <- peakCallResults[[cellPopulation]]
    
    # Extract matrices of samples by peak tileIDs with TotalIntensity
    samplePeakIntensityMat <- RaggedExperiment::compactAssay(
        peaksExperiment, i="TotalIntensity")
    
    # Replace NAs with zeroes for zero-inflated hypothesis testing 
    samplePeakIntensityMat[is.na(samplePeakIntensityMat)] <- 0
    
    # Filter to just peaks in the given reproduciblePeaks
    samplePeakIntensityMat <- samplePeakIntensityMat[reproduciblePeaks, ]
    
    return(samplePeakIntensityMat)
        
#     ## Extract Peak list from Sample Specific Peaks Object 
#     sample_names <-  sapply(sample_specific_peaks, function(x) names(x))
    
#     ## assign sample names to intensity matrix
#     for(i in 1:length(sample_specific_peaks)){

#         sample_specific_peaks[[i]][[1]]$Sample <- sample_names[i]

#     }
    
#     ## extract matrices
#     sample_specific_peaks <- lapply(sample_specific_peaks,
#                         function(x) x[[1]]
#                         )
    
#     ## concatenate matrix 
#     fullMatrix <- rbindlist(sample_specific_peaks)
    
#     ## Filter and rename feature names
#     fullMatrix <-  fullMatrix[,c('TotalIntensity','tileID','Sample')]
#     colnames(fullMatrix)[1] <- 'lambda1'

#     ## Reshape long to wide to create 
#     ## Sample X Peak matrix. 
#     sample_peak_matrix <- data.table::dcast(fullMatrix, 
#                                             formula=tileID ~ Sample, 
#                                             value.var='lambda1')
    
#     ## Replace NAs with zeroes for zero-inflated
#     ## hypothesis testing 
#     sample_peak_matrix[is.na(sample_peak_matrix)] <- 0 
    
#     if(is.null(union_peaks)){
        
#         return(sample_peak_matrix)
        
#     } else { 
        
#         return(sample_peak_matrix[tileID %in% union_peaks])
#     }
        
}