#' @title \code{bulkLSI}
#'
#' @description \code{bulkLSI} generates LSI 
#'
#' @param TSAM_Object The SummarizedExperiment object output from getSampleTileMatrix
#' @param cellPopulations vector of strings. Cell subsets for which to call
#'   peaks. This list of group names must be identical to names that appear in
#'   the SampleTileObj.  Optional, if cellPopulations='ALL', then peak
#'   calling is done on all cell populations. Default is 'ALL'.
#' @param componentNumber integer. Number of components to include in LSI. 
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#'
#' @return LSI_SE a SummarizedExperiment containing PC components from the LSI and metadata from the TSAM_Object
#' 
#' @examples
#' \dontrun{
#' LSIse <- MOCHA::bulkLSI(TSAM, cellType = 'CD16 Mono')
#' 
#' @export
#'
#' 
bulkLSI <- function(TSAM_Object, cellType = 'All', componentNumber = 30, verbose = FALSE){
    ## code adapted from https://github.com/GreenleafLab/10x-scATAC-2019/blob/master/code/02_Get_Peak_Set_hg19_v2.R
    
    allCellTypes = names(SummarizedExperiment::assays(TSAM_Object))
    if(all(tolower(cellType) == 'all')){

        fullObj <- combineSampleTileMatrix(TSAM_Object)
        countMat <- SummarizedExperiment::assays(fullObj)[[1]]

    }else if(all(cellType %in% allCellTypes)){
        newTSAM <- subsetMOCHAObject(TSAM_Object, subsetBy = 'celltype',
                                    groupList = cellType, subsetPeaks = TRUE,
                                    verbose = verbose)
        fullObj <- combineSampleTileMatrix(newTSAM)
        countMat <- SummarizedExperiment::assays(fullObj)[[1]]

    }else{
        stop('Cell type names not found.')
    }

    #TF-IDF step
    freqs <- t(t(countMat)/Matrix::colSums(countMat))
    idf   <- log(1 + ncol(countMat) / Matrix::rowSums(countMat))
    tfidf <- Matrix::Diagonal(x=as.vector(idf)) %*% freqs

    #SVD step, using 30 components
    svd <- irlba::irlba(tfidf,componentNumber, componentNumber)
    svdDiag <- matrix(0, nrow=componentNumber, ncol=componentNumber)
    diag(svdDiag) <- svd$d
    matSVD <- t(svdDiag %*% t(svd$v))
    rownames(matSVD) <- colnames(countMat)
    colnames(matSVD) <- paste0("PC",seq_len(ncol(matSVD)))   

    assayList1 <- list(t(matSVD))
    names(assayList1) = 'LSI'
    newExp <- SummarizedExperiment(
            assayList1,
            metadata = fullObj@metadata
            colData = SummarizedExperiment::colData(fullObj)
    )
    return(newExp)
}

#' @title \code{bulkUMAP}
#'
#' @description \code{bulkUMAP} generates UMAP from pseudobulk LSI_SE object, and merges in metadata.
#'
#' @param LSI_SE The SummarizedExperiment object output from bulkLSI. 
#' 
#' @param componentNumber vector of integers. Number of components to include in LSI (1:30 typically)
#' 
#' @return data.frame of UMAP values with metadata attached. 
#' 
#' @examples
#' \dontrun{
#' LSIse <- MOCHA::bulkLSI(TSAM, cellType = 'CD16 Mono')
#' 
#' @export
#'
#' 

bulkUMAP <- function(LSI_SE, components = c(1:30)){

    countMat <- SummarizedExperiment::assays(LSI_SE)[[1]]
    metaData <- as.data.frame(SummarizedExperiment::colData(LSI_SE))

    if(all(components %in% seq_along(colnames(countMat)))){
        stop('Component list does not align with number of LSI components.')
    }
    
    subUMAP <- as.data.frame(uwot::umap(t(countMat)[,components]))

    colnames(subUMAP) <- c('UMAP1', 'UMAP2')
    subUMAP$Sample <- rownames(subUMAP)

    fullUMAP <- dplyr::full_join(subUMAP, as.data.frame(metadata), 
                    by = 'Sample') 
    
    return(fullUMAP)
}
