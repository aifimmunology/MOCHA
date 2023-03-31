#' @title \code{bulkLSI}
#'
#' @description \code{bulkLSI} generates LSI (Latent semantic indexing)
#'
#' @param SampleTileObj The SummarizedExperiment object output from getSampleTileMatrix
#' @param cellPopulations vector of strings. Cell subsets for which to call
#'   peaks. This list of group names must be identical to names that appear in
#'   the SampleTileObj.  Optional, if cellPopulations='ALL', then peak
#'   calling is done on all cell populations. Default is 'ALL'.
#' @param componentNumber integer. Number of components to include in LSI. This must be strictly less than
#' the number of samples times the number of cellTypes in your SampleTileObj.
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#'
#' @return LSIObj a SummarizedExperiment containing PC components from the LSI and metadata from the SampleTileObj
#' 
#' @examples
#' \dontrun{
#' LSIObj <- MOCHA::bulkLSI(SampleTileObj, cellType = 'CD16_Mono')
#' }
#' @export
#'
#' 
bulkLSI <- function(SampleTileObj, cellType = 'All', componentNumber = 30, verbose = FALSE){
    
    allCellTypes = names(SummarizedExperiment::assays(SampleTileObj))
    if(all(tolower(cellType) == 'all')){

        fullObj <- combineSampleTileMatrix(SampleTileObj)
        countMat <- SummarizedExperiment::assays(fullObj)[[1]]

    }else if(all(cellType %in% allCellTypes)){
        newTSAM <- subsetMOCHAObject(SampleTileObj, subsetBy = 'celltype',
                                    groupList = cellType, subsetPeaks = TRUE,
                                    verbose = verbose)
        fullObj <- combineSampleTileMatrix(newTSAM)
        countMat <- SummarizedExperiment::assays(fullObj)[[1]]

    }else{
        stop('cellType not found. SampleTileObj must contain the given cellType.')
    }

    #TF-IDF step
    freqs <- t(t(countMat)/Matrix::colSums(countMat))
    idf   <- log(1 + ncol(countMat) / Matrix::rowSums(countMat))
    tfidf <- Matrix::Diagonal(x=as.vector(idf)) %*% freqs

    # SVD step, using 30 components
    # max(nu, nv) must be strictly less than min(nrow(A), ncol(A))
    tryCatch(
      {svd <- irlba::irlba(tfidf, componentNumber, componentNumber)},
      error=function(cond){
        message(cond)
        error("Columns containing all NAs may be present in SampleTileObj")
      }
    )
    svdDiag <- matrix(0, nrow=componentNumber, ncol=componentNumber)
    diag(svdDiag) <- svd$d
    matSVD <- t(svdDiag %*% t(svd$v))
    rownames(matSVD) <- colnames(countMat)
    colnames(matSVD) <- paste0("PC",seq_len(ncol(matSVD)))   

    assayList1 <- list(t(matSVD))
    names(assayList1) = 'LSI'
    LSIObj <- SummarizedExperiment::SummarizedExperiment(
            assayList1,
            metadata = fullObj@metadata,
            colData = SummarizedExperiment::colData(fullObj)
    )
    return(LSIObj)
}

#' @title \code{bulkUMAP}
#'
#' @description \code{bulkUMAP} generates UMAP from pseudobulk LSIObj object, and merges in metadata.
#'
#' @param LSIObj The SummarizedExperiment object output from bulkLSI. 
#' 
#' @param components A vector of integers. Number of components to include in LSI (1:30 typically).
#' @param n_neighbors See  \link[uwot]{umap}. The size of local neighborhood (in terms of number of
#'           neighboring sample points) used for manifold approximation. Default is 15.
#' 
#' @return fullUMAP data.frame of UMAP values with metadata attached. 
#' 
#' @examples
#' \dontrun{
#' UMAPvalues <- MOCHA::bulkUMAP(LSIObj)
#' }
#' @export
#'
bulkUMAP <- function(LSIObj, components = c(1:30), n_neighbors = 15){

    countMat <- t(SummarizedExperiment::assays(LSIObj)[[1]])

    if(!all(components %in% seq_along(colnames(countMat)))){
      stop('Component list does not align with number of LSI components.')
    }
    
    subUMAP <- as.data.frame(
      uwot::umap(countMat[,components], n_neighbors=n_neighbors)
    )

    colnames(subUMAP) <- c('UMAP1', 'UMAP2')
    subUMAP$Sample <- rownames(subUMAP)

    fullUMAP <- dplyr::full_join(
      subUMAP, 
      as.data.frame(SummarizedExperiment::colData(LSIObj)),
      by = 'Sample',
      copy=TRUE
    ) 
    
    return(fullUMAP)
}
