#' @title \code{pseudobulkUMAP}
#'
#' @description \code{pseudobulkUMAP} runs UMAP dimensionality reduction on a SampleTileMatrix, and returns a 
#'  SummarizedExperiment 
#'
#' @param Object A MultiAssayExperiment or RangedSummarizedExperiment,
#' @param subsetBy the variable to subset by. Can either be 'celltype', or a
#'   column from the sample metadata (see colData(Object))
#' @param groupList the list of cell type names or sample-associated data that
#'   should be used to subset the Object
#' @param na.rm removes groups that are NA if set to true. If set to false, then
#'   you filter for everything in the groupList and also NA values.
#' @param subsetPeaks If subsetting by cell types, then you need to decide if you want to subset the tile set down to tiles 
#'   only called in those cell types. The default is TRUE. 
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#'
#' @return Object the input Object, filtered down to either the cell type or
#'   samples desired.
#'
#'
#' @export

pseudobulkUMAP <- function(STObject,
                              cellTypeSpecific = FALSE,
                            columnName = '',
                            numCores = 1) {

     cellTypes <- names(SummarizedExperiment::assays(STObject))

     if(cellTypeSpecific){

       
        cl <- parallel::makeCluster(numCores)
        
        suppressMessages(allUMAPs <- pbapply::pblapply(cellTypes, function(x){

                accMat <- getCellPopMatrix(STObject,  cellPopulation = x)
                uwot::umap(t(accMat))

        }, cl = cl), classes = "message")

        parallel::stopCluster(cl)

        names(allUMAPs) <- cellTypes

        meta1 <- SummarizedExperiment::colData(STObject)

     }else{


        allMat <- do.call('cbind', as.list(SummarizedExperiment::assays(STObject)))

        colnames(allMatrices) <- apply(expand.grid(colnames(STObject),cellTypes), 1, 
                                paste, collapse="__") %>% gsub(" ", "", .)

        allUMAPs <- list(uwot::umap(t(allMatrices)))

        names(allUMAPs) <- 'AllCellTypes'

        meta1 <- lapply(cellTypes, function(x) {
            
                    tmpMeta <- colData(STObject)
                    tmpMeta$Sample <- paste(tmpMeta$Sample, x, sep = "__") %>% gsub(" ", "", .)
                    tmpMeta$CellTypeAnnot <- rep(x, dim(tmpMeta)[1])
                    rownames(tmpMeta) <- tmpMeta$Sample

                }) %>% do.call('rbind',.)

     }

    umapObject <- SummarizedExperiment(assays = allUMAPs, colData = meta1)

    return(umapObject)

}