#' @title \code{linearModeling}
#'
#' @description \code{linearModeling} Runs linearModeling on MOCHA TSAM
#'
#' @param Obj A RangedSummarizedExperment generated from getSampleTileMatrix
#' @param formula Formula used for linear modeling. 
#' @param CellType The name of the celltype you wish to model. Should align with assayNames of the Obj. 
#' @param threshold A number greater than 0 and less then or equal to 1. The threshold used to determine whether to model a region or not, based on the fraction of non-zero measurements across samples at that location. 
#' @param NAtoZero Boolean, whether to convert NA region (no accessibility measurement) to zero. 
#' @param numCores Number of threads to parallelize modeling over. Default is 1. 
#'
#' @return A list of lmer model objects
#' 
#' 
#' @examples
#' \dontrun{
#'      vd_mat <- varDecomp(SampleTileObj, Donor_col = 'PTID', Time_col = 'Days')
#' )
#' }
#'
#' @export
#' 
#' 

linearModeling <- function(Obj, formula, CellType, threshold = 0, NAtoZero = FALSE, numCores = 1){

    . <- Sample <- Sample2 <- NULL
    meta1 <- as.data.frame(SummarizedExperiment::colData(Obj))

    if(length(CellType) == 1){

        mat1 <- MOCHA::getCellPopMatrix(Obj,CellType,NAtoZero = NAtoZero)

        meta <- meta1[meta1$Sample %in% colnames(mat1),]
        rowsToKeep <- rowSums(!is.na(mat1))/dim(mat1)[2] > threshold

        mat1 <- mat1[rowsToKeep,]

    }else{

        allMatrices <- do.call('cbind', SummarizedExperiment::assays(Obj))
        if(NAtoZero){
           allMatrices[is.na(allMatrices)] = 0
        }


        colnames(allMatrices) <- apply(expand.grid(names(SummarizedExperiment::assays(Obj)), unique(colnames(allMatrices))), 1, 
                                paste, collapse="__") %>% gsub(" ", "_", .)
        mat1 <- allMatrices[,colSums(allMatrices) != 0]

        meta <- parallel::mclapply(1:length(SummarizedExperiment::assays(Obj)), function(x){
        
                    meta1 %>% as.data.frame() %>%
                        dplyr::mutate(Sample2 = paste(names(SummarizedExperiment::assays(Obj))[x], Sample, sep = "__"),
                                    CellType = names(SummarizedExperiment::assays(Obj))[x]) %>%
                        dplyr::mutate(Sample2 = gsub(" ","_", Sample2), 
                                    CellType = gsub(" ", "_", CellType))
        
        }, mc.cores = numCores) %>% do.call('rbind', .)

        meta <- meta[meta$Sample2 %in% colnames(mat1),]

        rowsToKeep <- rowSums(is.na(mat1))/dim(mat1)[2] > threshold

        mat1 <- mat1[rowsToKeep,]

    }

    suppressMessages(lmem_res <- pbapply::pblapply(c(1:nrow(mat1)),
        function(x) {
            df <- data.frame(exp = as.numeric(mat1[x, ]), 
                meta, stringsAsFactors = FALSE)
            lmerTest::lmer(formula = formula, data = df)
        }, cl = numCores), classes = "message")

    names(lmem_res) = rownames(mat1)
    return(lmem_res)

}
