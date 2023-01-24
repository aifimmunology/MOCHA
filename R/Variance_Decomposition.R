#' @title \code{varDecomp}
#'
#' @description \code{varDecomp} Runs PALMO to model variance decomposition along the variables needed. 
#'
#' @param Obj A RangedSummarizedExperment generated from getSampleTileMatrix
#' @param Donor_col The metadata column with donor information.
#' @param Time_col The metadata column with longitudinal information
#' @param Group_col The metadata column with group information. Default is NULL, in which case it will attempt to use celltype as a proxy. Do not run group and cell type at the same time. 
#' @param returnObj Boolean flag for whether or not you want the full PALMO object. Default is FALSE, in which case it returns the variance decomposition data.frame. 
#'
#' @return variance decomposition matrix
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

linearModeling <- function(Obj, formula, CellType, rowsToKeep = NA, NAtoZero = FALSE, numCores = 1){

    meta1 <- as.data.frame(SummarizedExperiment::colData(Obj))

    if(length(CellType) == 1){

        mat1 <- MOCHA::getCellPopMatrix(Obj,CellType,NAtoZero = NAtoZero)

        meta <- meta1[meta1$Sample %in% colnames(mat1),]
        
        if(!all(is.na(rowsToKeep))){

            mat1 <- mat1[rowsToKeep,]

        }

        

    }else{

        allMatrices <- do.call('cbind', SummarizedExperiment::assays(Obj))
        if(NAtoZero){
           allMatrices[is.na(allMatrices)] = 0
        }


        colnames(allMatrices) <- apply(expand.grid(colnames(Obj),names(SummarizedExperiment::assays(Obj))), 1, 
                                paste, collapse="__") %>% gsub(" ", "_", .)
        mat1 <- allMatrices[,colSums(allMatrices) != 0]
        rm(allMatrices)
        
        gc()

        cl <- parallel::makeCluster(numCores)
        parallel::clusterExport(cl, varlist = c('meta1', 'Obj'), envir = environment())

        meta <- pbapply::pblapply(1:length(SummarizedExperiment::assays(Obj)), function(x){
        
                    meta1 %>% as.data.frame() %>%
                        dplyr::mutate(Sample2 = paste(names(SummarizedExperiment::assays(Obj))[x], Sample, sep = "__"),
                                    CellType = names(SummarizedExperiment::assays(Obj))[x]) %>%
                        dplyr::mutate(Sample2 = gsub(" ","_", Sample2), 
                                    CellType = gsub(" ", "_", CellType))
        
        }, cl = cl) %>% do.call('rbind', .)

        parallel::stopCluster(cl)

        meta <- meta[meta$Sample2 %in% colnames(mat1),]

        if(!all(is.na(rowsToKeep))){

            mat1 <- mat1[rowsToKeep,]

        }

        mat1 <- mat1[rowsToKeep,]

    }

    cl <- parallel::makeCluster(numCores)

    parallel::clusterExport(cl=cl, varlist=c("meta", "formula", "mat1"), envir=environment())

    suppressMessages(lmem_res <- pbapply::pblapply(c(1:nrow(mat1)),
        function(x) {
            df <- data.frame(exp = as.numeric(mat1[x,]), 
                meta, stringsAsFactors = FALSE)
            lmerTest::lmer(formula = formula, data = df)
        }, cl = cl), classes = "message")

    names(lmem_res) = rownames(mat1)

    stopCluster(cl)

    return(lmem_res)

}


calculateVarDecomp <- function(Obj, CellType, variableList, rowsToKeep = NA, NAtoZero = FALSE, numCores = 1){

    varForm <- paste0(unlist(lapply(variableList, function(x) paste('(1|',x,')',sep=''))), collapse = ' + ')
    formula1 <- as.formula(paste('exp ~ ',varForm, sep = ''))

    linRes <- linearModeling(Obj, formula = formula1, CellType, rowsToKeep, NAtoZero, numCores)
    
    varDecomp <- getVarDecomp(linRes, numCores)

    return(varDecomp)
}



getVarDecomp <- function(lmemList, numCores = 1){

    variableList <- names(summary(lmemList[[1]])$ngrps)

    cl <- parallel::makeCluster(numCores)

    parallel::clusterExport(cl=cl, varlist=c("lmemList"), envir=environment())

    varDecompRes <- suppressMessages(lmem_res <- pbapply::pblapply(seq_along(lmemList),
        function(x) {

            extractVarDecomp(lmemList[[x]], variableList)

     }, cl = cl), classes = "message") %>% do.call('rbind',.)

    stopCluster(cl)

    gc()

    output_df <- cbind(data.frame(Tiles = names(lmemList)), varDecompRes)
    return(output_df)

}


extractVarDecomp <- function(lmem1, variableList){

    lmem_re <- as.data.frame(lme4::VarCorr(lmem1))
    row.names(lmem_re) <- c(lmem_re$grp)

    lmem_re <- lmem_re[c(variableList,'Residual'), ]
    fix_effect <- lme4::fixef(lmem1)  #get fixed effect
    lmem_re$CV <- lmem_re$sdcor/fix_effect

    normVar <- (lmem_re$vcov)/sum(lmem_re$vcov)
    names(normVar) <- row.names(lmem_re) 

    return(normVar)

}