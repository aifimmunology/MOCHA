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

linearModeling <- function(Obj, formula, CellType, numCores = 1){

    meta <- as.data.frame(SummarizedExperiment::colData(Obj))

    if(length(CellType) == 1){

        mat1 <- MOCHA::getCellPopMatrix(Obj,CellType)

        mat1[is.na(mat1)] = 0

    }else{



    }

    meta <- meta[meta$Sample %in% colnames(mat1),]

    suppressMessages(lmem_res <- pbapply::pblapply(c(1:nrow(mat1)),
        function(x) {
            df <- data.frame(exp = as.numeric(mat1[x, ]), 
                meta, stringsAsFactors = FALSE)
            lmerTest::lmer(formula = formula, data = df)
        }, cl = numCores), classes = "message")

    return(lmem_res)

}

lmeModeling <- function(Obj, formula1 = NULL, variableList = NULL, groupCol = NULL, percNonZero = 0, numCores = 1){

    if(is.null(variableList) & is.null(formula1)){
        
        error('No variables or formula give')

    }else if(!is.null(variableList)){

        ##check if variabeList includes random affects or interactions. 
        testForm <- lapply(paste(variableList, ")",sep =''), function(x) grepl(x, as.character(formula1))  
        
    }





    if(length(SummarizedExperiment::assays(Obj)) > 1 ){

        variableList <- append(variableList, 'CellType')

        allMatrices <- do.call('cbind', SummarizedExperiment::assays(Obj))
        allMatrices[is.na(allMatrices)] = 0

        colnames(allMatrices) <- apply(expand.grid(names(SummarizedExperiment::assays(Obj)), unique(colnames(allMatrices))), 1, 
                                paste, collapse="__") %>% gsub(" ", "_", .)
        subMatrix <- allMatrices[,colSums(allMatrices) != 0]
        meta <- parallel::mclapply(1:length(SummarizedExperiment::assays(Obj)), function(x){
        
                    SummarizedExperiment::colData(Obj) %>% as.data.frame() %>%
                        dplyr::mutate(Sample = paste(names(SummarizedExperiment::assays(Obj))[x], Sample, sep = "__"),
                                    CellType = names(SummarizedExperiment::assays(Obj))[x]) %>%
                        dplyr::mutate(Sample = gsub(" ","_", Sample), 
                                    CellType = gsub(" ", "_", CellType))
        
        }, mc.cores = numCores) %>% do.call('rbind', .)
        
    }else if(length(SummarizedExperiment::assays(Obj)) == 1 ){

        allMatrices <- MOCHA::getCellPopMatrix(Obj,names(assays(Obj))[1])
        subMatrix <- allMatrices[,colSums(allMatrices) != 0]

        meta <- SummarizedExperiment::colData(Obj)

    }

    meta_f <- meta[meta$Sample %in% colnames(subMatrix),]

    #Filter out rows that are too sparse for modeling, and/or may cause issues. 
    if(is.null(groupCol)){    

        subMatrix <- subMatrix[rowSums(subMatrix != 0)/ncol(subMatrix) > percNonZero,]

    }else if(groupCol %in% colnames(meta_f)){

        groupList <- unique(meta_f[,groupCol])
        keepRows <- mclapply(groupList, function(x){

            group_tmp <- meta_f$Sample[meta_f[,groupCol] == x]
            rowSums(subMatrix[,colnames(subMatrix) %in% group_tmp] != 0)/length(group_tmp) > percNonZero

        }, mc.cores = numCores) %>% do.call(cbind,.) %>% rowSums(.)

        keptTiles <- keepRows > 1

        subMatrix <- subMatrix[keptTiles,]

    }


    form <- paste(paste("(1|", variableList, ")", sep = ""), collapse = " + ")
    form1 <- as.formula(paste("exp ~ ", form, sep = ""))

    rowN <- dim(subMatrix)[1] 
    op <- pbapply::pboptions(type = "timer")
    
    suppressMessages(lmem_res <- pbapply::pblapply(c(1:rowN),
        function(x) {
            df <- data.frame(exp = as.numeric(subMatrix[x, ]), 
                meta_f, stringsAsFactors = FALSE)
            lmem <- lmerTest::lmer(formula = form1, data = df)


            getVarDecomp(lmem, variableList, rownames(subMatrix)[x])
           

        }, cl = numCores), classes = "message")
    pbapply::pboptions(op)
    lmem_res <- do.call(rbind, lmem_res)
    lmem_res[, c(variableList, 'Residual')] <- 100 * lmem_res[,variableList]
    return(lmem_res)


}

getVarDecomp <- function(lmem, variableList, TileName){

    lmem_re <- as.data.frame(lmerTest::VarCorr(lmem))
    row.names(lmem_re) <- c(lmem_re$grp, 'Residual')

    lmem_re <- lmem_re[c(variableList,'Residual'), ]
    fix_effect <- fixef(lmem)  #get fixed effect
    lmem_re$CV <- lmem_re$sdcor/fix_effect

    normVar <- (lmem_re$vcov)/sum(lmem_re$vcov)
    names(normVar) <- row.names(lmem_re) 

    output_df <- data.frame(Tiles = TileName, Mean = mean(df$exp, na.rm = TRUE),
                Median = median(df$exp, na.rm = TRUE), SD = sd(df$exp, na.rm = TRUE),
                Max = max(df$exp, na.rm = TRUE))
    return(cbind(output_df, t(as.data.frame(normVar))))

}


varDecomp <- function(Obj, Donor_col = NULL, Time_col = NULL, Group_col = NULL, otherCols = NULL, meanThreshold=0, returnObj = FALSE, numCores = 1){


     if(is.null(Group_col) & length(SummarizedExperiment::assays(Obj)) > 1 ){

        Group_col = 'CellType'

        allMatrices <- do.call('cbind', SummarizedExperiment::assays(Obj))
        allMatrices[is.na(allMatrices)] = 0

        colnames(allMatrices) <- apply(expand.grid(names(SummarizedExperiment::assays(Obj)), unique(colnames(allMatrices))), 1, 
                                paste, collapse="__") %>% gsub(" ", "_", .)
        subMatrix <- allMatrices[,colSums(allMatrices) != 0]
        meta <- parallel::mclapply(1:length(SummarizedExperiment::assays(Obj)), function(x){
        
                    SummarizedExperiment::colData(Obj) %>% as.data.frame() %>%
                        dplyr::mutate(Sample = paste(names(SummarizedExperiment::assays(Obj))[x], Sample, sep = "__"),
                                    CellType = names(SummarizedExperiment::assays(Obj))[x]) %>%
                        dplyr::mutate(Sample = gsub(" ","_", Sample), 
                                    CellType = gsub(" ", "_", CellType))
        
        }, mc.cores = numCores) %>% do.call('rbind', .)
        
    }else if(!is.null(Group_col) & length(SummarizedExperiment::assays(Obj)) == 1 ){

        allMatrices <- SummarizedExperiment::assays(Obj)[[1]]
        subMatrix <- allMatrices[,colSums(allMatrices) != 0]

        meta <- SummarizedExperiment::colData(Obj)

    }else{

        error("You cannot model cell type variance and group variance at the same time. Please choose one, and subset your SampleTileObj accordingly.")

    }

    meta_f <- meta[meta$Sample %in% colnames(subMatrix),]

    if(!is.null(Time_col) & !is.null(Donor_col)){

        featureSet1 <- append(otherCols, c(Donor_col, Time_col, Group_col))
        palmo_obj <- PALMO::createPALMOobject(anndata = meta_f, data= subMatrix)
        palmo_obj<- PALMO::annotateMetadata(data_object=palmo_obj, sample_column= "Sample", donor_column= Donor_col, time_column= Time_col, group_column = Group_col)

        palmo_obj <- PALMO::mergePALMOdata(data_object=palmo_obj, datatype="bulk")


    }else if(is.null(Time_col)  & !is.null(Donor_col)){

        featureSet1 <- append(otherCols, c(Donor_col, Group_col))
        meta_f$Time = seq(1,length(meta_f$Sample),1)
        palmo_obj <- PALMO::createPALMOobject(anndata = meta_f, data= subMatrix)
        Time_col = "Time"
        palmo_obj<- PALMO::annotateMetadata(data_object=palmo_obj, sample_column= "Sample", donor_column= Donor_col, time_column= Time_col, group_column = Group_col)
        palmo_obj <- PALMO::mergePALMOdata(data_object=palmo_obj, datatype="bulk")

    }else if(is.null(Time_col) & is.null(Donor_col)){

        featureSet1 <- append(otherCols, c(Group_col))
        meta_f$Time = seq(1,length(meta_f$Sample),1)
        meta_f$Donor = rep(2, seq(1,length(meta_f$Sample),1)
        return(meta_f)
        palmo_obj <- PALMO::createPALMOobject(anndata = meta_f, data= subMatrix)
        palmo_obj<- PALMO::annotateMetadata(data_object=palmo_obj, sample_column= "Sample", donor_column= "Donor", time_column= "Time", group_column = Group_col)
        palmo_obj <- PALMO::mergePALMOdata(data_object=palmo_obj, datatype="bulk")

    }

    palmo_obj <- lmeVariance(data_object=palmo_obj,
            featureSet=featureSet1,
            meanThreshold=0,
            cl=numCores,
            fileName='scATAC')

    if(returnObj){

        return(palmo_obj)

    }else{
        return(palmo_obj@result$variance_decomposition)
    }

}