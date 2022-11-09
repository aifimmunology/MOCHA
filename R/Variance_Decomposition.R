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

varDecomp <- function(Obj, Donor_col, Time_col = NULL, Group_col = NULL, meanThreshold=0, returnObj = FALSE, numCores = 1){


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

    if(!is.null(Time_col)){

        featureSet1 <- c(Donor_col, Time_col, Group_col) 
        palmo_obj <- PALMO::createPALMOobject(anndata = meta_f, data= subMatrix)
        palmo_obj<- PALMO::annotateMetadata(data_object=palmo_obj, sample_column= "Sample", donor_column= Donor_col, time_column= Time_col, group_column = Group_col)

        palmo_obj <- PALMO::mergePALMOdata(data_object=palmo_obj, datatype="bulk")


    }else{

        featureSet1 <- c(Donor_col, Group_col) 
        meta_f$Time = rep(1, length(meta_f$Sample))
        #return(list(meta_f, subMatrix))
        palmo_obj <- PALMO::createPALMOobject(anndata = meta_f, data= subMatrix)
        palmo_obj<- PALMO::annotateMetadata(data_object=palmo_obj, sample_column= "Sample", donor_column= Donor_col, time_column= "Time", group_column = Group_col)
        palmo_obj <- PALMO::mergePALMOdata(data_object=palmo_obj, datatype="bulk")

    }

    palmo_obj <- lmeVariance(data_object=palmo_obj,
            featureSet=featureSet1 ,
            meanThreshold=0,
            cl=numCores,
            fileName='scATAC')

    if(returnObj){

        return(palmo_obj)

    }else{
        return(palmo_obj@result$variance_decomposition)
    }

}