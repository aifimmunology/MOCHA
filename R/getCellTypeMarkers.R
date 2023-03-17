


## Code that both runs a linear model across variables, and extracts the variance. This helps keep memory lower. 


getCellTypeMarkers <- function(TSAM_Object, verbose = FALSE, numCores = 1){


    fullObj <- combineSampleTileMatrix(TSAM_Object, NAtoZero = TRUE, verbose = verbose)

    CountDF <- SummarizedExperiment::assays(fullObj)[[1]]

    MetaDF <- SummarizedExperiment::colData(fullObj)

    variableList <- c('CellType', 'Freq')

    #Generate formula
    varForm <- paste0(unlist(lapply(variableList, function(x) paste('(1|',x,')',sep=''))), collapse = ' + ')
    formula1 <- as.formula(paste('exp ~ ',varForm, sep = ''))

    MetaDF <- dplyr::filter(as.data.frame(MetaDF), Sample %in% colnames(CountDF))
    CountDF <- CountDF[,match(colnames(CountDF), MetaDF$Sample)]
    
    #Ensure that MetaDF and CountDF are aligned.
    if(!all(colnames(CountDF) == MetaDF$Sample)){
        stop('Samples in CountDF are not aligned with MetaDF')
    }

    #Ensure variable list is found in MetaDF
    if(!all(variableList %in% colnames(MetaDF))){
        stop('variableList includes variables not found in MetaDF.')
    }

    #Subset metadata to just the variables needed. This minimizes overhead for parallelization
    MetaDF <- MetaDF[, colnames(MetaDF) %in% c('Sample', variableList)]

    if(verbose){
        message('Modeling results & Extracting Variance Decomposition')
    }

    # Make your clusters for efficient parallelization
    #return(list(formula1, CountDF, MetaDF, getIndividualVariance))
    cl <- parallel::makeCluster(numCores)
    parallel::clusterExport(cl=cl, varlist = c('formula1', 'CountDF', 'MetaDF','getIndividualVariance'), 
                        envir = environment())
    decompList <- pbapply::pblapply(cl = cl, X = rownames(CountDF), getIndividualVariance)
    parallel::stopCluster(cl)
    return(decompList)

    varDecomp <- do.call('rbind', decompList)
    #colnames(varDecomp) <- c(variableList, 'Residual')
    output_df <- cbind(data.frame(Variables = rownames(CountDF)), varDecomp)

    return(output_df)
}


#' @title \code{runModeling}
#'
#' @description \code{runModeling} Runs linear modeling on data provided. Written for efficient parallelization.
#'
#' @param refList. A list where the first index is a data.frame to use for modeling, and the second is the formula for modeling.  
#'
#' @return A linear model

#' @export
#' 
#' 
getIndividualVariance <- function(x){

    df <-  data.frame(exp = as.numeric(CountDF[x,]), 
                MetaDF, stringsAsFactors = FALSE)

    modelRes <- tryCatch({lme4::lmer(formula = formula1, data =  df)}, 
                                error = function(e){NA})
    variableList <-  all.vars(formula1[-1])

    output_vector <- tryCatch({
        if(!is.na(modelRes)){
            ## Extract variance decomposition.
            lmem_re <- as.data.frame(lme4::VarCorr(modelRes))
            row.names(lmem_re) <- c(lmem_re$grp)
        
            lmem_re <- lmem_re[c(variableList,'Residual'), ]

            VarDecomp <- (lmem_re$vcov)/sum(lmem_re$vcov)

        }else{
            VarDecomp <- rep(NA, length(variableList) + 1)
        }

        VarDecomp
    }, 
      error = function(e){rep(NA, length(variableList) + 1)})
    return(output_vector)
}

tmpM <- getZIVariance(rownames(CountDF)[1])
tmp2 <- summary(tmpM)

getZIVariance <- function(x){

    df <-  data.frame(exp = as.numeric(CountDF[x,]), 
                MetaDF, stringsAsFactors = FALSE)

    modelRes <- tryCatch({
        glmmTMB::glmmTMB(formula1, 
            family = gaussian(),
            ziformula = ~ CellType + Freq,
            data =  df)
        }, 
        error = function(e){NA})
    return(modelRes)
    variableList <-  all.vars(formula1[-1])

    output_vector <- tryCatch({
        if(!is.na(modelRes)){
            ## Extract variance decomposition.
            lmem_re <- as.data.frame(glmmTMB::VarCorr(modelRes))
            row.names(lmem_re) <- c(lmem_re$grp)
        
            lmem_re <- lmem_re[c(variableList,'Residual'), ]

            VarDecomp <- (lmem_re$vcov)/sum(lmem_re$vcov)

        }else{
            VarDecomp <- rep(NA, length(variableList) + 1)
        }

        VarDecomp
    }, 
      error = function(e){rep(NA, length(variableList) + 1)})
    return(output_vector)
}

