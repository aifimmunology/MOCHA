## To do: Implement ZI options, and figure out how to handle Log2 or not log2. 

## Will take data (metadata and modeling data), and use it to model a given formula.
## It then returns a data.frame of intercept, slope, significance, and (residual?) 
## for each row of modelingData

## formula must be in the form exp ~ <variables>

## Initial sampling: When the model fails to converge, it needs to return a default data.frame with NAs that is the same size. 
## The initial sampling provides a blueprint for how big that NA data.frame should be. The default is that it'll run 5 models, and if they are all NA, through an error. 
## This function should only be used on continuous, non-zero inflated data. 

runLMEM <- function(TSAM_Object,
                        cellTypeName = NULL,
                        modelFormula = NULL, 
                        initialSampling = 5,
                        verbose = FALSE, 
                        numCores = 1){

    
    if (is.null(cellTypeName)) {
        stop("No cell type name was provided.")
    }else if(length(cellTypeName) > 1){
        stop("Please provide only one string within cellTypeName. If you want to run over multiple cell types, please use combineSampleTileMatrix() to generate a new object, and use that object instead, with cellTypeName = 'counts'")
    }else if(!cellTypeName %in% names(SummarizedExperiment::assays(TSAM_Object))){
        stop("cellTypeName was not found within TSAM_Object.")
    }

    modelingData <- as.data.frame(getCellPopMatrix(TSAM_Object, cellPopulation = cellTypeName, NAtoZero = TRUE))
    MetaDF <- as.data.frame(SummarizedExperiment::colData(TSAM_Object))

    if (is.null(modelFormula)) {
        stop("No formula was provided.")
    }else if(!all(all.vars(modelFormula) %in% c('~','exp',colnames(MetaDF)))){
        stop("Model formula is not in the correct format (exp ~ factors) or model factors are not found in column names of metadata within the TSAM_Object.")
    }

    variableList <- all.vars(modelFormula)[all.vars(modelFormula) != 'exp']

    MetaDF <- dplyr::filter(MetaDF, Sample %in% colnames(modelingData))
    modelingData <- modelingData[,match(colnames(modelingData), MetaDF$Sample)]
    

    #Subset metadata to just the variables needed. This minimizes overhead for parallelization
    MetaDF <- MetaDF[, colnames(MetaDF) %in% c('Sample', variableList)]

    
    if(verbose){
        message('Running a quick test.')
    }


    #Generate pilot data for the null data.frame
    pilotIndices <- sample(x = 1:dim(modelingData)[1], size= initialSampling, replace = FALSE)
    modelList <- pbapply::pblapply(X = pilotIndices, function(x){

        df <-  data.frame(exp = as.numeric(modelingData[x,]), 
                MetaDF, stringsAsFactors = FALSE)

        tryCatch({
            lmerTest::lmer(formula = modelFormula, data =  df)
        }, error = function(e){ NA})

    }, cl = NULL)

    if(all(is.na(unlist(modelList)))){
        stop('For the initial sampling, every test model failed. Reconsider modelFormula or increase initial sampling size.')
    }else{

        idx <- which(!is.na(unlist(modelList)))
        nullDF <- as.data.frame(summary(modelList[[idx[1]]])$coefficients)
        nullDF[!is.na(nullDF)] <- NA
        rm(modelList)
    }

    
    if(verbose){
        message('Modeling results.')
    }

    # Make your clusters for efficient parallelization

    cl <- parallel::makeCluster(numCores)
    parallel::clusterExport(cl=cl, varlist = c('modelFormula', 'modelingData', 'MetaDF','individualLMEM','nullDF'), 
                        envir = environment())
    coeffList <- pbapply::pblapply(cl = cl, X = rownames(modelingData), individualLMEM)
    parallel::stopCluster(cl)

    if(verbose){
        message('Reorganizing coefficients.')
    }

    slopes <-  do.call('rbind', pbapply::pblapply(X = coeffList, function(x){ 
                slope_tmp <- x$Estimate
                names(slope_tmp) <- rownames(x)
                slope_tmp

    },cl = NULL))

    significance <-  do.call('rbind', pbapply::pblapply(X = coeffList, function(x){ 
                sig_tmp <- x$'Pr(>|t|)'
                names(sig_tmp) <- rownames(x)
                sig_tmp

    },cl = NULL))
   
    stdError <-  do.call('rbind', pbapply::pblapply(X = coeffList, function(x){ 
                error_tmp <- x$'Std. Error'
                names(error_tmp) <- rownames(x)
                error_tmp

    },cl = NULL))

    rownames(stdError) <- rownames(significance) <- rownames(slopes) <- rownames(modelingData)


    output_list <- list('Slopes' = slopes, 'Significance' = significance, StdError = stdError)

    return(output_list)
}


#' @title \code{IndividualLMEM}
#'
#' @description \code{IndividualLMEM} Runs linear modeling on data provided. Written for efficient parallelization.
#'
#' @param refList. A list where the first index is a data.frame to use for modeling, and the second is the formula for modeling.  
#'
#' @return A linear model

#' @export
#' 
#' 
individualLMEM <- function(x){

    df <-  data.frame(exp = as.numeric(modelingData[x,]), 
                MetaDF, stringsAsFactors = FALSE)

    output_vector <- tryCatch({
        modelRes <- lmerTest::lmer(formula = modelFormula, data =  df)
        as.data.frame(summary(modelRes)$coefficients)
    }, error = function(e){
            nullDF
      })
    return(output_vector)
}


## Code for testing out formulas on the data. Runs a given formula on a subset of the data, and returns the model results. 
## This is meant to help during the model selection process. 

pilotLMEM <- function(TSAM_Object,
                        cellTypeName = NULL,
                        modelFormula = NULL, 
                        ZI = FALSE,
                        verbose = FALSE,
                        pilotIndices = 1:10){

    if (is.null(cellTypeName)) {
        stop("No cell type name was provided.")
    }else if(length(cellTypeName) > 1){
        stop("Please provide only one string within cellTypeName. If you want to run over multiple cell types, please use combineSampleTileMatrix() to generate a new object, and use that object instead, with cellTypeName = 'counts'")
    }else if(!cellTypeName %in% names(SummarizedExperiment::assays(TSAM_Object))){
        stop("cellTypeName was not found within TSAM_Object.")
    }

    modelingData <- as.data.frame(getCellPopMatrix(TSAM_Object, cellPopulation = cellTypeName, NAtoZero = TRUE))
    MetaDF <- as.data.frame(SummarizedExperiment::colData(TSAM_Object))

    if (is.null(modelFormula)) {
        stop("No formula was provided.")
    }else if(!all(all.vars(modelFormula) %in% c('~','exp',colnames(MetaDF)))){
        stop("Model formula is not in the correct format (exp ~ factors) or model factors are not found in column names of metadata within the TSAM_Object.")
    }

    variableList <- all.vars(modelFormula)[all.vars(modelFormula) != 'exp']

    MetaDF <- dplyr::filter(MetaDF, Sample %in% colnames(modelingData))
    modelingData <- modelingData[pilotIndices,match(colnames(modelingData), MetaDF$Sample)]
    

    #Subset metadata to just the variables needed. This minimizes overhead for parallelization
    MetaDF <- MetaDF[, colnames(MetaDF) %in% c('Sample', variableList)]

    if(verbose){
        message('Modeling results.')
    }

    # Make your clusters for efficient parallelization

    modelList <- pbapply::pblapply(X = pilotIndices, function(x){

        df <-  data.frame(exp = as.numeric(modelingData[x,]), 
                MetaDF, stringsAsFactors = FALSE)

        lmerTest::lmer(formula = modelFormula, data =  df)
    }, cl = NULL)

    return(modelList)
}