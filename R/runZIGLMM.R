#' @title Run Zero-inflated Generalized Linear Mixed Modeling on pseudobulked scATAC data
#'
#' @description \code{runZIGLMM} Runs linear mixed-effects modeling for
#'   zero-inflated data using \code{\link[glmmTMB]{glmmTMB}}. 
#' 
#' @param TSAM_Object A SummarizedExperiment object generated from
#'   getSampleTileMatrix. 
#' @param cellTypeName Name of a cell type(s), or 'all'. The function will combine the cell types mentioned into one matrix before running the model.
#' @param continuousFormula The formula for the continuous data that should be used within glmmTMB. It should be in the
#'   format (exp ~ factors). All factors must be found in column names
#'   of the TSAM_Object metadata, except for CellType, FragNumber and CellCount, which will be extracted from the TSAM_Object.
#'   modelFormula must start with 'exp' as the response.
#'   See \link[glmmTMB]{glmmTMB}.
#' @param ziformula The formula for the zero-inflated data that should be used within glmmTMB. It should be in the
#'   format ( ~ factors). All factors must be found in column names
#'   of the TSAM_Object colData metadata, except for CellType, FragNumber and CellCount, which will be extracted from the TSAM_Object.
#' @param zi_threshold Zero-inflated threshold ( range = 0-1), representing the fraction of samples with zeros. At or above this threshold, the zero-inflated modeling kicks in.
#' @param initialSampling Size of data to use for pilot
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @param numCores integer. Number of cores to parallelize across.
#'
#' @return results a SummarizedExperiment containing LMEM results
#'
#'
#'
#' @examples
#' \dontrun{
#'   modelList <- runZIGLMM(STM[c(1:1000),], 
#'                  cellTypeName = 'CD16 Mono',
#'                  continuousFormula = exp~ Age + Sex + days_since_symptoms + (1|PTID), 
#'                  ziformula = ~ FragNumber + Age, 
#'                  verbose = TRUE, 
#'                  numCores = 35 )
#' }
#'
#' @export
#' 
runZIGLMM <- function(TSAM_Object,
                      cellTypeName = 'all',
                      continuousFormula = NULL,
                      ziformula = NULL,
                      zi_threshold = 0,
                      initialSampling = 5,
                      verbose = FALSE,
                      numCores = 1) {
  if (any(c(class(continuousFormula), class(ziformula)) != "formula")) {
    stop("continuousFormula and/or ziformula was not provided as a formula.")
  }

  if (zi_threshold < 0 | zi_threshold > 1 | ! is.numeric(zi_threshold)) {
    stop("zi_threshold must be between 0 and 1.")
  }

  if (is.null(cellTypeName)) {
    stop("No cell type name was provided.")
  } else if (all(tolower(cellTypeName) == 'all')) {

    #Merge all together. 
    newObj <- combineSampleTileMatrix(TSAM_Object)

  } else if (all(cellTypeName %in% names(SummarizedExperiment::assays(TSAM_Object)))) {

    #Subset down to just those
    newObj <- combineSampleTileMatrix(subsetMOCHAObject(TSAM_Object, subsetBy = 'celltype', groupList = cellTypeName, subsetPeaks = TRUE))

  } else {
  
    stop("Error around cell type name. Some or all were not found within TSAM_Object.")

  }
  
  modelingData <- log2(SummarizedExperiment::assays(newObj)[['counts']]+1)
  MetaDF <- as.data.frame(SummarizedExperiment::colData(newObj))

  if (!all(all.vars(continuousFormula) %in% c("exp", colnames(MetaDF)))) {
    stop("Model formula is not in the correct format (exp ~ factors) or model factors are not found in column names of metadata within the TSAM_Object.")
  }

  if (!all(all.vars(ziformula) %in% colnames(MetaDF))) {
    stop("factors from the ziformula were not found in the metadata.")
  }

  variableList <- c(all.vars(continuousFormula)[all.vars(continuousFormula) != "exp"], all.vars(ziformula))

  MetaDF <- dplyr::filter(MetaDF, Sample %in% colnames(modelingData))
  modelingData <- modelingData[, match(colnames(modelingData), MetaDF$Sample)]


  # Subset metadata to just the variables needed. This minimizes overhead for parallelization
  MetaDF <- MetaDF[, colnames(MetaDF) %in% c("Sample", variableList)]


  if (verbose) {
    message("Running a quick test.")
  }


  # Generate pilot data for the null data.frame
  pilotIndices <- sample(x = 1:dim(modelingData)[1], size = initialSampling, replace = FALSE)
  modelList <- pbapply::pblapply(X = pilotIndices, function(x) {
    df <- data.frame(
      exp = as.numeric(modelingData[x, ]),
      MetaDF, stringsAsFactors = FALSE
    )

    tryCatch(
      {
        glmmTMB::glmmTMB(continuousFormula,
          ziformula = ziformula,
          data = df,
          family = stats::gaussian(),
          REML = TRUE
        )
      },
      error = function(e) {
        NA
      }
    )
  }, cl = NULL)

  if (all(is.na(unlist(modelList)))) {
    stop("For the initial sampling, every test model failed. Reconsider modelFormula or increase initial sampling size.")
  } else {
    idx <- which(!is.na(unlist(modelList)))
    coeffList <- lapply(1:length(idx), function(x){
      tryCatch(
      {
         summary(modelList[[x]])$coefficients
      }, error = function(e){ NA})
    })
    if (all(is.na(unlist(coeffList)))) {
      stop("For the initial sampling, every test model failed. Reconsider modelFormula or increase initial sampling size.")
    }
    coeff2 <- coeffList[!is.na(unlist(coeffList))]
    nullDF <- lapply(coeff2[[1]], as.data.frame)
    nullDF$cond[!is.na(nullDF$cond)] <- NA
    nullDF$zi[!is.na(nullDF$zi)] <- NA
    rm(modelList)
  }

  if (verbose) {
    message("Modeling results.")
  }


  # Make your clusters for efficient parallelization
  cl <- parallel::makeCluster(numCores)
  parallel::clusterEvalQ(cl, {
    library(glmmTMB)
  })
  parallel::clusterExport(
    cl = cl, varlist = c("continuousFormula", "ziformula", "modelingData", "MetaDF", "individualZIGLMM", "nullDF","zi_threshold"),
    envir = environment()
  )
  coeffList <- pbapply::pblapply(cl = cl, X = rownames(modelingData), individualZIGLMM)
  parallel::stopCluster(cl)

  output_list <- lapply(list("cond", "zi"), function(varType) {
    if (verbose) {
      message(stringr::str_interp("Extracting coefficients for the ${varType} component"))
    }
    valuesToExtract <- c("Estimate", "Pr(>|z|)", "Std. Error")

    dfList <- lapply(valuesToExtract, function(valType) {
      tmpDF <- do.call("rbind", pbapply::pblapply(X = coeffList, function(x) {
        extractVariable(x, varType, valType, nullDF)
      }, cl = NULL))
      if (!is.null(tmpDF)) {
        rownames(tmpDF) <- rownames(modelingData)
      }
      tmpDF
    })

    names(dfList) <- c("Slopes", "Significance", "StdError")
    dfList
  })

  names(output_list) <- c("cond", "zi")

  combinedList <- lapply(c("Slopes", "Significance", "StdError"), function(x) {
    cond1 <- output_list[["cond"]][[x]]
    zi1 <- output_list[["zi"]][[x]]

    colnames(zi1) <- paste("ZI_", colnames(zi1), sep = "")
    cbind(cond1, zi1)
  })

  names(combinedList) <- c("Slopes", "Significance", "StdError")

  results <- SummarizedExperiment::SummarizedExperiment(
    combinedList,
    rowRanges = SummarizedExperiment::rowRanges(newObj),
    metadata = newObj@metadata
  )

  return(results)
}

extractVariable <- function(varList, varType, variable, nullDF) {
  varDF <- varList[[varType]]
  if (dim(varDF)[2] > 0) {
    val_tmp <- unlist(varDF[, variable])
    names(val_tmp) <- rownames(varDF)
    val_tmp
  } else {
    varDF <- nullDF[[varType]]
    val_tmp <- unlist(varDF[, variable])
    names(val_tmp) <- rownames(varDF)
    val_tmp
  }
}


#' @title \code{IndividualZIGLMM}
#'
#' @description \code{IndividualZIGLMM} Runs zero-inflated linear modeling on data provided. Written for efficient parallelization.
#'
#' @param refList. A list where the first index is a data.frame to use for modeling, and the second is the formula for modeling.
#'
#' @return A linear model

#' @export
#'
#'
individualZIGLMM <- function(x) {
  df <- data.frame(
    exp = as.numeric(modelingData[x, ]),
    MetaDF, stringsAsFactors = FALSE
  )

  output_vector <- tryCatch(
    {
     
      if(sum(df$exp == 0)/length(df$exp) == 0){
        modelRes <- glmmTMB::glmmTMB(continuousFormula,
          ziformula = ~ 0,
          data = df,
          family = stats::gaussian(),
          REML = TRUE
        )
       
      }else if(sum(df$exp == 0)/length(df$exp) <= zi_threshold){
        df$exp[df$exp == 0] = NA
        modelRes <- glmmTMB::glmmTMB(continuousFormula,
          ziformula = ~ 0,
          data = df,
          family = stats::gaussian(),
          REML = TRUE
        )
       
      }else {
        modelRes <- glmmTMB::glmmTMB(continuousFormula,
          ziformula = ziformula,
          data = df,
          family = stats::gaussian(),
          REML = TRUE
        )
      }
      
      lapply(summary(modelRes)$coefficients, as.data.frame)
    },
    error = function(e) {
      nullDF
    }
  )
  return(output_vector)
}


#' @title Run a test case for Zero-inflated Generalized Linear Mixed Modeling on pseudobulked scATAC data
#'
#' @description \code{pilotZIGLMM} ode for testing out formulas on the data. Runs a given formula on a subset of the data, and returns the model results. This is meant to help during the model selection process. \code{\link[glmmTMB]{glmmTMB}}. 
#'
#' @param TSAM_Object A SummarizedExperiment object generated from
#'   getSampleTileMatrix. 
#' @param cellTypeName Name of a cell type(s), or 'all'. The function will combine the cell types mentioned into one matrix before running the model.
#' @param continuousFormula The formula for the continuous data that should be used within glmmTMB. It should be in the
#'   format (exp ~ factors). All factors must be found in column names
#'   of the TSAM_Object metadata, except for CellType, FragNumber and CellCount, which will be extracted from the TSAM_Object.
#'   modelFormula must start with 'exp' as the response.
#'   See \link[glmmTMB]{glmmTMB}.
#' @param ziformula The formula for the zero-inflated data that should be used within glmmTMB. It should be in the
#'   format ( ~ factors). All factors must be found in column names
#'   of the TSAM_Object colData metadata, except for CellType, FragNumber and CellCount, which will be extracted from the TSAM_Object.
#' @param zi_threshold Zero-inflated threshold ( range = 0-1), representing the fraction of samples with zeros. At or above this threshold, the zero-inflated modeling kicks in.
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @param pilotIndices integer. Specific locations to test within the peakset of the cell type chosen. 
#'
#' @return results a SummarizedExperiment containing LMEM results
#'
#'
#'
#' @examples
#' \dontrun{
#'   modelList <- pilotZIGLMM(STM, 
#'                            'CD16 Mono',
#'                            exp~ Age + Sex + days_since_symptoms + (1|PTID),
#'                             ~ FragNumber, verbose = TRUE )
#' }
#'
#' @export
#' 

pilotZIGLMM <- function(TSAM_Object,
                        cellTypeName = NULL,
                        continuousFormula = NULL,
                        ziformula = NULL,
                        zi_threshold = 0,
                        verbose = FALSE,
                        pilotIndices = 1:10) {
  if (any(c(class(continuousFormula), class(ziformula)) != "formula")) {
    stop("continuousFormula and/or ziformula was not provided as a formula.")
  }

  if (zi_threshold < 0 | zi_threshold > 1 | ! is.numeric(zi_threshold)) {
    stop("zi_threshold must be between 0 and 1.")
  }

  if (is.null(cellTypeName)) {
    stop("No cell type name was provided.")
  } else if (length(cellTypeName) > 1) {
    stop("Please provide only one string within cellTypeName. If you want to run over multiple cell types, please use combineSampleTileMatrix() to generate a new object, and use that object instead, with cellTypeName = 'counts'")
  } else if (!cellTypeName %in% names(SummarizedExperiment::assays(TSAM_Object))) {
    stop("Cell type name not found within TSAM_Object.")
  }

  newObj <- combineSampleTileMatrix(subsetMOCHAObject(TSAM_Object, subsetBy = 'celltype', groupList = cellTypeName, subsetPeaks = TRUE))
  modelingData <- log2(SummarizedExperiment::assays(newObj)[['counts']]+1)
  MetaDF <- as.data.frame(SummarizedExperiment::colData(newObj))

  if (!all(all.vars(continuousFormula) %in% c("exp", colnames(MetaDF)))) {
    stop("Model formula is not in the correct format (exp ~ factors) or model factors are not found in column names of metadata within the TSAM_Object.")
  }

  if (!all(all.vars(ziformula) %in% c(colnames(MetaDF))) & length(all.vars(ziformula)) > 0) {
    stop("factors from the ziformula were not found in the metadata.")
  }

  variableList <- c(all.vars(continuousFormula)[all.vars(continuousFormula) != "exp"], all.vars(ziformula))

  MetaDF <- dplyr::filter(MetaDF, Sample %in% colnames(modelingData))
  pilotNames <- rownames(modelingData)[pilotIndices]
  modelingData <- modelingData[pilotNames, match(colnames(modelingData), MetaDF$Sample)]

  # Subset metadata to just the variables needed. This minimizes overhead for parallelization
  MetaDF <- MetaDF[, colnames(MetaDF) %in% c("Sample", variableList)]

  if (verbose) {
    message("Modeling results.")
  }

  # Make your clusters for efficient parallelization
  modelList <- pbapply::pblapply(X = pilotNames, function(x) {
  
    df <- data.frame(
      exp = as.numeric(modelingData[x, ]),
      MetaDF, stringsAsFactors = FALSE
    )

    if(sum(df$exp == 0)/length(df$exp) == 0){
        modelRes <- glmmTMB::glmmTMB(continuousFormula,
          ziformula = ~ 0,
          data = df,
          family = stats::gaussian(),
          REML = TRUE
        )
       
      }else if(sum(df$exp == 0)/length(df$exp) <= zi_threshold){
        df$exp[df$exp == 0] = NA
        modelRes <- glmmTMB::glmmTMB(continuousFormula,
          ziformula = ~ 0,
          data = df,
          family = stats::gaussian(),
          REML = TRUE
        )
       
      }else {
        modelRes <- glmmTMB::glmmTMB(continuousFormula,
          ziformula = ziformula,
          data = df,
          family = stats::gaussian(),
          REML = TRUE
        )
      }

  }, cl = NULL)
  names(modelList) <- pilotNames

  return(modelList)
}


#' @title getModelValues from runZIGLMM output. 
#'
#' @description \code{getModelValues} Pull out a data.frame of model values (slope, significance, and std.error) for a given factor from the SummarizedExperiment output of runZIGLMM.
#' @param object A SummarizedExperiment object generated from runZIGLMM. 
#' @param specificVariable A string, describing the factor of influence. 
#'
#' @return A data.frame of slopes, significance, and standard error for one factor. 
#'
#'
#'
#' @examples
#' \dontrun{
#'   age_df <- getModelValues(runZIGLMM_output, 'Age')
#' }
#'
#' @export
#' 




getModelValues <- function(object, specificVariable){
    
    slopes = SummarizedExperiment::assays(object)[['Slopes']]
    significance = SummarizedExperiment::assays(object)[['Significance']]
    if(length(specificVariable) > 1){
        stop('Cannot provide more than one value to specificVariable.')
    }

    df <- data.frame('Element' = rownames(slopes),
                'Estimate' = slopes[,specificVariable],
                'PValue' = significance[,specificVariable])
    return(df)
}


