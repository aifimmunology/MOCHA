## To do: Implement ZI options, and figure out how to handle Log2 or not log2.

## Will take data (metadata and modeling data), and use it to model a given formula.
## It then returns a data.frame of intercept, slope, significance, and (residual?)
## for each row of modelingData

## formula must be in the form exp ~ <variables>


# example: runZIGLMM(STM[c(1:1000),], 'CD16 Mono',exp~ Age + Sex + days_since_symptoms + (1|PTID), ~ Age, verbose = TRUE, numCores = 35 )

runZIGLMM <- function(TSAM_Object,
                      cellTypeName = NULL,
                      continuousFormula = NULL,
                      ziformula = NULL,
                      initialSampling = 5,
                      verbose = FALSE,
                      numCores = 1) {
  if (any(c(class(continuousFormula), class(ziformula)) != "formula")) {
    stop("continuousFormula and/or ziformula was not provided as a formula.")
  }

  if (is.null(cellTypeName)) {
    stop("No cell type name was provided.")
  } else if (length(cellTypeName) > 1) {
    stop("Please provide only one string within cellTypeName. If you want to run over multiple cell types, please use combineSampleTileMatrix() to generate a new object, and use that object instead, with cellTypeName = 'counts'")
  } else if (!cellTypeName %in% names(SummarizedExperiment::assays(TSAM_Object))) {
    stop("No cell type name not found within TSAM_Object.")
  }

  newObj <- combineSampleTileMatrix(subsetMOCHAObject(TSAM_Object, subsetBy = 'celltype', groupList = cellTypeName, subsetPeaks = TRUE))
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
          family = gaussian(),
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
    nullDF <- lapply(summary(modelList[[idx[1]]])$coefficients, as.data.frame)
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
    cl = cl, varlist = c("continuousFormula", "ziformula", "modelingData", "MetaDF", "individualZIGLMM", "nullDF"),
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
        extractVariable(x[[varType]], valType)
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

    colnames(zi1) <- paste("ZI_", colnames(zi1), paste = "")
    cbind(cond1, zi1)
  })

  names(combinedList) <- c("Slopes", "Significance", "StdError")

  results <- SummarizedExperiment::SummarizedExperiment(
    combinedList,
    rowRanges = SummarizedExperiment::rowRanges(TSAM_Object),
    metadata = metadata(TSAM_Object)
  )

  return(results)
}

extractVariable <- function(varDF, variable) {
  if (dim(varDF)[2] > 0) {
    val_tmp <- unlist(varDF[, variable])
    names(val_tmp) <- rownames(varDF)
    val_tmp
  } else {
    NULL
  }
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
individualZIGLMM <- function(x) {
  df <- data.frame(
    exp = as.numeric(modelingData[x, ]),
    MetaDF, stringsAsFactors = FALSE
  )

  output_vector <- tryCatch(
    {
      modelRes <- glmmTMB::glmmTMB(continuousFormula,
        ziformula = ziformula,
        data = df,
        family = stats::gaussian(),
        REML = TRUE
      )
      lapply(summary(modelRes)$coefficients, as.data.frame)
    },
    error = function(e) {
      nullDF
    }
  )
  return(output_vector)
}


## Code for testing out formulas on the data. Runs a given formula on a subset of the data, and returns the model results.
## This is meant to help during the model selection process.

# example: pilotZIGLMM(STM, 'CD16 Mono',exp~ Age + Sex + days_since_symptoms + (1|PTID), ~ 0, verbose = TRUE )

pilotZIGLMM <- function(TSAM_Object,
                        cellTypeName = NULL,
                        continuousFormula = NULL,
                        ziformula = NULL,
                        verbose = FALSE,
                        pilotIndices = 1:10) {
  if (any(c(class(continuousFormula), class(ziformula)) != "formula")) {
    stop("continuousFormula and/or ziformula was not provided as a formula.")
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
  modelingData <- modelingData[pilotIndices, match(colnames(modelingData), MetaDF$Sample)]



  # Subset metadata to just the variables needed. This minimizes overhead for parallelization
  MetaDF <- MetaDF[, colnames(MetaDF) %in% c("Sample", variableList)]

  if (verbose) {
    message("Modeling results.")
  }

  # Make your clusters for efficient parallelization

  modelList <- pbapply::pblapply(X = pilotIndices, function(x) {
    df <- data.frame(
      exp = as.numeric(modelingData[x, ]),
      MetaDF, stringsAsFactors = FALSE
    )

    glmmTMB::glmmTMB(continuousFormula,
      ziformula = ziformula,
      data = df,
      family = stats::gaussian(),
      REML = TRUE
    )
  }, cl = NULL)

  return(modelList)
}
