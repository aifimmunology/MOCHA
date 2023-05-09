#' @title Run Zero-Inflated Generalized Linear Mixed-Modeling for zero inflated
#'   data
#'
#' @description \code{runZIGLMM} Runs linear mixed-effects modeling for
#'   zero inflated data using \code{\link[glmmTMB]{glmmTMB}}
#'
#' @param TSAM_Object A SummarizedExperiment object generated from
#'   getSampleTileMatrix, chromVAR, or other. It is expected to contain only one
#'   assay, or only the first assay will be used for the model. Data should not
#'   be zero-inflated.
#' @param cellPopulation A single cell population on which to run this model
#' @param continuousFormula The formula, see \code{\link[glmmTMB]{glmmTMB}}.
#'   Combined fixed and random effects formula, following lme4 syntax.
#' @param ziformula The zero-inflated formula, see
#'   \code{\link[glmmTMB]{glmmTMB}}. a one-sided (i.e., no response variable)
#'   formula for zero-inflation combining fixed and random effects: the default
#'   ~0 specifies no zero-inflation. Specifying ~. sets the zero-inflation
#'   formula identical to the right-hand side of formula (i.e., the conditional
#'   effects formula); terms can also be added or subtracted. When using ~. as
#'   the zero-inflation formula in models where the conditional effects formula
#'   contains an offset term, the offset term will automatically be dropped. The
#'   zero-inflation model uses a logit link.
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
#'   modelList <- runLMEM(ExperimentObj,
#'     modelFormula = NULL,
#'     initialSampling = 5,
#'     verbose = FALSE,
#'     numCores = 1
#'  )
#' }
#'
#' @export
runZIGLMM <- function(TSAM_Object,
                      cellPopulation = NULL,
                      continuousFormula = NULL,
                      ziformula = NULL,
                      initialSampling = 5,
                      verbose = FALSE,
                      numCores = 1) {
  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    stop(
      "Package 'glmmTMB' is required for runZIGLMM. ",
      "Please install 'glmmTMB' to proceed."
    )
  }
  Sample <- NULL
  
  if (any(c(class(continuousFormula), class(ziformula)) != "formula")) {
    stop("continuousFormula and/or ziformula was not provided as a formula.")
  }

  if (is.null(cellPopulation)) {
    stop("No cell type name was provided.")
  } else if (length(cellPopulation) > 1) {
    stop("Please provide only one string within cellPopulation. If you want to run over multiple cell types, please use combineSampleTileMatrix() to generate a new object, and use that object instead, with cellPopulation = 'counts'")
  } else if (!cellPopulation %in% names(SummarizedExperiment::assays(TSAM_Object))) {
    stop("No cell type name not found within TSAM_Object.")
  }

  newObj <- combineSampleTileMatrix(subsetMOCHAObject(TSAM_Object, subsetBy = 'celltype', groupList = cellPopulation, subsetPeaks = TRUE))
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
  iterList <- lapply(rownames(modelingData), function(x){
    list(x, continuousFormula, ziformula, modelingData, MetaDF, nullDF)
  })
  parallel::clusterEvalQ(cl, {
    library(glmmTMB)
  })
  parallel::clusterExport(
    cl = cl, varlist = c(
      iterList
      # "continuousFormula", 
      # "ziformula", 
      # "modelingData", 
      # "MetaDF", "individualZIGLMM", "nullDF"
    ),
    envir = environment()
  )
  coeffList <- pbapply::pblapply(
    cl = cl, 
    X = iterList, 
    individualZIGLMM
  )
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

extractVariable <- function(varDF, variable) {
  if (dim(varDF)[2] > 0) {
    val_tmp <- unlist(varDF[, variable])
    names(val_tmp) <- rownames(varDF)
    val_tmp
  } else {
    NULL
  }
}

#' @title Internal function to run linear modeling
#'
#' @description \code{IndividualLMEM} Runs linear modeling
#'   with glmmTMB::glmmTMB
#' @param iterList A list where the first index is a data.frame
#'   to use for modeling, and the second is the formula for modeling.
#'   list(x, continuousFormula, ziformula, modelingData, MetaDF, nullDF)
#' 
#' @return output_vector A linear model
#'
#' @noRd
individualZIGLMM <- function(iterList) {
  x <- iterList[[1]]
  continuousFormula <- iterList[[2]]
  ziformula <- iterList[[3]]
  modelingData <- iterList[[4]]
  MetaDF <- iterList[[5]]
  nullDF <- iterList[[6]]
  
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


#' @title Execute a pilot run of model on a subset of data
#'
#' @description \code{pilotLMEM} Runs linear mixed-effects modeling for
#'   zero inflated data using \code{\link[glmmTMB]{glmmTMB}}
#'
#' @param TSAM_Object A SummarizedExperiment object generated from
#'   getSampleTileMatrix, chromVAR, or other.
#' @param cellPopulation A single cell population on which to run this pilot 
#'   model
#' @param continuousFormula The formula, see \code{\link[glmmTMB]{glmmTMB}}.
#'   Combined fixed and random effects formula, following lme4 syntax.
#' @param ziformula The zero-inflated formula, see
#'   \code{\link[glmmTMB]{glmmTMB}}. a one-sided (i.e., no response variable)
#'   formula for zero-inflation combining fixed and random effects: the default
#'   ~0 specifies no zero-inflation. Specifying ~. sets the zero-inflation
#'   formula identical to the right-hand side of formula (i.e., the conditional
#'   effects formula); terms can also be added or subtracted. When using ~. as
#'   the zero-inflation formula in models where the conditional effects formula
#'   contains an offset term, the offset term will automatically be dropped. The
#'   zero-inflation model uses a logit link.
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @param pilotIndices A vector of integers defining the subset of
#'   the ExperimentObj matrix. Default is 1:10.
#'
#' @return modelList a list of outputs from glmmTMB::glmmTMB
#'
#'
#' @export
pilotZIGLMM <- function(TSAM_Object,
                        cellPopulation = NULL,
                        continuousFormula = NULL,
                        ziformula = NULL,
                        verbose = FALSE,
                        pilotIndices = 1:10) {
  Sample <- NULL
  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    stop(
      "Package 'glmmTMB' is required for pilotZIGLMM. ",
      "Please install 'glmmTMB' to proceed."
    )
  }
  if (any(c(class(continuousFormula), class(ziformula)) != "formula")) {
    stop("continuousFormula and/or ziformula was not provided as a formula.")
  }

  if (is.null(cellPopulation)) {
    stop("No cell type name was provided.")
  } else if (length(cellPopulation) > 1) {
    stop("Please provide only one string within cellPopulation. If you want to run over multiple cell types, please use combineSampleTileMatrix() to generate a new object, and use that object instead, with cellPopulation = 'counts'")
  } else if (!cellPopulation %in% names(SummarizedExperiment::assays(TSAM_Object))) {
    stop("Cell type name not found within TSAM_Object.")
  }

  newObj <- combineSampleTileMatrix(subsetMOCHAObject(TSAM_Object, subsetBy = 'celltype', groupList = cellPopulation, subsetPeaks = TRUE))
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
