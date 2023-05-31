#' @title Run Linear Mixed-Effects Modeling for continuous,
#'  non-zero inflated data
#'
#' @description \code{runLMEM} Runs linear mixed-effects modeling for
#'   continuous, non-zero inflated data using \code{\link[lmerTest]{lmer}}
#'
#' @param ExperimentObj A SummarizedExperiment object generated from
#'   getSampleTileMatrix, chromVAR, or other. It is expected to contain only
#'   one assay, or only the first assay will be used for the model.
#'   Data should not be zero-inflated.
#' @param modelFormula The formula to use with lmerTest::lmer, in the
#'   format (exp ~ factors). All factors must be found in column names
#'   of the ExperimentObj metadata. modelFormula must start with 'exp' as the response.
#'   See \link[lmerTest]{lmer}.
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
runLMEM <- function(ExperimentObj,
                    modelFormula = NULL,
                    initialSampling = 5,
                    verbose = FALSE,
                    numCores = 1) {
  Sample <- NULL
  if (!requireNamespace("lmerTest", quietly = TRUE)) {
    stop(
      "Package 'lmerTest' is required for runLMEM ",
      "Please install 'lmerTest' to proceed."
    )
  }
  
  modelingData <- as.data.frame(
    SummarizedExperiment::assays(ExperimentObj)[[1]]
  )
  MetaDF <- as.data.frame(SummarizedExperiment::colData(ExperimentObj))
  
  if (!methods::is(modelFormula, "formula")){
    stop("modelFormula is not a formula. modelFormula must be a formula in the format ",
        "(exp ~ factors)")   
  }
  
  if (!"exp" %in% all.vars(modelFormula)) {
    stop(
      "modelFormula is not in the format (exp ~ factors). ",
      "modelFormula must start with 'exp' as the response."
    )
  }

  if (
    !all(all.vars(modelFormula) %in% c("~", "exp", colnames(MetaDF)))
  ) {
    stop(
      "Model factors are not found ",
      "in the 'colData' of the ExperimentObj, or ",
      "modelFormula is not in the format ",
      "(exp ~ factors)."
    )
  }

  variableList <- all.vars(modelFormula)[all.vars(modelFormula) != "exp"]

  MetaDF <- dplyr::filter(MetaDF, Sample %in% colnames(modelingData))
  modelingData <- modelingData[
    , match(colnames(modelingData), MetaDF$Sample), drop=FALSE
  ]

  # Subset metadata to just the variables in modelFormula
  MetaDF <- MetaDF[, colnames(MetaDF) %in% c("Sample", variableList)]

  if (verbose) {
    message("Running an initial model on sample of data.")
  }

  # Generate pilot data for the null data.frame
  pilotIndices <- sample(
    x = seq_len(dim(modelingData)[1]),
    size = initialSampling,
    replace = FALSE
  )
  modelList <- pbapply::pblapply(X = pilotIndices, function(x) {
    df <- data.frame(
      exp = as.numeric(modelingData[x, ]),
      MetaDF, stringsAsFactors = FALSE
    )

    tryCatch(
      {
        lmerTest::lmer(formula = modelFormula, data = df)
      },
      error = function(e) {
        NA
      }
    )
  }, cl = NULL)
  # Initial sampling: When the model fails to converge, we
  # return a data.frame with NAs that is the same size.
  if (all(is.na(unlist(modelList)))) {
    stop(
      "For the initial sampling, every test model failed to converge. ",
      "Reconsider modelFormula or increase 'initialSampling'."
    )
  } else {
    idx <- which(!is.na(unlist(modelList)))
    nullDF <- as.data.frame(summary(modelList[[idx[1]]])$coefficients)
    nullDF[!is.na(nullDF)] <- NA 
    rm(modelList)
    # Why do we make and then export to cluster this nullDF if it is not used?
  }

  cl <- parallel::makeCluster(numCores)
  iterList <- lapply(rownames(modelingData), function(x){
    list(x, modelFormula, modelingData, MetaDF, nullDF)
  })
  parallel::clusterExport(
    cl = cl, varlist = c(
      iterList
      # "modelFormula", "modelingData",
      # "MetaDF", "individualLMEM", "nullDF" 
    ),
    envir = environment()
  )
  parallel::clusterEvalQ(cl, {
    library(lmerTest)
  })
  coeffList <- pbapply::pblapply(
    cl = cl,
    X = iterList,
    individualLMEM
  )
  parallel::stopCluster(cl)

  if (verbose) {
    message("Reorganizing coefficients.")
  }

  slopes <- do.call("rbind", pbapply::pblapply(X = coeffList, function(x) {
    slope_tmp <- x$Estimate
    names(slope_tmp) <- rownames(x)
    slope_tmp
  }, cl = NULL))

  significance <- do.call(
    "rbind", pbapply::pblapply(X = coeffList, function(x) {
      sig_tmp <- x$"Pr(>|t|)"
      names(sig_tmp) <- rownames(x)
      sig_tmp
    }, cl = NULL)
  )

  stdError <- do.call(
    "rbind", pbapply::pblapply(X = coeffList, function(x) {
      error_tmp <- x$"Std. Error"
      names(error_tmp) <- rownames(x)
      error_tmp
    }, cl = NULL)
  )

  rownames(slopes) <- rownames(modelingData)

  rownames(stdError) <- rownames(significance) <- rownames(modelingData)

  output_list <- list(
    "Slopes" = slopes,
    "Significance" = significance,
    "StdError" = stdError
  )

  results <- SummarizedExperiment::SummarizedExperiment(
    output_list,
    rowData = SummarizedExperiment::rowData(ExperimentObj),
    metadata = S4Vectors::metadata(ExperimentObj)
  )

  return(results)
}


#' @title Internal function to run linear modeling
#'
#' @description \code{IndividualLMEM} Runs linear modeling
#'   with lmerTest::lmer
#' @param iterList A list where the first index is a data.frame
#'   to use for modeling, and the second is the formula for modeling.
#'   list(x, modelFormula, modelingData, MetaDF, nullDF)
#' 
#' @return output_vector A linear model
#'
#' @noRd
individualLMEM <- function(iterList) {
  x <- iterList[[1]]
  modelFormula <- iterList[[2]]
  modelingData <- iterList[[3]]
  MetaDF <- iterList[[4]]
  nullDF <- iterList[[5]]
  
  df <- data.frame(
    exp = as.numeric(modelingData[x, ]),
    MetaDF, stringsAsFactors = FALSE
  )

  output_vector <- tryCatch(
    {
      modelRes <- lmerTest::lmer(formula = modelFormula, data = df)
      as.data.frame(summary(modelRes)$coefficients)
    },
    error = function(e) {
      nullDF
    }
  )
  return(output_vector)
}

#' @title Execute a pilot run of single linear model on a subset of data
#'
#' @description \code{pilotLMEM} Runs linear mixed-effects modeling for
#'   continuous, non-zero inflated data using \code{\link[lmerTest]{lmer}}
#'
#' @param ExperimentObj A SummarizedExperiment object generated from
#'   getSampleTileMatrix, chromVAR, or other.
#' @param cellPopulation A single cell population on which to run this pilot 
#'   model
#' @param modelFormula The formula to use with lmerTest::lmer, in the
#'   format (exp ~ factors). All factors must be found in column names
#'   of the ExperimentObj metadata.
#' @param pilotIndices A vector of integers defining the subset of
#'   the ExperimentObj matrix. Default is 1:10.
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#'
#' @return modelList a list of outputs from lmerTest::lmer
#'
#'
#' @export
pilotLMEM <- function(ExperimentObj,
                      cellPopulation = NULL,
                      modelFormula = NULL,
                      pilotIndices = 1:10,
                      verbose = FALSE) {
  Sample <- NULL
  if (!requireNamespace("lmerTest", quietly = TRUE)) {
    stop(
      "Package 'lmerTest' is required for pilotLMEM ",
      "Please install 'lmerTest' to proceed."
    )
  }
  if (length(cellPopulation) > 1) {
    stop(
      "More than one cell population was provided. ",
      "cellPopulation must be length 1. To run over multiple cell types, ",
      "run combineSampleTileMatrix() to produce ExperimentObj and set ",
      "cellPopulation = 'counts'."
    )
  } else if (
    !cellPopulation %in% names(SummarizedExperiment::assays(ExperimentObj))
  ) {
    stop("cellPopulation was not found within ExperimentObj.")
  }

  modelingData <- as.data.frame(
    MOCHA::getCellPopMatrix(
      ExperimentObj,
      cellPopulation = cellPopulation, NAtoZero = TRUE
    )
  )
  MetaDF <- as.data.frame(SummarizedExperiment::colData(ExperimentObj))
  
  if (!methods::is(modelFormula, "formula")){
    stop("modelFormula is not a formula. modelFormula must be a formula in the format ",
        "(exp ~ factors)")   
  }

  if (!"exp" %in% all.vars(modelFormula)) {
    stop(
      "modelFormula is not in the format (exp ~ factors). ",
      "modelFormula must start with 'exp' as the response."
    )
  }

  if (
    !all(all.vars(modelFormula) %in% c("~", "exp", colnames(MetaDF)))
  ) {
    stop(
      "Model factors are not found ",
      "in the 'colData' of the ExperimentObj, or ",
      "modelFormula is not in the format ",
      "(exp ~ factors)."
    )
  }


  variableList <- all.vars(modelFormula)[all.vars(modelFormula) != "exp"]

  MetaDF <- dplyr::filter(MetaDF, Sample %in% colnames(modelingData))
  modelingData <- modelingData[
    pilotIndices, match(colnames(modelingData), MetaDF$Sample), drop=FALSE
  ]

  # Subset metadata to just the variables in modelFormula
  MetaDF <- MetaDF[
    , colnames(MetaDF) %in% c("Sample", variableList), drop=FALSE
  ]
  
  modelList <- pbapply::pblapply(pilotIndices, function(x) {
    df <- data.frame(
      exp = as.numeric(modelingData[x, ]),
      MetaDF, stringsAsFactors = FALSE
    )

    lmerTest::lmer(formula = modelFormula, data = df)
  }, cl = NULL)

  return(modelList)
}
