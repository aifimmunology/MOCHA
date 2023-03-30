## Initial sampling: When the model fails to converge, it needs to return a default data.frame with NAs that is the same size.
## The initial sampling provides a blueprint for how big that NA data.frame should be. The default is that it'll run 5 models, and if they are all NA, through an error.
## This function should only be used on continuous, non-zero inflated data.

require(lmerTest)


#' @title Run Linear Mixed-Effects Modeling for continuous,
#'  non-zero inflated data
#'
#' @description \code{runLMEM} Runs linear mixed-effects modeling for
#'   continuous, non-zero inflated data using \code{\link[lmerTest]{lmer}}
#'
#' @param ExperimentObj A SummarizedExperiment object generated from
#'   getSampleTileMatrix, chromVAR, or other.
#' @param modelFormula The formula to use with lmerTest::lmer, in the
#'   format (exp ~ factors). All factors must be found in column names
#'   of the ExperimentObj metadata.
#' @param initialSampling Size of data to use for pilot 
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @param numCores integer. Number of cores to parallelize across.
#'
#' @return
#'
#'
#'
#' @examples
#' \dontrun{
#'
#' }
#'
#' @export
runLMEM <- function(ExperimentObj,
                    modelFormula = NULL,
                    initialSampling = 5,
                    verbose = FALSE,
                    numCores = 1) {
  Sample <- NULL
  modelingData <- as.data.frame(
    SummarizedExperiment::assays(ExperimentObj)[[1]]
  )
  MetaDF <- as.data.frame(SummarizedExperiment::colData(ExperimentObj))

  if (
    !all(all.vars(modelFormula) %in% c("~", "exp", colnames(MetaDF)))
  ) {
    stop(
      "Model formula is not in the correct format ",
      "(exp ~ factors) or model factors are not found ",
      "in column names of metadata within the ExperimentObj."
    )
  }

  variableList <- all.vars(modelFormula)[all.vars(modelFormula) != "exp"]

  MetaDF <- dplyr::filter(MetaDF, Sample %in% colnames(modelingData))
  modelingData <- modelingData[, match(colnames(modelingData), MetaDF$Sample)]

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

  if (all(is.na(unlist(modelList)))) {
    stop(
      "For the initial sampling, every test model failed.",
      "Reconsider modelFormula or increase 'initialSampling'."
    )
  } else {
    idx <- which(!is.na(unlist(modelList)))
    nullDF <- as.data.frame(summary(modelList[[idx[1]]])$coefficients)
    nullDF[!is.na(nullDF)] <- NA
    rm(modelList)
  }

  cl <- parallel::makeCluster(numCores)
  parallel::clusterExport(
    cl = cl, varlist = c(
      "modelFormula", "modelingData",
      "MetaDF", "individualLMEM", "nullDF"
    ),
    envir = environment()
  )
  parallel::clusterEvalQ(cl, {
    library(lmerTest)
  })
  coeffList <- pbapply::pblapply(
    cl = cl,
    X = rownames(modelingData),
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

  output_list <- list("Slopes" = slopes,
                      "Significance" = significance,
                      "StdError" = stdError)

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
#'   on data provided. Written for efficient parallelization.
#' @param refList. A list where the first index is a data.frame 
#'   to use for modeling, and the second is the formula for modeling.
#' @return A linear model
#'
#' @noRd
individualLMEM <- function(x) {
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


## Code for testing out formulas on the data. Runs a given formula
## on a subset of the data, and returns the model results.
## This is meant to help during the model selection process.

#' @title Execute a pilot run of single linear model on a subset of data
#'
#' @description \code{pilotLMEM} Runs linear mixed-effects modeling for
#'   continuous, non-zero inflated data using \code{\link[lmerTest]{lmer}}
#'
#' @param ExperimentObj A SummarizedExperiment object generated from
#'   getSampleTileMatrix, chromVAR, or other.
#' @param modelFormula The formula to use with lmerTest::lmer, in the
#'   format (exp ~ factors). All factors must be found in column names
#'   of the ExperimentObj metadata.
#' @param pilotIndices A vector of integers defining the subset of
#'   the ExperimentObj matrix. Default is 1:10.
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @param numCores integer. Number of cores to parallelize across.
#'
#' @return
#'
#'
#'
#' @examples
#' \dontrun{
#'
#' }
#'
#' @export
pilotLMEM <- function(ExperimentObj,
                      cellPopulation = NULL,
                      modelFormula = NULL,
                      pilotIndices = 1:10,
                      verbose = FALSE) {
  if (length(cellPopulation) > 1) {
    stop("More than one cell population was provided. ",
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
      ExperimentObj, cellPopulation = cellPopulation, NAtoZero = TRUE)
  )
  MetaDF <- as.data.frame(SummarizedExperiment::colData(ExperimentObj))

  if (
    !all(all.vars(modelFormula) %in% c("~", "exp", colnames(MetaDF)))
  ) {
    stop(
      "Model formula is not in the correct format ",
      "(exp ~ factors) or model factors are not found ",
      "in column names of metadata within the ExperimentObj."
    )
  }

  variableList <- all.vars(modelFormula)[all.vars(modelFormula) != "exp"]

  MetaDF <- dplyr::filter(MetaDF, Sample %in% colnames(modelingData))
  modelingData <- modelingData[pilotIndices, match(colnames(modelingData), MetaDF$Sample)]

  # Subset metadata to just the variables in modelFormula
  MetaDF <- MetaDF[, colnames(MetaDF) %in% c("Sample", variableList)]

  modelList <- pbapply::pblapply(X = pilotIndices, function(x) {
    df <- data.frame(
      exp = as.numeric(modelingData[x, ]),
      MetaDF, stringsAsFactors = FALSE
    )

    lmerTest::lmer(formula = modelFormula, data = df)
  }, cl = NULL)

  return(modelList)
}
