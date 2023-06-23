
#' @title pilot_scATAC
#'
#' @description \code{pilot_s} Model for testing out ZI-GLMM formulas on the data from MOCHA. Runs a given formula on a subset of the data, and returns the model results. This is meant to help during the model selection process. \code{\link[glmmTMB]{glmmTMB}}. 
#'
#' @param TSAM_Object A SummarizedExperiment object generated from
#'   getSampleTileMatrix. 
#' @param cellPopulation Name of a cell type within the TSAM_Object
#' @param continuousFormula The formula for the continuous data that should be used within glmmTMB. It should be in the
#'   format (exp ~ factors). All factors must be found in column names
#'   of the TSAM_Object metadata, except for CellType, FragNumber and CellCount, which will be extracted from the TSAM_Object.
#'   modelFormula must start with 'exp' as the response.
#'   See \link[glmmTMB]{glmmTMB}.
#' @param ziformula The formula for the zero-inflated data that should be used within glmmTMB. It should be in the
#'   format ( ~ factors). All factors must be found in column names
#'   of the TSAM_Object colData metadata, except for CellType, FragNumber and CellCount, which will be extracted from the TSAM_Object.
#'   FragNumber and CellCounts will be log10 normalized within the function. 
#' @param zi_threshold Zero-inflated threshold ( range = 0-1), representing the fraction of samples with zeros. At or above this threshold, the zero-inflated modeling kicks in.
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @param pilotIndices integer. Specific locations to test within the peakset of the cell type chosen. 
#'
#' @return results model results
#'
#'
#'
#' @examples
#' \dontrun{
#'   modelList <- pilot_TilesGeneLinks(STM, 
#'                            'CD16 Mono',
#'                            exp~ Age + Sex + days_since_symptoms + (1|PTID),
#'                             ~ FragNumber, verbose = TRUE )
#' }
#'
#' @export
#' 

pilot_TilesGeneLinks <- function(TSAM_Object,
                        cellPopulation = NULL,
                        continuousFormula = NULL,
                        ziformula = NULL,
                        zi_threshold = 0,
                        verbose = FALSE,
                        numCores = 1,
                        pilotIndices = 1:10) {
  if (any(c(class(continuousFormula), class(ziformula)) != "formula")) {
    stop("continuousFormula and/or ziformula was not provided as a formula.")
  }

  if (zi_threshold < 0 | zi_threshold > 1 | ! is.numeric(zi_threshold)) {
    stop("zi_threshold must be between 0 and 1.")
  }

 if (length(cellPopulation) > 1) {
    stop(
      "More than one cell population was provided. ",
      "cellPopulation must be length 1. To run over multiple cell types, ",
      "run combineSampleTileMatrix() to produce a new combined TSAM_Object and set ",
      "cellPopulation = 'counts'."
    )
  } else if (
    (!cellPopulation %in% names(SummarizedExperiment::assays(TSAM_Object)))
  ) {
    stop("cellPopulation was not found within TSAM_Object.")
  } else if(cellPopulation == 'counts'){
    newObj <- TSAM_Object
  }else{
    newObj <- combineSampleTileMatrix(subsetMOCHAObject(TSAM_Object, subsetBy = 'celltype', groupList = cellPopulation, subsetPeaks = TRUE))
  }

  modelList <- .pilotModels_generic(SE_Object = newObj,
                        continuousFormula = continuousFormula,
                        ziformula = ziformula,
                        zi_threshold = zi_threshold,
                        verbose = verbose,
                        numCores = numCores,
                        pilotIndices = pilotIndices)
  
  return(modelList)
}


#' @title pilot_scRNA
#'
#' @description \code{pilot_scRNA} Model for testing out GLM formulas on the pseudobulked RNA via ChAI.
#'   Runs a given formula on a subset of the data, and returns the model results. This is meant to help during the model selection process. \code{\link[glmmTMB]{glmmTMB}}. 
#'
#' @param scRNA_Object A SummarizedExperiment object generated from
#'   makePseudobulkRNA. 
#' @param cellPopulation Name of a cell type within the scRNA_Object
#' @param continuousFormula The formula for the continuous data that should be used within glmmTMB. It should be in the
#'   format (exp ~ factors). All factors must be found in column names
#'   of the scRNA_Object metadata, except for CellType, FragNumber and CellCount, which will be extracted from the scRNA_Object.
#'   modelFormula must start with 'exp' as the response.
#'   See \link[glmmTMB]{glmmTMB}.
#' @param family String. Can be 'negativeBinomial1', 'negativeBinomial2', or 'poisson'. Default is 'poisson', which handles the zeros best. 
#' @param zi_threshold Zero-inflated threshold ( range = 0-1), representing the fraction of samples with zeros. At or above this threshold, the zero-inflated modeling kicks in.
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @param pilotIndices integer. Specific locations to test within the peakset of the cell type chosen. 
#'
#' @return results model results
#'
#'
#'
#' @examples
#' \dontrun{
#'   modelList <- pilot_scRNA(scRNA_Object, 
#'                            'CD16 Mono',
#'                            exp~ Age + Sex + days_since_symptoms + (1|PTID),
#'                              verbose = TRUE )
#' }
#'
#' @export
#' 

pilot_scRNA <- function(scRNA_Object,
                        cellPopulation = NULL,
                        continuousFormula = NULL,
                        verbose = FALSE,
                        family = 'poisson',
                        numCores = 1,
                        pilotIndices = 1:10) {
  if (any(c(class(continuousFormula), class(ziformula)) != "formula")) {
    stop("continuousFormula and/or ziformula was not provided as a formula.")
  }

  if (zi_threshold < 0 | zi_threshold > 1 | ! is.numeric(zi_threshold)) {
    stop("zi_threshold must be between 0 and 1.")
  }

 if (length(cellPopulation) > 1) {
    stop(
      "More than one cell population was provided. ",
      "cellPopulation must be length 1. To run over multiple cell types, ",
      "run combinePseudobulkRNA() to produce a new combined scRNA_Object and set ",
      "cellPopulation = 'counts'."
    )
  } else if (
    (!cellPopulation %in% names(SummarizedExperiment::assays(scRNA_Object)))
  ) {
    stop("cellPopulation was not found within scRNA_Object.")
  } else if(cellPopulation == 'counts'){
    newObj <- scRNA_Object
  }else{
    #newObj <- combinePseudobulkRNA(scRNA_Object)
    #Needs to be generated, analoguous to combineSampleTileMatrix
    newObj <- scRNA_Object
    SummarizedExperiment::assays(newObj) <- SummarizedExperiment::assays(scRNA_Object)['CD16_Mono']
    names(SummarizedExperiment::assays(newObj) ) <- 'counts'
  }

  if(tolower(family) == 'negativebinomial2'){ 
    family = glmmTMB::nbinom2()
  }else if(tolower(family)== 'negativebinomial1'){
    family = glmmTMB::nbinom1()
  }else if(tolower(family) == 'poisson'){
    family = stats::poisson()
  }else{
    stop('family not recognized.')
  }

  modelList <- .pilotModels_generic(SE_Object = newObj,
                        continuousFormula = continuousFormula,
                        ziformula = ~0,
                        family = family,
                        zi_threshold = 0,
                        verbose = verbose,
                        numCores = numCores,
                        pilotIndices = pilotIndices)
  
  return(modelList)
}






#' @title Execute a pilot run of single linear model on a subset of data
#'
#' @description \code{pilot_General} Runs linear mixed-effects modeling for
#'   continuous, non-zero inflated data using \code{\link[lmerTest]{lmer}}
#'
#' @param ExperimentObj A SummarizedExperiment-type object generated from
#'   makeChromVAR, makePseudobulkRNA, or other. Objects from getSampleTileMatrix can work, 
#'   but we recommend runZIGLMM for those objects, not runLMEM>
#' @param assayName a character string, matching the name of an assay within the SummarizedExperiment. 
#'   The assay named will be used for modeling. 
#' @param modelFormula The formula to use with lmerTest::lmer, in the
#'   format (exp ~ factors). All factors must be found in column names
#'   of the ExperimentObj metadata.
#' @param family distribution family parameter, passed to glmmTMB to describe the data's distribution.
#'     Default is normal (gaussian()). See  \link[glmmTMB]{glmmTMB}.
#' @param pilotIndices A vector of integers defining the subset of
#'   the ExperimentObj matrix. Default is 1:10.
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @param numCores integer. Number of cores to parallelize across.
#'
#' @return modelList a list of outputs from lmerTest::lmer
#'
#'
#' @export
pilot_General <- function(ExperimentObj,
                      assayName = NULL,
                      modelFormula = NULL,
                      numCores = 1,
                      family = stats::gaussian(),
                      pilotIndices = 1:10,
                      verbose = FALSE) {
  if (length(assayName) > 1) {
    stop(
      "More than one assay was provided. ",
      "assayName must be length 1. To run over multiple assays/cell types, ",
      "combine the matrices into one and wrap them in a SummarizedExperiment. Then set assayName to ",
      "the name of that matrix in the obejct."
    )
  } else if (
    !assayName %in% names(SummarizedExperiment::assays(ExperimentObj))
  ) {
    stop("assayName was not found within ExperimentObj.")
  }

  names(SummarizedExperiment::assays(ExperimentObj))[names(SummarizedExperiment::assays(ExperimentObj)) == assayName] = 'counts'

  modelList <- .pilotModels_generic(SE_Object = ExperimentObj,
                        continuousFormula = modelFormula,
                        ziformula = ~0,
                        zi_threshold = 0,
                         family = family,
                        verbose = verbose,
                        numCores = numCores,
                        pilotIndices = pilotIndices)
  return(modelList)
  
}



#' @title pilotModels_generic 
#'
#' @description \code{mpilotModels_generic} Internal generic funciton for modeling. 
#' @param SE_Object A SummarizedExperiment object generated from
#' @param continuousFormula The formula for the continuous data that should be used within glmmTMB. It should be in the
#'   format (exp ~ factors). All factors must be found in column names of the SE_Object colData (i.e. sample metadata).
#'   See \link[glmmTMB]{glmmTMB}.
#' @param ziformula The formula for the zero-inflated data that should be used within glmmTMB. It should be in the
#'   format ( ~ factors). All factors must be found in column names
#'   of SE_Object colData (i.e. sample metadata).
#' @param family distribution family parameter, passed to glmmTMB to describe the data's distribution.
#'     Default is normal (gaussian()). See  \link[glmmTMB]{glmmTMB}.
#' @param zi_threshold Zero-inflated threshold ( range = 0-1), representing the fraction of samples with zeros.
#'           At or above this threshold, the zero-inflated modeling kicks in.
#' @param pilotIndices integer. Specific locations to test within the peakset of the cell type chosen. 
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @param numCores integer. Number of cores to parallelize across.
#'
#' @return a list of model objects from glmmTMB
#'
#' @noRd
#' 
.pilotModels_generic <- function(SE_Object,
                        continuousFormula = NULL,
                        ziformula = NULL,
                        family = stats::gaussian(),
                        zi_threshold = 0,
                        verbose = FALSE,
                        numCores = 1,
                        pilotIndices = 1:10) {

  modelingData <- log2(SummarizedExperiment::assays(SE_Object)[['counts']]+1)
  MetaDF <- as.data.frame(SummarizedExperiment::colData(SE_Object))

  if (!all(all.vars(continuousFormula) %in% c("exp", colnames(MetaDF)))) {
    stop("Model formula is not in the correct format (exp ~ factors) or model factors are not found in sample-level metadata.")
  }

  if (!all(all.vars(ziformula) %in% c(colnames(MetaDF))) & length(all.vars(ziformula)) > 0) {
    stop("factors from the ziformula were not found in the metadata.")
  }

  variableList <- c(all.vars(continuousFormula)[all.vars(continuousFormula) != "exp"], all.vars(ziformula))

  MetaDF <- dplyr::filter(MetaDF, Sample %in% colnames(modelingData))
  pilotNames <- rownames(modelingData)[pilotIndices]

  ## Log transform the FragmentNumbers so as to stabilize the model. But only if FragNumber is in the model. Same for CellCounts.
  if(any(colnames(MetaDF) %in% c('FragNumber'))){
    MetaDF$rawFragNumber = MetaDF$FragNumber
    MetaDF$FragNumber = log10(MetaDF$FragNumber)
  }
  if(any(colnames(MetaDF) %in% c('CellCounts'))){
    MetaDF$rawCellCounts = MetaDF$CellCounts
    MetaDF$CellCounts = log10(MetaDF$CellCounts)
  }

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
    tryCatch({
    if(sum(df$exp == 0)/length(df$exp) == 0){
        modelRes <- glmmTMB::glmmTMB(continuousFormula,
          ziformula = ~ 0,
          data = df,
          family = family,
          REML = TRUE,
          control = glmmTMB::glmmTMBControl(parallel = numCores)
        )
       
      }else if(sum(df$exp == 0)/length(df$exp) <= zi_threshold){
        df$exp[df$exp == 0] = NA
        modelRes <- glmmTMB::glmmTMB(continuousFormula,
          ziformula = ~ 0,
          data = df,
          family = family,
          REML = TRUE,
          control = glmmTMB::glmmTMBControl(parallel = numCores)
        )
       
      }else {
        modelRes <- glmmTMB::glmmTMB(continuousFormula,
          ziformula = ziformula,
          data = df,
          family = family,
          REML = TRUE,
          control = glmmTMB::glmmTMBControl(parallel = numCores)
        )
      }
    }, error = function(e){
      list(e,x, df)
    })

  }, cl = NULL)
  names(modelList) <- pilotNames

  return(modelList)
}



#' @title plotZIModels
#'
#' @description \code{plotZIModels} plots the data and the generalized slope from the output of pilotZIGLMM
#' @param modelList a list of models, output from pilotZIGLMM or pilotLMEM
#' @param x The variable that you want to plot the x-axis over. 
#' @param group the group column name, for ggplot line plots.
#' @param colour the colour parameter for ggplot
#' @param returnMatrix Returns the data.frame of values and predictions, without plotting. 
#' 
#' @return list of gpplots, or a list of data.frames. 
#'
#'
#'
#' @examples
#' \dontrun{
#'   predictedDF <- plotZIModels(modelList,returnMatrix = TRUE)
#' }
#'
#' @export
#' 

plotModels <- function(modelList, x_var = 'days_since_symptoms', group ='PTID',  colour ='PTID', returnMatrix = FALSE, rowNames = NULL){


    if(any(grepl('list|List',class(modelList)[1]))){

      pList <- lapply(seq_along(modelList), function(i){
          
          if(class(modelList[[i]]) == 'glmmTMB'){
            #Extract the original values
            df <- as.data.frame(modelList[[i]]$frame)
            tile <- names(modelList)[i]
            tryCatch({
              #Extract model coefficients
              sum_fit = summary(modelList[[i]])
              coefs = sum_fit$coefficients$cond[,1]
              numericCoefs <- names(coefs)[names(coefs) %in% colnames(df) & names(coefs) != x_var]

              df$Prediction = coefs[1] + sum(unlist(lapply(numericCoefs, function(z) { coefs[z]* mean(df[,z])}))) + coefs[x_var]*df[,x_var]

              if(returnMatrix){
                  
                  df
                  
              }else{
                  ggplot(df, aes_string(x=x_var,y='exp', colour=colour, group = group)) + geom_point(size=0.8, alpha=0.4)+
                      geom_line(data=df[df$exp>0,],aes_string(x=x_var, y='exp',alpha=0.3), linewidth=0.5)+
                      ggtitle(paste(tile, 'Suceeded',sep = ' '))+ 
                      geom_line(data=df,aes_string(x=x_var,y='Prediction'),linewidth=1, col='black') + theme_minimal() + 
                      theme(legend.position = 'none')
                      
                }
            },error=function(e){
                df <- as.data.frame(modelList[[i]]$frame)
                tile <- names(modelList)[i]
                if(returnMatrix){
                    
                    df
                    
                }else{
                    ggplot(df, aes_string(x=x_var,y='exp', colour=colour, group = group)) + geom_point(size=0.8, alpha=0.4)+
                        geom_line(data=df[df$exp>0,],aes_string(x=x_var, y='exp',alpha=0.3), linewidth=0.5)+
                        ggtitle(paste(tile,e,sep = ' '))+ 
                        theme_minimal() + 
                        theme(legend.position = 'none')
                        
                }
            })

          }else{
            df <- as.data.frame(modelList[[i]][[3]])
            tile = modelList[[i]][[2]]
            if(returnMatrix){
                    
                    df
                    
            }else{

              ggplot(df, aes_string(x=x_var,y='exp', colour=colour, group = group)) + geom_point(size=0.8, alpha=0.4)+
                        geom_line(data=df[df$exp>0,],aes_string(x=x_var, y='exp',alpha=0.3), linewidth=0.5)+
                        ggtitle(paste(tile, modelList[[i]][[1]],sep = ' '))+ 
                        theme_minimal() + 
                        theme(legend.position = 'none')
            }
          }
          
      })
      return(pList)
    }else{
      stop('modelList type not recognized.')
    }
  }



#' @title getPilotCoefficients
#'
#' @description \code{getPilotCoefficients} Attempts to pull coefficients from a list of models from pilotZIGLMM or pilotLMEM. 
#'    Returns a list of either the coefficients for each
#'    model, or the error generated when attempting to get coefficients. 
#' @param  pilotModelList 
#' @return modelList a list of outputs from lmerTest::lmer
#'
#'
#' @export

getPilotCoefficients <- function(pilotModelList){

   coeffList <-  lapply(pilotModelList, function(x) {
                    tryCatch({
                        summary(x)$coefficients
                      },
                      error=function(e){
                        list(e, x$frame)
                      })
                  })
   return(coeffList)
}
