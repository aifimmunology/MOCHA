#' @title Run Zero-inflated Generalized Linear Mixed Modeling on pseudobulked scATAC data
#'
#' @description \code{runZIGLMM} Runs linear mixed-effects modeling for
#'   zero-inflated data using \code{\link[glmmTMB]{glmmTMB}}. 
#' 
#' @param TSAM_Object A SummarizedExperiment object generated from
#'   getSampleTileMatrix. 
#' @param cellPopulation Name of a cell type(s), or 'all'. The function will combine the cell types mentioned into one matrix before running the model.
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
#' @return results a SummarizedExperiment containing ZI-GLMM results
#'
#'
#'
#' @examples
#' \dontrun{
#'   modelList <- runZIGLMM(STM[c(1:1000),], 
#'                  cellPopulation = 'CD16 Mono',
#'                  continuousFormula = exp~ Age + Sex + days_since_symptoms + (1|PTID), 
#'                  ziformula = ~ FragNumber + Age, 
#'                  verbose = TRUE, 
#'                  numCores = 35 )
#' }
#'
#' @export
#' 
runZIGLMM <- function(TSAM_Object,
                      cellPopulation = 'all',
                      continuousFormula = NULL,
                      ziformula = NULL,
                      zi_threshold = 0,
                      initialSampling = 5,
                      verbose = FALSE,
                      numCores = 2) {
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

    ## Log transform the FragmentNumbers so as to stabilize the model. But only if FragNumber is in the model. Same for CellCounts.
  if(any(colnames(MetaDF) %in% c('FragNumber'))){
    MetaDF$rawFragNumber = MetaDF$FragNumber
    MetaDF$FragNumber <- log10(MetaDF$FragNumber)
  }
  if(any(colnames(MetaDF) %in% c('CellCounts'))){
    MetaDF$rawCellCounts = MetaDF$CellCounts
    MetaDF$CellCounts <- log10(MetaDF$CellCounts)
  }


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
  NAList <- unlist(lapply(modelList, function(x){
                    ## First identify all the pilot models that ran at all without errors. 
                   tmp1 <- all(is.na(x))
                   if(tmp1){ 
                    return(FALSE)
                   }else{
                     ## Then identify all the models that converged. 
                    tryCatch(
                      {
                        tmp1 <- summary(x)$coefficients
                        return(TRUE)
                      }, error = function(e){ return(FALSE)})
                   }
            }))
   
  if (all(!NAList)) {
    stop("For the initial sampling, every test model failed. Reconsider modelFormula or increase initial sampling size.")
  } else {
    idx <- which(NAList)
    # Did any of the models have both zero-inflated and continous portions?
    bothZI_Cont <- unlist(lapply(idx, function(x){
      ziDF <- summary(modelList[[x]])$coefficients$zi
      !is.null(ziDF) | sum(dim(ziDF)) != 0

    }))
    if(all(!bothZI_Cont)){
      warning('No working models using the zero-inflated formula. Do you need to modify the zero-inflated formula?')
    }
    #Extract the first representative model. And use it to create a null template. 
    modelRes <- modelList[[idx[which(bothZI_Cont)[1]]]]
    coeff2 <- summary(modelRes)$coefficients
    coeff <- lapply(coeff2, as.data.frame)
    coeff$cond[!is.na(coeff$cond)] <- NA
    rownames(coeff$cond)[grepl('(Intercept)',rownames(coeff$cond))] = 'Intercept'
    if(all(!bothZI_Cont)){
      combinedCoeff <- coeff$cond
    }else{
      coeff$zi[!is.na(coeff$zi)] <- NA
      rownames(coeff$zi)[grepl('(Intercept)',rownames(coeff$zi))] = 'Intercept'
      rownames(coeff$zi) <- paste('ZI', rownames(coeff$zi), sep ='_')
      combinedCoeff <- rbind(coeff$cond, coeff$zi)
    }
    
    Resid <- stats::resid(modelRes)
    
    if(!all(MetaDF$Sample %in% names(Resid))){
      NA_samples <- rep(NA, sum(!MetaDF$Sample %in% names(Resid)))
      names(NA_samples) <- MetaDF$Sample[!MetaDF$Sample %in% names(Resid)]
      Resid <- c(Resid, NA_samples)
    }
    Resid <- Resid[match(names(Resid),MetaDF$Sample)]
    Resid[!is.na(Resid)] <- NA
    varCorrObj <- glmmTMB::VarCorr(modelRes)
    cond_other = unlist(varCorrObj$cond)
    names(cond_other) = paste('Cond', names(cond_other), sep = "_")
    residual = as.vector(attr(varCorrObj$cond, "sc")^2)
    names(residual) = 'Residual'
    if(!is.null(varCorrObj$zi)){
      zi_other = unlist(varCorrObj$zi)
      names(zi_other) = paste('ZI', names(zi_other), sep = "_")
      varcor_df <- c(cond_other, zi_other,residual)
    }else{
      varcor_df <- c(cond_other, residual)
    }
    varcor_df[!is.na(varcor_df)] = NA
    nullDFList <- list('Coeff' = combinedCoeff, 'Resid' = Resid, 'VCov'= varcor_df)

    rm(modelList)
  }

  if (verbose) {
    message("Modeling results.")
  }

  if(numCores <= 1){
    stop('numCores must be greater than 1. This method is meant to be parallelized.')
  }
  # Make your clusters for efficient parallelization
  cl <- parallel::makeCluster(numCores)
  parallel::clusterEvalQ(cl, {
    library(glmmTMB)
  })
  parallel::clusterExport(
    cl = cl, varlist = c("continuousFormula", "ziformula", "modelingData", "MetaDF", "individualZIGLMM", "nullDFList","zi_threshold"),
    envir = environment()
  )
  coeffList <- pbapply::pblapply(cl = cl, X = rownames(modelingData), individualZIGLMM)
  parallel::stopCluster(cl)

  
  if (verbose) {
      message("Reorganizing residuals and random effect variance.")
  }


  processedOuts <- processModelOutputs(modelOutputList = coeffList, 
                                        nullDFList = nullDFList, 
                                        rownamesList = rownames(modelingData),
                                        ranged = TRUE,
                                        SummarizedExperimentObj = newObj
                                        )
  processedOuts@metadata = append(processedOuts@metadata, list('Type' = 'scATAC'))
  return(processedOuts)
}


#' @title \code{IndividualZIGLMM}
#'
#' @description \code{IndividualZIGLMM} Runs zero-inflated linear modeling on data provided. Written for efficient parallelization.
#'
#' @param refList. A list where the first index is a data.frame to use for modeling, and the second is the formula for modeling.
#'
#' @return A linear model
#' 
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
        modelRes <- glmmTMB::glmmTMB(as.formula(continuousFormula),
          ziformula = ~ 0,
          data = df,
          family = stats::gaussian(),
          REML = TRUE,
          control = glmmTMB::glmmTMBControl(parallel = 1)
        )
       
      }else if(sum(df$exp == 0)/length(df$exp) <= zi_threshold){
        df$exp[df$exp == 0] = NA
        modelRes <- glmmTMB::glmmTMB(as.formula(continuousFormula),
          ziformula = ~ 0,
          data = df,
          family = stats::gaussian(),
          REML = TRUE,
          control = glmmTMB::glmmTMBControl(parallel = 1)
        )
       
      }else {
        modelRes <- glmmTMB::glmmTMB(as.formula(continuousFormula),
          ziformula = as.formula(ziformula),
          data = df,
          family = stats::gaussian(),
          REML = TRUE,
          control = glmmTMB::glmmTMBControl(parallel = 1)
        )
      }

      Coeff <- lapply(summary(modelRes)$coefficients, as.data.frame)
      rownames(Coeff[['cond']])[grepl('(Intercept)',rownames(Coeff[['cond']]))] = 'Intercept'
      if(sum(dim(Coeff$zi)) != 0){
        rownames(Coeff$zi)[grepl('(Intercept)',rownames(Coeff$zi))] = 'Intercept'
        rownames(Coeff$zi) <- paste('ZI', rownames(Coeff$zi), sep ='_')
      }else{
        Coeff$zi = nullDFList$Coeff[grepl('ZI',rownames(nullDFList$Coeff)),]
      }
      combinedCoeff <- rbind(Coeff$cond, Coeff$zi)

      Resid <- stats::resid(modelRes)
      
      #Clean up residuals. Data that is NA will be removed from the residuals, so we add them back in as NAs for the sake of completion. 
      if(!all(MetaDF$Sample %in% names(Resid))){
        NA_samples <- rep(NA, sum(!MetaDF$Sample %in% names(Resid)))
        names(NA_samples) <- MetaDF$Sample[!MetaDF$Sample %in% names(Resid)]
        Resid <- c(Resid, NA_samples)
      }
      Resid <- Resid[match(names(Resid),MetaDF$Sample)]

      #Now process the variance from random effects. 
      varcor_df <- tryCatch({
        varCorrObj <- glmmTMB::VarCorr(modelRes)
        cond_other = unlist(varCorrObj$cond)
        names(cond_other) = paste('Cond', names(cond_other), sep = "_")
        residual = as.vector(attr(varCorrObj$cond, "sc")^2)
        names(residual) = 'Residual'

        #Process variance
        if(!is.null(varCorrObj$zi)){
          zi_other = unlist(varCorrObj$zi)
          names(zi_other) = paste('ZI', names(zi_other), sep = "_")
          varcor_df <- c(cond_other, zi_other,residual)
        }else if(all(df$exp !=0)){
          subNull = nullDFList[grepl('ZI_', names(nullDFList))]
          zi_other = rep(0, length(subNull))
          names(zi_other) = names(subNull)
          varcor_df <- c(cond_other, zi_other,residual)
        }else {
          varcor_df <- c(cond_other, residual)
        }
        varcor_df
      }, error = function(e){
        nullDFList$VCov
      })
      
      list('Coeff' = combinedCoeff, 'Resid' = Resid, 'VCov'= varcor_df)
    },
    error = function(e) {
      nullDFList
    }
  )
  return(output_vector)
}


#' @title Run a test case for Zero-inflated Generalized Linear Mixed Modeling on pseudobulked scATAC data
#'
#' @description \code{pilotZIGLMM} Model for testing out formulas on the data. Runs a given formula on a subset of the data, and returns the model results. This is meant to help during the model selection process. \code{\link[glmmTMB]{glmmTMB}}. 
#'
#' @param TSAM_Object A SummarizedExperiment object generated from
#'   getSampleTileMatrix. 
#' @param cellPopulation Name of a cell type(s), or 'all'. The function will combine the cell types mentioned into one matrix before running the model.
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
#'   modelList <- pilotZIGLMM(STM, 
#'                            'CD16 Mono',
#'                            exp~ Age + Sex + days_since_symptoms + (1|PTID),
#'                             ~ FragNumber, verbose = TRUE )
#' }
#'
#' @export
#' 

pilotZIGLMM <- function(TSAM_Object,
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
          family = stats::gaussian(),
          REML = TRUE,
          control = glmmTMB::glmmTMBControl(parallel = numCores)
        )
       
      }else if(sum(df$exp == 0)/length(df$exp) <= zi_threshold){
        df$exp[df$exp == 0] = NA
        modelRes <- glmmTMB::glmmTMB(continuousFormula,
          ziformula = ~ 0,
          data = df,
          family = stats::gaussian(),
          REML = TRUE,
          control = glmmTMB::glmmTMBControl(parallel = numCores)
        )
       
      }else {
        modelRes <- glmmTMB::glmmTMB(continuousFormula,
          ziformula = ziformula,
          data = df,
          family = stats::gaussian(),
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


#' @title getModelValues from runZIGLMM or runLMEM output. 
#'
#' @description \code{getModelValues} Pull out a data.frame of model values for a particular row.
#' @param object A SummarizedExperiment object generated from runZIGLMM. 
#' @param rowName A string, describing the row you want to analyze. 
#'
#' @return A data.frame coefficient info by factor. 
#'
#'
#'
#' @examples
#' \dontrun{
#'   age_df <- getModelValues(runZIGLMM_output, 'Chr1:500-999')
#' }
#'
#' @export
#' 

getModelValues <- function(object, rowName){
    
    if(length(rowName) != 1){
      stop('Please provide just one string for rowName, not a list or vector.')
    }
    if(!rowName %in% rownames(object)){
      stop('rowName not found within the SummarizedObject provided.')
    }
    newDF <- do.call('rbind', lapply(as.list(SummarizedExperiment::assays(object)), function(x){
      x[rowName,, drop = FALSE]
    }))
    rownames(newDF) <- names(SummarizedExperiment::assays(object))
    return(newDF)
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

plotZIModels <- function(modelList, x_var = 'days_since_symptoms', group ='PTID',  colour ='PTID', returnMatrix = FALSE, rowNames = NULL){


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

  
#' @title plotModelPredictions
#'
#' @description \code{plotModelPredictions} Uses the raw data object and the model prediction from runZIGLMM or runLMEM to generate a data.frame for each row for plotting.
#' @param modelObject a SummarizedOutput object from runZIGLMM or runLMEM. 
#' @param dataObject The SummarizedExperiment object
#' @param specVariable The variable of interest. If other fixed effects are present, they will be adjusted for.
#' @param rowNames rownames of the specific 
#' @return a list of data.frames, including the original data, the adjusted data, and the predicted group trend. 
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
#' 
plotModelPredictions <- function(modelObject, dataObject, specVariable = 'days_since_symptoms', rowNames = NULL){

  if(grepl('SummarizedExperiment', class(modelObject)[1]) & grepl('SummarizedExperiment', class(dataObject)[1])){
      if(is.null(rowNames)){
        stop('Please provide the rownames that you want to plot via rowNames.')
      }else if(!all(rowNames %in% rownames(modelObject))){
        stop('Some rowNames not found within modelObject.')
      }else if(!all(rowNames %in% rownames(dataObject))){
        stop('Some rowNames not found within dataObject.')
      }else if(! specVariable %in% names(SummarizedExperiment::assays(modelObject))){
        stop('specVariable does not appear to be a fixed effect within the modelObject.',
                ' Please doublecheck the specVariable. It should be the name of an assay within the modelObject, which comes from runLMEM or runZIGLMM.')
      }else if(! specVariable %in% names(SummarizedExperiment::assays(modelObject))){
        stop('specVariable does not appear to be a fixed effect within the modelObject.',
                ' Please doublecheck the specVariable. It should be the name of an assay within the modelObject, which comes from runLMEM or runZIGLMM.')
      }

      allModels <- do.call('rbind', lapply(rowNames, function(x){
        tmpValues <- getModelValues(modelObject, x)
        t(tmpValues[,'Estimate', drop= FALSE])
      }))
      rownames(allModels) <- rowNames

      metaData
      
      adjustedDFs <- lapply(allModels, function(x){
           if(all(is.na(x))){
              return(NA)
            }else{
              df <- x
              
              numericCoefs <- names(coefs)[names(coefs) %in% colnames(df) & names(coefs) != x_var]

              df$Prediction = coefs[1] + sum(unlist(lapply(numericCoefs, function(z) { coefs[z]* mean(df[,z])}))) + coefs[x_var]*df[,x_var]

              p1 <- ggplot(df, aes_string(x=x_var,y='exp', colour=colour, group = group)) + geom_point(size=0.8, alpha=0.4)+
                  geom_line(data=df[df$exp>0,],aes_string(x=x_var, y='exp',alpha=0.3), linewidth=0.5)+
                  ggtitle(paste(tile, 'Suceeded',sep = ' '))+ 
                  geom_line(data=df,aes_string(x=x_var,y='Prediction'),linewidth=1, col='black') + theme_minimal() + 
                  theme(legend.position = 'none')
              return(p1)
            }

        }) 

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

