#' @title Run Zero-inflated Generalized Linear Mixed Modeling on pseudobulked scATAC data
#'
#' @description \code{model_scATAC} Runs linear mixed-effects modeling for
#'   zero-inflated data using \code{\link[glmmTMB]{glmmTMB}}. 
#' 
#' @param TSAM_Object A SummarizedExperiment object generated from
#'   getSampleTileMatrix. 
#' @param cellPopulation Name of a cell type. 
#' @param continuousFormula The formula for the continuous data that should be used within glmmTMB. It should be in the
#'   format (exp ~ factors). All factors must be found in column names
#'   of the TSAM_Object metadata, except for CellType, FragNumber and CellCount, which will be extracted from the TSAM_Object.
#'   modelFormula must start with 'exp' as the response.
#'   See \link[glmmTMB]{glmmTMB}.
#' @param ziFormula The formula for the zero-inflated data that should be used within glmmTMB. It should be in the
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
#'   modelList <- model_scATAC(STM[c(1:1000),], 
#'                  cellPopulation = 'CD16 Mono',
#'                  continuousFormula = exp~ Age + Sex + days_since_symptoms + (1|PTID), 
#'                  ziFormula = ~ FragNumber + Age, 
#'                  verbose = TRUE, 
#'                  numCores = 35 )
#' }
#'
#' @export
#' 
model_scATAC <- function(TSAM_Object,
                      cellPopulation,
                      continuousFormula = NULL,
                      ziFormula = NULL,
                      zi_threshold = 0,
                      initialSampling = 5,
                      verbose = FALSE,
                      numCores = 2) {

  if (class(continuousFormula) == 'character') {
    continuousFormula <- as.formula(continuousFormula)
  }

  
  if (class(ziFormula) == 'character') {
    ziFormula <- as.formula(ziFormula)
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

  exp <- .model_pseudobulk_default(SE_Object = newObj,
                      continuousFormula = continuousFormula,
                      ziFormula = ziFormula,
                      zi_threshold = zi_threshold,
                      initialSampling = initialSampling,
                      family = stats::gaussian(),
                      modality = 'scATAC',
                      verbose = verbose,
                      numCores = numCores)
  return(exp)

}

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
#' @param assayName The name of the assay to model within the SummarizedExperiment. 
#' @param modelFormula The formula to use with lmerTest::lmer, in the
#'   format (exp ~ factors). All factors must be found in column names
#'   of the ExperimentObj metadata. modelFormula must start with 'exp' as the response.
#'   See \link[lmerTest]{lmer}.
#' @param initialSampling Size of data to use for pilot
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @param numCores integer. Number of cores to parallelize across.
#'
#' @return results a SummarizedExperiment containing LMEM results. Assays are metrics related to the model coefficients,
#'          including the Estimate, Std_Error, df, t_value, p_value. Within each assay, each row corresponds to each row of
#'          the SummarizedExperiment and columns correspond to each fixed effect variable within the model.
#'          Any row metadata from the ExperimentObject (see rowData(ExperimentObj)) is preserved in the output. 
#'          The Residual matrix and the variance of the random effects are saved in the metadata slot of the output. 
#'
#'
#'
#' @examples
#' \dontrun{
#'   modelList <- runLMEM(ExperimentObj,
#'    assayName = 'z'
#'     modelFormula = NULL,
#'     initialSampling = 5,
#'     verbose = FALSE,
#'     numCores = 1
#'  )
#' }
#'
#' @export
runLMEM <- function(ExperimentObj,
                    assayName = NULL,
                    modelFormula = NULL,
                    initialSampling = 5,
                    verbose = FALSE,
                    numCores = 2) {

  

  if (verbose) {
    message("Reorganizing coefficients.")
  }

  processedOuts <- processModelOutputs(modelOutputList = coeffList, 
                                        nullDFList = nullDFList, 
                                        rownamesList = rownames(modelingData),
                                        SummarizedExperimentObj = ExperimentObj
                                        )

  return(processedOuts)
}


#' @title Run Linear Mixed-Effects Modeling for continuous,
#'  non-zero inflated data
#'
#' @description \code{model_General} Runs linear mixed-effects modeling for
#'   continuous, non-zero inflated data using \code{\link[glmmTMB]{glmmTMB}}
#'
#' @param ExperimentObj A SummarizedExperiment object generated from chromVAR, or other.
#'   It is expected to contain normally distributed data without zero-inflation. 
#' @param assayName The name of the assay to model within the SummarizedExperiment. 
#' @param modelFormula The formula to use, in the
#'   format (exp ~ factors). All factors must be found in column names
#'   of the ExperimentObj metadata. modelFormula must start with 'exp' as the response.
#'   See \link[glmmTMB]{glmmTMB}.
#' @param family distribution family parameter, passed to glmmTMB to describe the data's distribution.
#'     Default is normal (gaussian()). See  \link[glmmTMB]{glmmTMB}.
#' @param initialSampling Size of data to use for pilot
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @param numCores integer. Number of cores to parallelize across.
#'
#' @return results a SummarizedExperiment containing LMEM results. Assays are metrics related to the model coefficients,
#'          including the Estimate, Std_Error, df, t_value, p_value. Within each assay, each row corresponds to each row of
#'          the SummarizedExperiment and columns correspond to each fixed effect variable within the model.
#'          Any row metadata from the ExperimentObject (see rowData(ExperimentObj)) is preserved in the output. 
#'          The Residual matrix and the variance of the random effects are saved in the metadata slot of the output. 
#'
#'
#'
#' @examples
#' \dontrun{
#'   modelList <- model_General(ExperimentObj,
#'    assayName = 'z'
#'     modelFormula = NULL,
#'     initialSampling = 5,
#'     verbose = FALSE,
#'     numCores = 1
#'  )
#' }
#'
#' @export
model_General <- function(ExperimentObj,
                    assayName = NULL,
                    modelFormula = NULL,
                    family = stats::gaussian(),
                    initialSampling = 5,
                    verbose = FALSE,
                    numCores = 2) {
  
  if(!any(names(SummarizedExperiment::assays(ExperimentObj)) %in% assayName)){
    stop('ExperimentObj does not contain an assay that matches the assayName input variable.')
  }
  SummarizedExperiment::assays(ExperimentObj) = SummarizedExperiment::assays(ExperimentObj)[assayName]

  exp <- .model_pseudobulk_default(SE_Object = ExperimentObj,
                      continuousFormula = as.formula(modelFormula),
                      ziFormula = ~0,
                      zi_threshold = 0,
                      family = family,
                      initialSampling = initialSampling,
                      modality = 'General',
                      verbose = verbose,
                      numCores = numCores)
  return(exp)

}




#' @title model_pseudobulk_default 
#'
#' @description \code{model_pseudobulk_default} Internal generic funciton for modeling. 
#' @param SE_Object A SummarizedExperiment object generated from
#' @param continuousFormula The formula for the continuous data that should be used within glmmTMB. It should be in the
#'   format (exp ~ factors). All factors must be found in column names of the SE_Object colData (i.e. sample metadata).
#'   See \link[glmmTMB]{glmmTMB}.
#' @param ziFormula The formula for the zero-inflated data that should be used within glmmTMB. It should be in the
#'   format ( ~ factors). All factors must be found in column names
#'   of SE_Object colData (i.e. sample metadata).
#' @param zi_threshold Zero-inflated threshold ( range = 0-1), representing the fraction of samples with zeros.
#'           At or above this threshold, the zero-inflated modeling kicks in.
#' @param initialSampling Size of data to use for pilot
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @param numCores integer. Number of cores to parallelize across.
#'
#' @return results a SummarizedExperiment containing ZI-GLMM results
#'
#' @noRd
#' 
.model_pseudobulk_default <- function(SE_Object,
                      continuousFormula = NULL,
                      ziFormula = NULL,
                      zi_threshold = 0,
                      initialSampling = 5,
                      family = stats::gaussian(),
                      modality = 'General',
                      verbose = FALSE,
                      numCores = 2) {

  if (c(class(continuousFormula)) != "formula") {
    stop("continuousFormula was not provided as a formula.")
  }
  if (c(class(ziFormula)) != "formula" & !is.null(ziFormula) ) {
    stop("ziFormula was not provided as a formula.")
  }

  if(modality == 'scATAC'){
    modelingData <- log2(SummarizedExperiment::assays(SE_Object)[[1]]+1)
  }else{
    modelingData <- SummarizedExperiment::assays(SE_Object)[[1]]
  }

  MetaDF <- as.data.frame(SummarizedExperiment::colData(SE_Object))
  if (!all(all.vars(continuousFormula) %in% c("exp", colnames(MetaDF)))) {
    stop("Model formula is not in the correct format (exp ~ factors) or model factors are not found in column names of metadata.")
  }

  if (!all(all.vars(ziFormula) %in% colnames(MetaDF))) {
    stop("factors from the ziFormula were not found in the metadata.")
  }

  variableList <- c(all.vars(continuousFormula)[all.vars(continuousFormula) != "exp"], all.vars(ziFormula))

  MetaDF <- dplyr::filter(MetaDF, Sample %in% colnames(modelingData))
  modelingData <- modelingData[, match(colnames(modelingData), MetaDF$Sample)]

  # Subset metadata to just the variables needed. This minimizes overhead for parallelization
  MetaDF <- MetaDF[, colnames(MetaDF) %in% c("Sample", variableList)]

  #Transform in character strings for multithreading. 
  continuousFormula <- deparse(continuousFormula)
  ziFormula <- deparse(ziFormula)

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
  nullDFList <- generateNULL(modelingData, MetaDF,continuousFormula, ziFormula, family, modality, initialSampling)
  if (verbose) {
    message("Modeling results.")
  }

  if(numCores <= 1){
    stop('numCores must be greater than 1. This method is meant to be parallelized.')
  }
  # Make your clusters for efficient parallelization

  cl <- parallel::makeCluster(numCores)
  parallel::clusterExport(
      cl = cl, varlist = c("continuousFormula", "ziFormula", "modelingData", "MetaDF",
                              "nullDFList","zi_threshold",'family'),
      envir = environment()
    )
  ## If it's scATAC, use glmmTMB. If it's general, use lmerTest instead.
  if(modality == 'scATAC'){
    parallel::clusterEvalQ(cl, {
        library('glmmTMB')
      })
    
    coeffList <- pbapply::pblapply(cl = cl, X = rownames(modelingData), individualZIGLMM)
    parallel::stopCluster(cl)

  }else if(modality == 'General'){
    parallel::clusterEvalQ(cl, {
        library('lmerTest')
      })

    coeffList <- pbapply::pblapply(cl = cl, X = rownames(modelingData), individualLMEM)
    parallel::stopCluster(cl)
  }

  if (verbose) {
      message("Reorganizing residuals and random effect variance.")
  }

  if(modality %in% c('scATAC','scRNA')){
    ranged = TRUE
  }else{
    ranged= FALSE
  }
  processedOuts <- processModelOutputs(modelOutputList = coeffList, 
                                        nullDFList = nullDFList, 
                                        rownamesList = rownames(modelingData),
                                        ranged = ranged,
                                        SummarizedExperimentObj = SE_Object
                                        )
  processedOuts@metadata = append(processedOuts@metadata, list('Type' = modality))
  return(processedOuts)
}


generateNULL <- function(modelingData, MetaDF, continuousFormula, ziFormula, family, modality, initialSampling){

  pilotIndices <- sample(x = 1:dim(modelingData)[1], size = initialSampling, replace = FALSE)
  modelList <- pbapply::pblapply(X = pilotIndices, function(x) {
    df <- data.frame(
      exp = as.numeric(modelingData[x, ]),
      MetaDF, stringsAsFactors = FALSE
    )

    tryCatch(
      {
        if(modality == 'scATAC'){
           glmmTMB::glmmTMB(as.formula(continuousFormula),
            ziformula = as.formula(ziFormula),
            data = df,
            family = family,
            REML = TRUE
          )
        }else if(modality == 'General'){
          lmerTest::lmer(formula = as.formula(continuousFormula), data = df)
        }else{
          NA
        }
       
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
    if(modality == 'scATAC'){
      # Did any of the models have both zero-inflated and continous portions?
      bothZI_Cont <- unlist(lapply(idx, function(x){
        ziDF <- summary(modelList[[x]])$coefficients$zi
        !is.null(ziDF) | sum(dim(ziDF)) != 0

      }))

      if(all(!bothZI_Cont) & !ziFormula %in% c('~ 0', '~0')){
        warning('No working models using the zero-inflated formula. Do you need to modify the zero-inflated formula?')
        modelRes <- modelList[[idx[1]]]
      }else{
        modelRes <- modelList[[idx[which(bothZI_Cont)[1]]]]
      }

      #Extract the first representative model. And use it to create a null template. 
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



    } else if(modality == 'General'){


      modelRes <- modelList[[idx[1]]]

      combinedCoeff <- as.data.frame(summary(modelRes)$coefficients)
      combinedCoeff[!is.na(combinedCoeff)] <- NA 
      Resid <- stats::resid(modelRes)  
      Resid[!is.na(Resid)] <- NA 

      #Add in missing residuals 
      if(!all(MetaDF$Sample %in% names(Resid))){
        NA_samples <- rep(NA, sum(!MetaDF$Sample %in% names(Resid)))
        names(NA_samples) <- MetaDF$Sample[!MetaDF$Sample %in% names(Resid)]
        Resid <- c(Resid, NA_samples)
      }
      Resid <- Resid[match(names(Resid),MetaDF$Sample)]

      #Get variance
      varcor_df <- as.data.frame(lme4::VarCorr(modelRes))$vcov
      names(varcor_df) <- as.data.frame(lme4::VarCorr(modelRes))$grp
      varcor_df[!is.na(varcor_df)] <- NA 
    }

    nullDFList <- list('Coeff' = combinedCoeff, 'Resid' = Resid, 'VCov'= varcor_df)

    rm(modelList)
  }
  return(nullDFList)
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
individualLMEM <- function(x) {
 
  
  df <- data.frame(
    exp = as.numeric(modelingData[x, ]),
    MetaDF, stringsAsFactors = FALSE
  )

   output_vector <- tryCatch(
    {
      modelRes <- lmerTest::lmer(formula = as.formula(continuousFormula), data = df)
      Coeff <- as.data.frame(summary(modelRes)$coefficients)
      Resid <- stats::resid(modelRes)
      
      if(!all(MetaDF$Sample %in% names(Resid))){
        NA_samples <- rep(NA, sum(!MetaDF$Sample %in% names(Resid)))
        names(NA_samples) <- MetaDF$Sample[!MetaDF$Sample %in% names(Resid)]
        Resid <- c(Resid, NA_samples)
      }
      Resid <- Resid[match(names(Resid),MetaDF$Sample)]
      vcov <- as.data.frame(lme4::VarCorr(modelRes))$vcov
      names(vcov) <- as.data.frame(lme4::VarCorr(modelRes))$grp
      list('Coeff' = Coeff, 'Resid' = Resid, 'VCov'= vcov)
    },
    error = function(e) {
      nullDFList
    }
  )
  return(output_vector)
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


individualZIGLMM <- function(x) {
  df <- data.frame(
    exp = as.numeric(modelingData[x, ]),
    MetaDF, stringsAsFactors = FALSE
  )

  output_vector <- tryCatch(
    {
     
      if(sum(df$exp == 0, na.rm = TRUE)/length(df$exp) == 0){
        modelRes <- glmmTMB::glmmTMB(as.formula(continuousFormula),
          ziformula = ~ 0,
          data = df,
          family = family,
          REML = TRUE,
          control = glmmTMB::glmmTMBControl(parallel = 1)
        )
       
      }else if(sum(df$exp == 0, na.rm = TRUE)/length(df$exp) <= zi_threshold){
        df$exp[df$exp == 0] = NA
        modelRes <- glmmTMB::glmmTMB(as.formula(continuousFormula),
          ziformula = ~ 0,
          data = df,
          family = family,
          REML = TRUE,
          control = glmmTMB::glmmTMBControl(parallel = 1)
        )
       
      }else {
        modelRes <- glmmTMB::glmmTMB(as.formula(continuousFormula),
          ziformula = as.formula(ziFormula),
          data = df,
          family = family,
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
        }else if(all(df$exp !=0, na.rm = TRUE)){
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



#' @title Internal function to processing model outputs
#'
#' @description \code{processModelOutputs} 
#' @param modelOutputList. A list of modeloutputs, processed by either individualLMEM or individualZIGLMM. 
#'        The first output is the coefficient data.frame, then the residuals, and the Variance. 
#' @param nullDFList A null templates for the model outputs
#' @param rownamesList a name of all the rows that were interate over. 
#' @param SummarizedExperiment SummarizedExperiment object used for modeling. From this, rowData and colData will be preserved with the residuals.
#' @return A SummarizedExperiment object, that captures the models performance. Each fixed effect will be one assay, columns will be the 
#'        statistics for the fixed effect and measurements (Estimate, Error, p-value, etc..). Residuals and Variance will be saved in the object's metadata.
#'
#' @noRd
processModelOutputs <- function(modelOutputList, nullDFList, rownamesList, ranged = FALSE,
                                  SummarizedExperimentObj, returnList = FALSE) {

    coeffNames <- rownames(nullDFList$Coeff)
    newColumnNames <- gsub('Pr\\(>\\|.\\|)','p_value', gsub(' |\\. ','_',colnames(nullDFList$Coeff)))
    output_list <- lapply(coeffNames, function(z){
      tmpCoef <- do.call("rbind", pbapply::pblapply(X = modelOutputList, function(x) {
            tmpDf <- x[['Coeff']][z,]
            colnames(tmpDf) <- newColumnNames
            tmpDf$FDR <- p.adjust(tmpDf$p_value, 'fdr')
            tmpDf
          }, cl = NULL))
      rownames(tmpCoef) <- rownamesList
      tmpCoef
    })
    names(output_list) <- gsub('Pr\\(>\\|.\\|)','p_value', gsub(' |\\. ','_',coeffNames))

    residual_tmp <- do.call(
      "rbind", pbapply::pblapply(X = modelOutputList, function(x) {
        x[['Resid']]
      }, cl = NULL)
    )
    vcov_tmp <- do.call(
      "rbind", pbapply::pblapply(X = modelOutputList, function(x) {
        x[['VCov']]
    }, cl = NULL)
    )
    rownames(residual_tmp) <- rownames(vcov_tmp) <- rownamesList
    residual_tmp <- residual_tmp[,match(rownames(SummarizedExperiment::colData(SummarizedExperimentObj)), colnames(residual_tmp))]

    if(returnList){
      return(list('output' = output_list , 'Resid' = residual_tmp , 'Variance' = vcov_tmp))
    }
    #Repackage Residuals into a SummarizedExperiment
    if(ranged){
      ResidualSE <- SummarizedExperiment::SummarizedExperiment(
                      list('Residual' =  residual_tmp),
                      colData = SummarizedExperiment::colData(SummarizedExperimentObj),
                      rowRanges = SummarizedExperiment::rowRanges(SummarizedExperimentObj),
                      metadata = S4Vectors::metadata(SummarizedExperimentObj)
              )
      #Package up the metadata list. 
      metaDataList <- list('Residuals' = ResidualSE,
                  'RandomEffectVariance' = vcov_tmp)

      results <- SummarizedExperiment::SummarizedExperiment(
          output_list,
          rowRanges = SummarizedExperiment::rowRanges(SummarizedExperimentObj),
          metadata = metaDataList
      )
    }else{
      ResidualSE <- SummarizedExperiment::SummarizedExperiment(
                  list('Residual' =  residual_tmp),
                  colData = SummarizedExperiment::colData(SummarizedExperimentObj),
                  rowData = SummarizedExperiment::rowData(SummarizedExperimentObj),
                  metadata = S4Vectors::metadata(SummarizedExperimentObj)
          )
      #Package up the metadata list. 
      metaDataList <- list('Residuals' = ResidualSE,
                        'RandomEffectVariance' = vcov_tmp)

     results <- SummarizedExperiment::SummarizedExperiment(
       output_list,
       rowData = SummarizedExperiment::rowData(SummarizedExperimentObj),
       metadata = metaDataList
     )
    
    }

              
    
    return(results)
}
