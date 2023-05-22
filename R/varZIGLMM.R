#' @title Zero-inflated Variance Decomposition for pseudobulked scATAC data
#'
#' @description \code{varZIGLMM} Identified variance decomposition on a given cell type across both zero-inflated and continous
#'  space using a zero-inflated general linear mixed model \code{\link[glmmTMB]{glmmTMB}}
#'
#' @param TSAM_Object A SummarizedExperiment object generated from
#'   getSampleTileMatrix. 
#' @param cellTypeName Name of a cell type. Should match up exactly with the assay name within SummarizedExperiment. 
#' @param continuousRandom Random effects to test in the continuous portion. All factors must be found in column names
#'   of the TSAM_Object metadata, except for FragNumber and CellCount, which will be extracted from the TSAM_Object's metadata.
#' @param ziformula Random effects to test in the zero-inflated portion. All factors must be found in column names
#'   of the TSAM_Object colData metadata, except for FragNumber and CellCount, which will be extracted from the TSAM_Object's metadata.
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @param numCores integer. Number of cores to parallelize across.
#'
#' @return results a SummarizedExperiment containing results from ZIGLMM (Fixed effect estiamtes, P-values, and Std Error)
#'
#'
#'
#' @examples
#' \dontrun{
#'   modelList <- runZIGLMM(STM[c(1:1000),], 
#'                  cellTypeName = 'CD16 Mono',
#'                  continuousRandom = c('Age', 'Sex', 'Days'),
#'                  ziRando = c('FragNumber', 'Days'),
#'                  verbose = TRUE, 
#'                  numCores = 35 )
#' }
#'
#' @export
#' 
varZIGLMM <- function(TSAM_Object,
                      cellTypeName = NULL,
                      continuousRandom = NULL,
                      ziRandom = NULL,
                      verbose = FALSE,
                      numCores = 1) {


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
  
  if (!all(continuousRandom %in% c("exp", colnames(MetaDF)))) {
    stop("Random continuous effects are not found in metadata.")
  }

  if (!all(ziRandom %in% colnames(MetaDF) | all(ziRandom == 0))) {
    stop("Random Zero-inflated effects are not found in metadata.")
  }

  CondvarForm <- paste0(unlist(lapply(continuousRandom, function(x) paste('(1|',x,')',sep=''))), collapse = ' + ')
  continuousFormula <- paste('exp ~ ',CondvarForm, sep = '')

  if(all(ziRandom != 0)){
    zi_form <- paste0(unlist(lapply(ziRandom, function(x) paste('(1|',x,')',sep=''))), collapse = ' + ')
    variableList <- c(continuousRandom, ziRandom)
  }else{
    zi_form = ziRandom
    variableList <- continuousRandom
  }
  ziformula = paste("~ ", zi_form)

  MetaDF <- dplyr::filter(MetaDF, Sample %in% colnames(modelingData))
  modelingData <- modelingData[, match(colnames(modelingData), MetaDF$Sample)]

  # Subset metadata to just the variables needed. This minimizes overhead for parallelization
  MetaDF <- MetaDF[, colnames(MetaDF) %in% c("Sample", variableList)]

  if (any(is.na(MetaDF))) {
    stop("NAs are included in the MetaDF. Please remove them and try again.")
  }


  #Generate null results.
  if(all(ziRandom != 0)){
    nullDF <- rep(NA, length(continuousRandom) + length(ziRandom) + 1)
    names(nullDF) <- c(paste('Cond', continuousRandom, sep="_"), paste('ZI', ziRandom, sep="_"), 'Residual')
  }else{
    nullDF <- rep(NA, length(continuousRandom) + 1)
    names(nullDF) <- c(paste('Cond', continuousRandom, sep="_"), 'Residual')
  }
  
  # Subset metadata to just the variables needed. This minimizes overhead for parallelization
  MetaDF <- MetaDF[, colnames(MetaDF) %in% c("Sample", variableList)]

  if (any(is.na(MetaDF))) {
    stop("NAs are included in the MetaDF. Please remove them and try again.")
  }
  
  # Make your clusters for efficient parallelization
  if(numCores > 1){
      cl <- parallel::makeCluster(numCores)
      parallel::clusterEvalQ(cl, {
                library(glmmTMB)
        })
  
    parallel::clusterExport(
    cl = cl, varlist = c("continuousFormula", "ziformula", "modelingData", "MetaDF", "individualVarZIGLMM", "nullDF"),
    envir = environment()
  )
      
  }else{

    cl = NULL
  }
   
  browser()
  varDecompList <- pbapply::pblapply(cl = cl, X = rownames(modelingData), individualVarZIGLMM)

  if(!is.null(cl)){
    parallel::stopCluster(cl)
  }
                                    
  results <- do.call('rbind', varDecompList)
  rownames(results) <- rownames(modelingData)

  return(results)
}


#' @title \code{individualVarZIGLMM }
#'
#' @description \code{individualVarZIGLMM } Runs ZI-GLMM on data provided and returns variance decomposition across random effects. Wirtten for parallelization.
#'
#' @param refList. A list where the first index is a data.frame to use for modeling, and the second is the formula for modeling.
#'
#' @return A linear model

#' @noRd
#'
#'

individualVarZIGLMM <- function(x) {
    
    
    df <- data.frame(
        exp = as.numeric(modelingData[x, ]),
        MetaDF, stringsAsFactors = FALSE
      )

    
  output_vector <- tryCatch(
    {
      if(all(df$exp !=0)){
        modelRes <- glmmTMB::glmmTMB(as.formula(continuousFormula),
          ziformula = ~ 0,
          data = df,
          family = stats::gaussian(),
          REML = TRUE
        )
      }else{
        modelRes <- glmmTMB::glmmTMB(as.formula(continuousFormula),
          ziformula = as.formula(ziformula),
          data = df,
          family = stats::gaussian(),
          REML = TRUE
        )
      }

      if(!modelRes$sdr$pdHess){
        return(nullDF)
      }


      cond_other = unlist(glmmTMB::VarCorr(modelRes)$cond)
      names(cond_other) = paste('Cond', names(cond_other), sep = "_")
      residual = as.vector(attr(glmmTMB::VarCorr(modelRes)$cond, "sc")^2)
      names(residual) = 'Residual'

      if(!is.null(glmmTMB::VarCorr(modelRes)$zi)){
        zi_other = unlist(glmmTMB::VarCorr(modelRes)$zi)
        names(zi_other) = paste('ZI', names(zi_other), sep = "_")
        varcor_df <- c(cond_other, zi_other,residual)
      }else if(all(df$exp !=0)){
        subNull = nullDF[grepl('ZI_', names(nullDF))]
        zi_other = rep(0, length(subNull))
        names(zi_other) = names(subNull)
        varcor_df <- c(cond_other, zi_other,residual)
      }else {
        varcor_df <- c(cond_other, residual)
      }

      varDecomp <- varcor_df/sum(varcor_df)
      return(varDecomp)

    },
    error = function(e) {
      nullDF
    }
  )
  return(output_vector)
}
