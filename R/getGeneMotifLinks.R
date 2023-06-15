#' @title \code{getGeneMotifLinks}
#'
#' @description \code{getGeneMotifLinks} models which chromVAR Z-scores predict gene expression, according to the formula variable. 
#' @param rnaSE SummarizedExperiment object of pseudobulked RNA counts (normalized). Designed to be downstream of makePseodulkRNA.
#' @param motifSE SummarizedExperiment object for ChromVAR z-scores. Will work for any normally distributed dataset put into a SummarizedExperiment format (e.g. Olink measurements).
#' @param cellPopulation A string denoting the cell population of interest, which must be present in both rnaSE and motifSE
#' @param sampleColumn A string denoting the name of the column in each SummarizedExperiment-type objec that has sample data.
#' @param formula A formula (or string describing a formula) that will be used for testing. Formula should be in the form Gene ~ Motif + OtherVariables + (1|RandomEffects).
#'                The sample metadata from rnaSE cannot contain a column called Gene or Motif. These are hardcoded to be data from rnaSE and motifSE.
#' @param sigGenes A list of gene names to test. Must aligned with gene names from the rnaSE. Should be differential.
#' @param sigMotifs a list of ChromVAR Z-scores to test. 
#' @param pilot A boolean to describe whether to model all combinations or run a pilot case of 10 combinations, for the purpose of testing the model.
#' @param numCores Optional, the number of cores to use with multiprocessing. Default is 1.
#' @return foregroundDF A data.table correlation matrix
#'
#' @details The technical details of the zero-inflated correlation can be
#'          found here:
#'
#'               Pimentel, Ronald Silva, "Kendall's Tau and Spearman's Rho
#'               for Zero-Inflated Data" (2009). Dissertations.
#'
#'          while the implementation (scHOT R package), can be found here:
#'               http://www.bioconductor.org/packages/release/bioc/html/scHOT.html
#'
#'
#' @export

getGeneMotifLinks <- function(rnaSE, motifSE, cellPopulation,sampleColumn, formula, rnaSig, motifSig, pilot = FALSE, numCores) {
  
   #Check whether samples align.     
  if(!all(rnaSE@colData[,sampleColumn] %in% motifSE@colData[,sampleColumn])) {
      stop('sampleColumns are not the same. Please ensure that sample names match in rnaSE and motifSE')
  }else if(!all(rnaSE@colData[,sampleColumn] == motifSE@colData[,sampleColumn])){
      warning('Reordering rnaSE sample names to match motifSE.')
      motifSE <- motifSE[,match(colnames(rnaSE), colnames(motifSE))]
  }
  
  geneMat <- SummarizedExperiment::assays(rnaSE)[[cellPopulation]]
  motifMat <- SummarizedExperiment::assays(motifSE)[[cellPopulation]]

  #Test whether rnaSig and motifSig are found in their assays.
  if(all(motifSig %in% rownames(motifMat))){
     stop('motifSig not found within motifSE. Please read documentation.')
  }
    
  #Test whether rnaSig and motifSig are found in their assays.
  if(all(rnaSig %in% rownames(geneMat))){
     stop('rnaSig not found within rnaSE. Please read documentation.')
  }
  geneMat <- geneMat[rownames(geneMat) %in% rnaSig,,drop = FALSE]
  motifMat <- motifMat[rownames(motifMat) %in% motifSig,,drop = FALSE]
  metadata1 <- SummarizedExperiment::colData(rnaSE)
  
 
  # make sure the formula is a string.
  formula = as.character(formula)
    
  # initialize parallel processing
  cl <- parallel::makeCluster(numCores)

  #Find all combinations of rnaSig and motifSig to test. 
  allCombo <- as.data.frame(expand.grid(as.character(motifSig),as.character(rnaSig)), stringsAsFactors = FALSE)
    
  if(pilot){
      
      allComboList <- lapply(1:10, function(x){
                   rev(as.character(unlist(allCombo[x,])))
      })
      
      pilotModels <- lapply(allComboList, function(x){
              df <- data.frame(Gene = unlist(geneMat[x[1],]), Motif = unlist(motifMat[x[2],]), metadata1, stringsAsFactors= FALSE)
               modelRes <- glmmTMB::glmmTMB(as.formula(continuousFormula),
                ziformula = ~0,
                data = df,
                family = glmmTMB::nbinom2(),
                REML = TRUE,
                control = glmmTMB::glmmTMBControl(parallel = numCores)
              )
        })

      names(pilotModels) <- unlist(lapply(allComboList, function(x) { paste(x[1],x[2],sep = '_')}))
      
      return(pilotModels)
    
    }

  # Generate all combinations of rnaSig and motifSig as a list to iterate over.
  allComboList <- pbapply::pblapply(cl = NULL, X = 1:dim(allCombo)[1], function(x){
                   rev(as.character(unlist(allCombo[x,])))
      })
  rownamesList = unlist(pbapply::pblapply(cl = NULL, X = allComboList, function(x){ paste0(unlist(x), collapse = '_')}))

  #Export data to each core
  parallel::clusterExport(cl, c("geneMat", "motifMat", "metadata1", "formula","individualGeneMotifModel"), envir = environment())
    
  #Test for each individual combination. 
  allRes <- pbapply::pblapply(cl = cl, X =allComboList, individualGeneMotifModel)
  
  parallel::stopCluster(cl)

  processedOuts <- processModelOutputs(processModelOutputs = allRes, 
                                        nullDFList = nullDFList, 
                                        rownamesList = rownamesList,
                                        SummarizedExperimentObj = newObj
                                        )

  return(processedOuts)
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

individualGeneMotifModel <- function(x) {
    
  df <- data.frame(Gene = unlist(geneMat[x[1],]), Motif = unlist(motifMat[x[2],]), metadata1, stringsAsFactors= FALSE)  

  output_vector <- tryCatch(
    {
      
      modelRes <- glmmTMB::glmmTMB(as.formula(continuousFormula),
          ziformula = ~0,
          data = df,
          family = glmmTMB::nbinom2(),
          REML = TRUE
        )
    
      if(!modelRes$sdr$pdHess){
        return(nullDF)
      }

      Coeff <-  as.data.frame(summary(modelRes)$coefficients$Cond)
      rownames(Coeff$cond)[grepl('(Intercept)',rownames(Coeff$cond))] = 'Intercept'
      Resid <- stats::resid(modelRes)
      
      #Clean up residuals. Data that is NA will be removed from the residuals, so we add them back in as NAs for the sake of completion. 
      if(!all(MetaDF$Sample %in% names(Resid))){
        NA_samples <- rep(NA, sum(!MetaDF$Sample %in% names(Resid)))
        names(NA_samples) <- MetaDF$Sample[!MetaDF$Sample %in% names(Resid)]
        Resid <- c(Resid, NA_samples)
      }
      Resid <- Resid[match(names(Resid),MetaDF$Sample)]

      cond_other = unlist(glmmTMB::VarCorr(modelRes)$cond)
      names(cond_other) = paste('Cond', names(cond_other), sep = "_")
      residual = as.vector(attr(glmmTMB::VarCorr(modelRes)$cond, "sc")^2)
      names(residual) = 'Residual'
      varcor_df <- c(cond_other, residual)

      return(list('Ceoff' = Coeff, 'Resid' = Resid, 'VCov' = varcor_df))

    },
    error = function(e) {
      nullDF
    }
  )
  return(output_vector)
}