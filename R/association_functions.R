#' @title \code{GeneTileAssociations}
#'
#' @description \code{GeneTileAssociations} leverages zero-inflated generalized linear mixed effect modeling to associate accessibility
#                       and gene expression.
#' @param rnaSE SummarizedExperiment object from getPseudobulkRNA
#' @param atacSE A SampleTileMatrix object from getSampleTileMatrix or combinedSampleTileMatrix
#' @param cellPopulatioj - name of the assay (cell population) to extract and use from rnaSE and atacSE
#' @param sampleColumn - A string: the name of the column with Sample names from the metadata. 
#'                  Must be the same from atacSE and rnaSE. This is used for aligning atacSE and rnaSE
#' @param continuousFormula : Formula used for modeling. Must be in the form Tile ~ Gene + other factors + (1|RandomEffect)
#' @param ziFormula list of  to tests from atacSE to test. 
#' @param geneList a list of genes from rnaSE to test. 
#' @param distance A number
#' @param pilot a boolean to test just 5 random genes form geneList and return full models for diagnosis, or
#'      test all options and return the coefficients. Default is FALSE, which will test all genes. 
#' @param numCores Optional, the number of cores to use with multiprocessing. Default is 1.
#' 
#' @return A summarized experiment summarizing the output. 
#'
#' @details
#'
#'
#' @export

GeneTileAssociations <- function(atacSE, rnaSE, cellPopulation, sampleColumn, 
                                    continuousFormula, 
                                    ziFormula, 
                                    geneList, 
                                    distance = 10^6,
                                    pilot = FALSE, numCores) {
  
   #Check whether samples align.     
  if(!all(atacSE@colData[,sampleColumn] %in% rnaSE@colData[,sampleColumn])) {
      stop('sampleColumns are not the same. Please ensure that sample names match in atacSE and rnaSE')
  }else if(!all(atacSE@colData[,sampleColumn] == rnaSE@colData[,sampleColumn])){
      warning('Reording atacSE sample names to match rnaSE.')
      rnaSE <- rnaSE[,match(colnames(atacSE), colnames(rnaSE))]
  }

  if(any(c(colnames(SummarizedExperiment::colData(atacSE)), 
      colnames(SummarizedExperiment::colData(rnaSE))) %in% c('exp1','exp2'))){

    stop('metadata of atacSE and/or rnaSE contains a column that contains the name exp1 or exp2.',
    'exp1 and exp2 are hardcoded to represent the data from atacSE and rnaSE, not the metadata.
    Please remove these and try again.')

  }

  geneRanges <- plyranges::filter(SummarizedExperiment::rowRanges(fullSE), GeneSymbol == geneList)


                                    }
  mat1 <- MOCHA::getCellPopMatrix(atacSE, cellPopulation, NAtoZero = TRUE)
  mat2 <- SummarizedExperiment::assays(rnaSE)[[cellPopulation]]
  #Test whether every gene in the geneList are found within the rnaSE
  if(!all(geneList %in% rownames(mat2))){
     stop('geneList not found within rnaSE. Please read documentation.')
  }
  
  #Find all tiles that overlap within the window of the geneTSS. 
  GenomicRanges::promoters(SummarizedExperiment::rowRanges(rnaSE))

  #Test whether the formula is in the right format
  if(!(all.vars(as.formula(formula))[1] =='exp1' & any(all.vars(as.formula(formula)) %in% 'exp2'))){
    stop('formula is not in the correct format. Data from rnaSE will be used to predict atacSE. Please format formula as exp1 ~ exp2 + OtherFixedEffects + RandomEffect.')
  }

  mat1 <- mat1[rownames(mat1) %in% sig1,,drop=FALSE]
  mat2 <- mat2[rownames(mat2) %in% sig2,,drop=FALSE]
  metadata1 <- SummarizedExperiment::colData(atacSE)
  
  # make sure the formula is a string.
  formula = as.character(formula)
    
  # initialize parallel processing
  cl <- parallel::makeCluster(numCores)

  #Find all combinations of sig1 and sig2 to test. 
  allCombo <- as.data.frame(expand.grid(as.character(sig2),as.character(sig1)), stringsAsFactors = FALSE)
    
  if(pilot){
      
      allComboList <- lapply(1:10, function(x){
                   rev(as.character(unlist(allCombo[x,])))
      })
      
      pilotModels <- lapply(allComboList, function(x){
              df <- data.frame(exp1 = unlist(mat1[x[1],]), exp2 = unlist(mat2[x[2],]), metadata1, stringsAsFactors= FALSE)
              lmerTest::lmer(formula = as.formula(formula) , data = df)
        })
      names(pilotModels) <- unlist(lapply(allComboList, function(x) { paste(x[1],x[2],sep = '_')}))
      
      return(pilotModels)
    
    }

  # Generate all combinations of sig1 and sig2 as a list to iterate over.
  allComboList <- pbapply::pblapply(cl = NULL, X = c(1:dim(allCombo)[1]), function(x){
                   rev(as.character(unlist(allCombo[x,])))
      })
    
  individualAssociations <- function(x){
    
    df <- data.frame(exp1 = unlist(mat1[x[[1]],]), exp2 = unlist(mat2[x[[2]],]), metadata1)
    modelRes <- lmerTest::lmer(as.formula(formula), data = df)
    outputDF <- as.data.frame(summary(modelRes)$coefficients)
    rownames(outputDF)[grepl('Intercept', rownames(outputDF))] = 'Intercept'
    outputDF$Modality1 = x[[1]]
    outputDF$Modality2 = x[[2]]

    Resid <- stats::resid(modelRes)
    
    if(!all(metadata1[,sampleColumn] %in% names(Resid))){
      NA_samples <- rep(NA, sum(!metadata1[,sampleColumn] %in% names(Resid)))
      names(NA_samples) <- metadata1[!metadata1[,sampleColumn] %in% names(Resid), sampleColumn]
      Resid <- c(Resid, NA_samples)
    }
    Resid <- Resid[match(names(Resid),metadata1[,sampleColumn])]
    vcov <- as.data.frame(lme4::VarCorr(modelRes))$vcov
    names(vcov) <- as.data.frame(lme4::VarCorr(modelRes))$grp
    return(list('Coeff' = outputDF, 'Resid' = Resid, 'VCov'= vcov))
    
    }
    
    
  #Export data to each core
  parallel::clusterExport(cl, c("mat1", "mat2", "metadata1", "formula","sampleColumn","individualAssociations"), envir = environment())
    
  #Test for each individual combination. 
  allRes <- pbapply::pblapply(cl = cl, X =allComboList, individualAssociations)

  ## Identify fixed effect variables to save, and generate nullDFList for processModelOutputs
  nullDFList <- allRes[[1]]
  nullDFList$Coeff[!is.na(nullDFList$Coeff)] = NA
  nullDFList$Resid[!is.na(nullDFList$Resid)] = NA
  nullDFList$VCov[!is.na(nullDFList$VCov)] = NA
  
  outputList <- processModelOutputs(modelOutputList=allRes, nullDFList = nullDFList, 
                      rownamesList =paste(allCombo$Var2,  allCombo$Var1, sep = "_"),
                                  SummarizedExperimentObj = atacSE, returnList = TRUE) 


  # Process rowData. 
  rowData1 <- as.data.frame(SummarizedExperiment::rowData(atacSE))
  rowData2 <- as.data.frame(SummarizedExperiment::rowData(rnaSE))
  if(any(dim(rowData1) ==0) & any(dim(rowData1) == 0)){
    #No row data available
    allRowData = NULL
  }else if(any(dim(rowData1) ==0)){
    #Only rowData from rnaSE
    allRowData <- rowData1[allCombo$Var2,]
    allRowData$Obj2 = allCombo$Var2
    rownames(allRowData) = paste(allCombo$Var2,  allCombo$Var1, sep = "_")
      
  }else if(any(dim(rowData2) ==0)){
    #Only rowData from atacSE
    allRowData <- rowData1[allCombo$Var1,]
    allRowData$Obj1 = allCombo$Var1
    rownames(allRowData) = paste(allCombo$Var2,  allCombo$Var1, sep = "_")

  }else{

    allRowData1 <- rowData1[allCombo$Var1,]
    allRowData1$Obj1 = allCombo$Var1
    rownames(allRowData1) = paste(allCombo$Var2,  allCombo$Var1, sep = "_")
    allRowData2 <- rowData2[allCombo$Var2,]
    allRowData2$Obj2 = allCombo$Var2
    rownames(allRowData2) = paste(allCombo$Var2,  allCombo$Var1, sep = "_")
    allRowData <- cbind(allRowData1, allRowData2)
  }
  residualDF <- outputList[[2]]
  colnames(residualDF)[is.na(colnames(residualDF))] = rownames(SummarizedExperiment::colData(atacSE))[
      !rownames(SummarizedExperiment::colData(atacSE)) %in% colnames(residualDF)]
  #Create ResidualObject, and metadata list.
  ResidualSE <- SummarizedExperiment::SummarizedExperiment(
                      list('Residual' =  residualDF),
                      colData = SummarizedExperiment::colData(atacSE),
                      rowData = allRowData,
                      metadata = S4Vectors::metadata(atacSE)
              )
              
  #Package up the metadata list. 
  metaDataList <- list('Residuals' = ResidualSE,
                      'RandomEffectVariance' = outputList[[3]])

  #Bind results into one data.frame.
  resDF <- SummarizedExperiment::SummarizedExperiment(outputList[[1]],
                    rowData = allRowData,
                    metadata= metaDataList)


  # stop parallel processing
  parallel::stopCluster(cl)
  # concatenate results and return
  return(resDF)
