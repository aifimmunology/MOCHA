#' @title \code{multiModalModeling}
#'
#' @description \code{multiModalModeling} allows for generalized linear mixed effect modeling across multiple modalities,
#'  as long as samples are aligned. This is only designed for integrating datatypes with normal distributions.
#' @param SE1 SummarizedExperiment object for measurement-type 1. 
#' @param SE2 SummarizedExperiment object for measurement-type 2. 
#' @param assay1 - name of the assay to extract and use from SE1
#' @param assay2 - name of the assay to extract and use from SE2
#' @param sampleColumn - A string: the name of the column with Sample names from the metadata. 
#'                  Must be the same from SE1 and SE2. This is used for aligning SE1 and SE2
#' @param formula : Formula used for modeling. Must be in the form exp1 ~ exp2 + other factors + (1|RandomEffect)
#' @param sig1 A list of rows from SE1 to test. 
#' @param sig2 a list of rows from SE2 to test. 
#' @param pilot a boolean to test just 10 combinations and return full models for diagnosis, or
#'      test all options and return the coefficients. 
#' @param numCores Optional, the number of cores to use with multiprocessing. Default is 1.
#' 
#' @return A summarized experiment summarizing the output. 
#'
#' @details
#'
#'
#' @export

multiModalModeling <- function(SE1, SE2, assay1, assay2, sampleColumn, formula, sig1, sig2, pilot = FALSE, numCores) {
  
   #Check whether samples align.     
  if(!all(SE1@colData[,sampleColumn] %in% SE2@colData[,sampleColumn])) {
      stop('sampleColumns are not the same. Please ensure that sample names match in SE1 and SE2')
  }else if(!all(SE1@colData[,sampleColumn] == SE2@colData[,sampleColumn])){
      warning('Reording SE1 sample names to match SE2.')
      SE2 <- SE2[,match(colnames(SE1), colnames(SE2))]
  }

  if(any(c(colnames(SummarizedExperiment::colData(SE1)), 
      colnames(SummarizedExperiment::colData(SE2))) %in% c('exp1','exp2'))){

    stop('metadata of SE1 and/or SE2 contains a column that contains the name exp1 or exp2.',
    'exp1 and exp2 are hardcoded to represent the data from SE1 and SE2, not the metadata.
    Please remove these and try again.')

  }
  mat1 <- SummarizedExperiment::assays(SE1)[[assay1]]
  mat2 <- SummarizedExperiment::assays(SE2)[[assay2]]
   #Test whether sig1 and sig2 are found in their assays.
  if(!all(sig2 %in% rownames(mat2))){
     stop('sig2 not found within SE2. Please read documentation.')
  }
    
  #Test whether sig1 and sig2 are found in their assays.
  if(!all(sig1 %in% rownames(mat1))){
     stop('sig1 not found within SE1. Please read documentation.')
  }

  #Test whether the formula is in the right format
  if(!(all.vars(as.formula(formula))[1] =='exp1' & any(all.vars(as.formula(formula)) %in% 'exp2'))){
    stop('formula is not in the correct format. Data from SE2 will be used to predict SE1. Please format formula as exp1 ~ exp2 + OtherFixedEffects + RandomEffect.')
  }

  mat1 <- mat1[rownames(mat1) %in% sig1,,drop=FALSE]
  mat2 <- mat2[rownames(mat2) %in% sig2,,drop=FALSE]
  metadata1 <- SummarizedExperiment::colData(SE1)
  
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
                                  SummarizedExperimentObj = SE1, returnList = TRUE) 


  # Process rowData. 
  rowData1 <- as.data.frame(SummarizedExperiment::rowData(SE1))
  rowData2 <- as.data.frame(SummarizedExperiment::rowData(SE2))
  if(any(dim(rowData1) ==0) & any(dim(rowData1) == 0)){
    #No row data available
    allRowData = NULL
  }else if(any(dim(rowData1) ==0)){
    #Only rowData from SE2
    allRowData <- rowData1[allCombo$Var2,]
    allRowData$Obj2 = allCombo$Var2
    rownames(allRowData) = paste(allCombo$Var2,  allCombo$Var1, sep = "_")
      
  }else if(any(dim(rowData2) ==0)){
    #Only rowData from SE1
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
  colnames(residualDF)[is.na(colnames(residualDF))] = rownames(SummarizedExperiment::colData(SE1))[
      !rownames(SummarizedExperiment::colData(SE1)) %in% colnames(residualDF)]
  #Create ResidualObject, and metadata list.
  ResidualSE <- SummarizedExperiment::SummarizedExperiment(
                      list('Residual' =  residualDF),
                      colData = SummarizedExperiment::colData(SE1),
                      rowData = allRowData,
                      metadata = S4Vectors::metadata(SE1)
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
}





#' @title \code{visualizeAssociations}
#'
#' @description \code{visualizeAssociations} Takes the output of multiModalModeling (or other) and visualizes the top interaction features
#'  for one of the modalities. 
#' @param interAct A dataframe of coefficients related to a modality, extracted from multiModalModeling or other. This should be prefiltered for significance (FDR < 0.1)
#' @param dataType1 Boolean flag. This determines whether it should look at the Modality1 or Modality2 column (i.e. the first or second assay given to multiModalIntegration)
#' @param dataTypeName A string to label the modality in the plots generated. 
#' @param threshold A number, describing how many points to label. By default, it will label the top 20 points. 
#' @param max.overlaps The maximum number of overlapping labels, a parameter which is passed to ggrepel when labeling points on the plot. 
#' @param verbose Boolean flag to determine verbosity. 
#' 
#' @return A ggplot object, showing number of interactions on the y-axis and the rank order of measurements on the x-axis 
#'    (Measurements rank-ordered by the number of interactions they have with the other modality.)
#'
#' @details
#'
#'
#' @export

visualizeAssociations <- function(interAct, dataType1 = TRUE, dataTypeName = 'TF', threshold = 20,  max.overlaps = 10,
                                 verbose = FALSE){
    if(dataType1){
        factor1 = 'Modality1'
    }else{
        factor1 = 'Modality2'
    }
    df1 <- as.data.frame(table(interAct[,factor1]))
    colnames(df1) = c(dataTypeName,'Interactions')
    df1 <- df1[order(df1$Interactions,decreasing =T),]
    df1$Rank <- as.numeric(c(1:length(df1$Interactions)))
    df1$Factor <- df1[,dataTypeName]
    df1$Factor[df1$Rank > threshold] = NA
    if(verbose){
       head(df1)
    }
    ggplot(df1, aes(x = Rank, y = Interactions, label = Factor)) + ggrastr::geom_point_rast() + theme_minimal() +
    ggtitle(factor1) + 
    xlab('Rank') + ylab('Number of Interactions')  + ggrepel::geom_text_repel(max.overlaps = max.overlaps)
}


#' @title \code{volcanoPairs}
#'
#' @description \code{volcanoPairs} Takes the output of multiModalModeling (or other) and generates a volcano plot of Estimate vs -log10(FDR) for that pair
#' @param interAct A dataframe of coefficients related to a modality, extracted from multiModalModeling or other. This should be not be filtered. 
#' @param factor1 Name for the 
#' @param dataTypeName A string to label the modality in the plots generated. 
#' @param FDR_threshold A number, which thresholds which points are considered significant. It will also add a red line to the plot to mark significance.
#' @param topLimit  A number, limiting how many points to label. By default, it will label the top 20 points by FDR, if there are more than 20 points.
#' @param max.overlaps The maximum number of overlapping labels, a parameter which is passed to ggrepel when labeling points on the plot. 
#' @param addSignificanceLine Boolean flag to determine whether or not to add a line at the FDR significance threshold. 
#' 
#' @return A ggplot object, showing number of interactions on the y-axis and the rank order of measurements on the x-axis 
#'    (Measurements rank-ordered by the number of interactions they have with the other modality.)
#'
#' @details
#'
#'
#' @export

volcanoPairs <- function(interAct, FDR_threshold = 0.05,  topLimit = 20, max.overlaps = 10, addSignificanceLine = TRUE){

    #Generate range for plotting
    xRange <- range(interAct$Estimate)*0.10
    yRange <- range(-log10(interAct$FDR))*0.10
    #Generate labels for each plot
    interAct$Pair = paste(interAct[,'Modality1'], interAct[,'Modality2'], sep = ":")
    #Filter out which pairs to label. First by FDR threshold, then by the topLimit parameter. 
    interAct$Pair[interAct$FDR >= FDR_threshold] = NA
    if(sum(interAct$FDR <= FDR_threshold)){
        topPairs <-  interAct$Pair[order(interAct$FDR)[c(1:topLimit)]]
        interAct$Pair[!interAct$Pair %in% topPairs] = NA
    }

    #Generate ggplot.  
    p1 <- ggplot(interAct, aes(x = Estimate, y = -log10(FDR), label = Pair)) + 
           geom_point() +   coord_cartesian(clip = 'off') + 
                theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
            ggrepel::geom_text_repel(max.overlaps = max.overlaps, force = 4,
                            max.iter = 100000) +
    theme_minimal() 
    if(addSignificanceLine){
      p1 <- p1+ 
            geom_hline(yintercept = -log10(FDR_threshold), color = 'red', alpha = 0.25) 
    }
    return(p1)

}