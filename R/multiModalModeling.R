#' @title \code{multiModalModeling}
#'
#' @description \code{multiModalModeling} allows for generalized linear mixed effect modeling across multiple modalities,
#'  as long as samples are aligned. This is only designed for integrating datatypes with normal distributions.
#' @param SE1 SummarizedExperiment object for measure-type 1. 
#' @param SE2 SummarizedExperiment object for measure-type 2. 
#' @param cellPopulation A string denoting the cell population of interest, which must be present in SampleTileObj
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

  # Generate all combinations of sig1 and sig2 as a list to iterate over.
  allComboList <- pbapply::pblapply(cl = NULL, X = c(1:dim(allCombo)[1]), function(x){
                   rev(as.character(unlist(allCombo[x,])))
      })
    
  individualAssociations <- function(x){
    
    df <- data.frame(exp1 = unlist(mat1[x[[1]],]), exp2 = unlist(mat2[x[[2]],]), metadata1)
    modelRes <- lmerTest::lmer(formula, data = df)
    outputDF <- as.data.frame(summary(modelRes)$coefficients)
    outputDF$Object1 = x[[1]]
    outputDF$Object2 = x[[2]]
    return(outputDF)
    
    }
    
    
  #Export data to each core
  parallel::clusterExport(cl, c("mat1", "mat2", "metadata1", "formula","individualAssociations"), envir = environment())
    
  #Test for each individual combination. 
  allRes <- pbapply::pblapply(cl = cl, X =allComboList, individualAssociations)
  
  variablesList <- all.vars(lme4::nobars(as.formula(formula)))
  variablesList <- c('(Intercept)',variablesList[!variablesList %in% 'exp1'])
  
  resDFList <- pbapply::pblapply(cl = NULL, X= variablesList, function(x){

      oneVar <- as.data.frame(do.call('rbind',lapply(allRes, function(y){
        tmp1 <- y[x,]
        colnames(tmp1) <- c('Estimate','StdError','df','tValue','pValue','Obj1','Obj2')
        tmp1$FDR <- p.adjust(tmp1$pValue,'fdr')
        tmp1
      })))
      rownames(oneVar) <- paste(oneVar$Obj1, oneVar$Obj2, sep = "_")
      oneVar
  })
  names(resDFList) = variablesList
  #Bind results into one data.frame.
  resDF <- SummarizedExperiment::SummarizedExperiment(resDFList)
  # stop parallel processing
  parallel::stopCluster(cl)
  # concatenate results and return
  return(resDF)
}




visualizeAssociations <- function(interAct, factor1, threshold = 20,  max.overlaps = 10){
    df1 <- as.data.frame(table(interAct[,factor1]))
    colnames(df1) = c('Factor','Interactions')
    df1 <- df1[order(df1$Interactions,decreasing =T),]
    df1$Rank <- as.numeric(c(1:length(df1$Interactions)))
    df1$Factor[df1$Rank > threshold] = NA
    print(head(df1))
    ggplot(df1, aes(x = Rank, y = Interactions, label = Factor)) + ggrastr::geom_point_rast() + theme_minimal() +
    ggtitle(factor1) + 
    xlab('Rank') + ylab('Number of Interactions')  + ggrepel::geom_text_repel(max.overlaps = max.overlaps)
}



volcanoPairs <- function(interAct, factor1 = 'TF', factor2 = 'Protein', 
                         FDR_threshold = 0.05,  topLimit = 20, max.overlaps = 10){
    
     xRange <- range(interAct$Estimate)*0.10
     yRange <- range(-log10(interAct$FDR))*0.10
    interAct$Pair = paste(interAct[,factor1], interAct[,factor2], sep = ":")
    interAct$Pair[interAct$FDR >= FDR_threshold] = NA
    if(sum(interAct$FDR <= FDR_threshold)){
        topPairs <-  interAct$Pair[order(interAct$FDR)[c(1:topLimit)]]
        interAct$Pair[!interAct$Pair %in% topPairs] = NA
    }
    
    ggplot(interAct, aes(x = Estimate, y = -log10(FDR), label = Pair)) + 
    geom_point() +   coord_cartesian(clip = 'off') + 
                theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
    ggrepel::geom_text_repel(max.overlaps = max.overlaps, force = 4,
                            max.iter = 100000) +
    theme_minimal() + 
               
    geom_hline(yintercept = -log10(FDR_threshold), color = 'red', alpha = 0.25) 
}
