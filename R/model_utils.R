



#' @title Identify regions that are strongly related to the sequencing depth of a given pseudobulked cell type.
#'
#' @description \code{identifyDropOut} 
#' 
#' @param TSAM_Object A SummarizedExperiment object generated from
#'   getSampleTileMatrix. 
#' @param cellPopulation Name of a cell type within the TSAM_Object
#' @param zi_threshold Zero-inflated threshold ( range = 0-1), representing the fraction of samples with zeros. 
#'         This function will only test for technical dropout if the function is above this threshold. 
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @param numCores integer. Number of cores to parallelize across.
#'
#' @return results a SummarizedExperiment containing ZI-GLMM results
                                       
identifyDropOut <- function(TSAM_Object, cellPopulation, zi_threshold = 0.1, verbose = FALSE, numCores = 2){

    mat1 <- getCellPopMatrix(TSAM_Object, cellPopulation)
    trueRows <- rownames(mat1)[rowSums(mat1 == 0)/dim(mat1)[2] >= zi_threshold]

    if(length(trueRows) == 0){

        stop('No locations passed the zi_threshold')

    }

    output <- runZIGLMM(TSAM_Object[trueRows,], cellPopulation, 
                       continuousFormula = exp ~ (1|Sample),
                       ziformula = ~ 0 + FragNumber,
                       verbose = verbose,
                       numCores = numCores)
    annotatedRanges = rowRanges(output) 
    rowValues <- getModelValues(output, 'ZI_FragNumber')    

    annotateRanges$DropOutFDR = rowValues$FDR
    annotateRanges$DropOutPValue = rowValues$p_value
    fullRange = join_by_overlap(rowRanges(STM), annotatedRanges)

    return(fullRange)                
}

modelPredictions <- function(SE1, assay1, measurement, variable, sampleColumn, modelObject, adjust = TRUE){

     #Check whether sampleColumn is present in both SE objects. 
    if(!(sampleColumn %in% colnames(SE1@colData))){
        stop('sampleColumn not found in SE1.')
    }

    #Test whether the measurement is found within their assays
    if(!all(measurement %in% rownames(SE1))){
     stop('measurement not found within SE1. Please read documentation.')
    }

    mat1 <- SummarizedExperiment::assays(SE1[measurement,])[[assay1]]
    mat1[is.na(mat1)] = 0
    assayList <- SummarizedExperiment::assays(modelObject)
    metaData <-  SummarizedExperiment::colData(SE1)


    allVariables <- names(assayList)[! names(assayList) %in% c('Intercept') & !grepl('ZI_',names(assayList))]
    numericVariables <- names(assayList)[names(assayList) %in% colnames(metaData)]

    remainingVariables <- allVariables[!allVariables %in% numericVariables]

    #Create a metadata column for each categorical variable, one-hot encoding them. 
    if(length(remainingVariables) > 1){

        nextVariables <- lapply(remainingVariables, function(x) matchCategorical(metaData, x))
        newMetaData = do.call('cbind', nextVariables)


    }else if(length(remainingVariables) > 0){
        newMetaData <- data.frame(newVar = matchCategorical(metaData, remainingVariables))
    }

    colnames(newMetaData) <- remainingVariables
    metaData = cbind(metaData, newMetaData)

    subMeta <- metaData[, colnames(metaData) %in% allVariables]
    if(any(is.na(subMeta))){

        stop('Metadata in SE1 is NA. SE1 is likely not the object used for modeling. Please provide the object used for modeling.')

    }

    if(modelObject@metadata$Type == 'scATAC'){
        logTransform = TRUE
    }else{
        logTransform = FALSE
    }

    allPredictions <- pbapply::pblapply(cl = NULL, seq_along(measurement), function(x){
        
        modelVals <- getModelValues(modelObject, measurement[[x]])[,'Estimate', drop=FALSE]
        if(logTransform){
            newData <- data.frame('exp1' =  log2(unlist(mat1[measurement[[x]],])+1))
        }else{
            newData <- data.frame('exp1' =  unlist(mat1[measurement[[x]],]))
        }

        rownames(newData) = colnames(mat1)
        
        if(adjust){
            newData$orig_exp1 = newData$exp1
            allAdjusts = do.call('cbind',lapply(allVariables[allVariables != variable], function(x){

                            as.data.frame( modelVals[x,]*as.numeric(subMeta[,x]))

                            }))
            newData$exp1 = newData$exp1 - rowSums(allAdjusts)
            newData$Prediction = modelVals['Intercept',] + modelVals[variable,]*subMeta[,variable]
           
        
        }else{
    
            numericCalls <- sum(unlist(lapply(numericVariables[numericVariables != variable], function(x) {
                                modelVals[x,]*mean(subMeta[,x])
                })))
            newData$Prediction = modelVals['Intercept',] + modelVals[variable,]*subMeta[,variable] + numericCalls
        }

        if(modelObject@metadata$Type == 'scATAC'){
                newData$exp1[newData$orig_exp1 ==0] = NA
                newData$Prediction[newData$orig_exp1 ==0] = NA
        }

        newData <- cbind(newData, metaData)
        newData
        
    })
        
    names(allPredictions) = measurement
                                
    exp = SummarizedExperiment::SummarizedExperiment(allPredictions, metadata = SummarizedExperiment::rowData(SE1[measurement,]))
   
    return(exp)
}


modelInterPredictions <- function(SE1, SE2, assay1, assay2, sampleColumn, interModel, firstPair, secondPair, adjust = TRUE){

     #Check whether sampleColumn is present in both SE objects. 
    if(!(sampleColumn %in% colnames(SE1@colData) & sampleColumn %in% colnames(SE2@colData))){
        stop('sampleColumn not found in SE1 and/or SE2.')
    }

    #Check whether samples align.     
    if(!all(SE1@colData[,sampleColumn] %in% SE2@colData[,sampleColumn])) {
      stop('sampleColumns are not the same. Please ensure that sample names match in SE1 and SE2')
    }else if(!all(SE1@colData[,sampleColumn] == SE2@colData[,sampleColumn])){
      warning('Reording SE1 sample names to match SE2.')
      SE2 <- SE2[,match(colnames(SE1), colnames(SE2))]
    }

    #Test whether the pair is found within their assays
    if(!all(firstPair %in% rownames(SE1))){
     stop('firstPair not found within SE1. Please read documentation.')
    }
    if(!all(secondPair %in% rownames(SE2))){
     stop('secondPair not found within SE2. Please read documentation.')
    }

    pair = paste(firstPair, secondPair, sep ='_')
    #Test whether the pair was actually modeled
    if(!all(pair %in% rownames(interModel))){
     stop('Pair of interacting features not found within interModel. This pair of features was not actually modeled.')
    }

   # if(length(pair) == 1){
   #     
   #     pair = list(pair)
   #    firstPair = list(firstPair)
    #    secondPair = list(secondPair)
   # }

    if(any(c(colnames(SummarizedExperiment::colData(SE1)), 
      colnames(SummarizedExperiment::colData(SE2))) %in% c('exp2'))){

    stop('metadata of SE1 and/or SE2 contains a column that contains the name exp2.',
    'exp2 are hardcoded to represent the data from SE1 and SE2, not the metadata.
    Please remove these and try again.')

    }
    mat1 <- SummarizedExperiment::assays(SE1[unique(firstPair),])[[assay1]]
    mat2 <- SummarizedExperiment::assays(SE2[unique(secondPair),])[[assay2]]

    assayList <- SummarizedExperiment::assays(interModel)
    metaData <-  SummarizedExperiment::colData(SE1)
    allVariables <- names(assayList)[! names(assayList) %in% c('exp2', 'Intercept')]
    numericVariables <- names(assayList)[names(assayList) %in% colnames(metaData)]

    remainingVariables <- allVariables[!allVariables %in% numericVariables]

    #Create a metadata column for each categorical variable, one-hot encoding them. 
    if(length(remainingVariables) > 1){

        nextVariables <- lapply(remainingVariables, function(x) matchCategorical(metaData, x))
        newMetaData = do.call('cbind', nextVariables)


    }else if(length(remainingVariables) > 0){
        newMetaData <- data.frame(newVar = matchCategorical(metaData, remainingVariables))
    }
   
    colnames(newMetaData) <- remainingVariables
    metaData = cbind(metaData, newMetaData)

    subMeta <- metaData[, colnames(metaData) %in% allVariables]
                         
    allPredictions <- pbapply::pblapply(cl = NULL, seq_along(pair), function(x){
        
        modelVals <- getModelValues(interModel, pair[[x]])[,'Estimate', drop=FALSE]

        bothData <- data.frame('exp1' =  unlist(mat1[firstPair[[x]],]), 'exp2' = unlist(mat2[secondPair[[x]],]))
        rownames(bothData) = colnames(mat1)
        
        if(adjust){

            bothData$orig_exp1 = bothData$exp1
            allAdjusts = do.call('cbind',lapply(allVariables, function(x){

                            as.data.frame( modelVals[x,]*as.numeric(subMeta[,x]))

                            }))

            bothData$exp1 = bothData$exp1 - rowSums(allAdjusts)
            bothData$Prediction = modelVals['Intercept',] + modelVals['exp2',]*bothData$exp2
        
        }else{
    
            numericCalls <- sum(unlist(lapply(numericVariables, function(x) {
                                modelVals[x,]*mean(subMeta[,x])
                })))
            bothData$Prediction = modelVals['Intercept',] + modelVals['exp2',]*bothData$exp2 + numericCalls
            
        }
        bothData <- cbind(bothData, metaData)
        bothData
        
    })
        
    names(allPredictions) = pair
                                
    exp = SummarizedExperiment::SummarizedExperiment(allPredictions)
   
    return(exp)
}
                    
matchCategorical <- function(metaData, variable, spacer = '.'){

    whichMatch = which(unlist(lapply(colnames(metaData), function(x) grepl(x, variable))))
        
    if(length(whichMatch) > 1){

       #Iterate over all possible matching columns, and see if any combinations of column name and 
       #values match the variable. 
       allMatches <- lapply(whichMatch, function(x){

            values = gsub(" ",spacer,unlist(metaData[, whichMatch]))
            valuesList <- paste(colnames(metaData)[whichMatch], values, sep ='')
            if(!any(valuesList == variable)){

               NA

            }else{

              valuesList == variable
            } 
           
        })
        #See if multiple columns match. If multiple match, then see if they are all the same or different.
        # if different, throw an error. If the same or only one match, then use it. 
        whichAllMatch <- unlist(lapply(allMatches, function(x) all(is.na(x))))
        if(sum(!whichAllMatch) > 1){

           combMatch <- do.call('rbind', allMatches[!whichAllMatch])
           if(any(unlist(lengths(apply(combMatch, 2, unique))) > 1)){
           
                stop(paste('The categorical variable ', variable,
                    ' could match with multiple different metadata columns. Please clean up metadata.', sep=''))
        
           }else{
               
                return(allMatches[[1]])
                
            }

        }else if(sum(!whichAllMatch) == 1){
            
          return(allMatches[[!whichAllMatch]])
            
        }else{
            
          stop(paste('The categorical variable ', variable, ' could not be matched to metadata', sep=''))
            
        }


   }else if(length(whichMatch) == 1){
       
       #If only one column matches, then pull that one out. 
        values = gsub(" ",spacer,unlist(metaData[, whichMatch]))
        valuesList <- paste(colnames(metaData)[whichMatch], values, sep ='')
       if(!any(valuesList == variable)){
            
           stop(paste('The categorical variable ', variable, ' could not be matched to metadata', sep=''))
           
       }else{
           
           outputList <-  valuesList == variable
           return(outputList)
       }
       
       
   }else{
       
       stop('Variable not found in metadata. Please ensure all variables for the model are within the metadata of SE1.')
       
   }
    
}