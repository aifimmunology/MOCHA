#' @title \code{runChromVAR}
#'
#' @description \code{runChromVAR} Runs ChromVAR on a MOCHA Tile-Sample Accessibility Matrix object
#'
#'
#' @param TSAM_Object RangedSummarizedExperiment object, containing the TSAM from getSampleTileMatrix, with a motif set added via addMotifSet()
#' @param motifSetName name of the motifList in the SE_object's metadata. Must be present to run chromVAR
#' @param groupList the list of cell type names or sample-associated data that
#'   should be used to subset the Object
#' @param na.rm removes groups that are NA if set to true. If set to false, then
#'   you filter for everything in the groupList and also NA values.
#'
#' @return chromVAR object
#'
#' @details This is a wrapper for basic SummarizedExperiment-based subsetting. 
#'
#'
#' @export
#' 

runChromVAR <- function(TSAM_Object,
                              motifSetName,
                              cellTypeSpecific = FALSE){

    if(!cellTypeSpecific){

        if(length(assays(TSAM_Object)) > 1){
            assayList <- do.call('cbind', SummarizedExperiment::assays(TSAM_Object))
            newSamplesNames <- unlist(lapply(names(SummarizedExperiment::assays(TSAM_Object)), function(x){
        
                            paste(x, colnames(TSAM_Object), sep = "_")
                    
            }))

            colnames(assayList) <- newSamplesNames

            nonEmptySamples <- apply(assayList, 2, function(x) !all(is.na(x)))
            assayList <- assayList[,nonEmptySamples]

            allSampleData <- do.call('rbind', lapply(1:length(SummarizedExperiment::assays(TSAM_Object)), function(x) colData(TSAM_Object)))

            allSampleData$Sample = newSamplesNames

            allSampleData <-  allSampleData[nonEmptySamples,]

            Obj1 <- SummarizedExperiment(
                assays = list('counts' = assayList),
                colData = allSampleData,
                rowRanges = rowRanges(TSAM_Object),
                metadata = TSAM_Object@metadata
            )
        }else{

            Obj1 <- TSAM_Object
            names(assays(Obj1)) <- 'counts'
            assays(Obj1)$counts[is.na(assays(Obj1)$counts)] = 0

        }

        CisbpAnno <- chromVAR::getAnnotations(TSAM_Object@metadata[[motifSetName]],
                                            rowRanges = rowRanges(Obj1))

        Obj1 <- chromVAR::addGCBias(Obj1, genome = metadata(TSAM_Object)$Genome)

        if(any(is.na(SummarizedExperiment::rowData(Obj1)$bias))){

            naList <- is.na(SummarizedExperiment::rowData(Obj1)$bias)
            message(paste(sum(naList), "NaNs found within GC Bias", sep =" "))

            SummarizedExperiment::rowData(Obj1)$bias[which(naList)] = mean(rowData(Obj1)$bias, na.rm = TRUE)

        }
                                    
        dev <- chromVAR::computeDeviations(object = Obj1, 
                                annotations = CisbpAnno)
        return(dev)

    }else{

        devList <- NULL

        motifList <- S4Vectors::metadata(TSAM_Object)[[motifSetName]]

        cellTypes <- names(SummarizedExperiment::assays(TSAM_Object))
        cellTypes <- cellTypes[cellTypes %in% colnames(GenomicRanges::mcols(SummarizedExperiment::rowRanges(TSAM_Object)))]

        for(i in cellTypes){

            message(paste('Running ChromVAR on ', i, sep =''))

            accMat <- MOCHA::getCellPopMatrix(TSAM_Object, cellPopulation = i, dropSamples = TRUE, NAtoZero = TRUE)
            if(dim(accMat)[2] <= 3){
                message('Three or less samples found. Skipping this cell type')
                
                devList <- append(devList, list(NULL))

            }else{

                sampleData <- SummarizedExperiment::colData(TSAM_Object)

                sampleData_filtered <- sampleData[rownames(sampleData) %in% colnames(accMat), ]

                subRanges <- MOCHA::StringsToGRanges(rownames(accMat))

                Obj1 <- SummarizedExperiment::SummarizedExperiment(
                    assays = list('counts' = accMat),
                    colData = sampleData_filtered,
                    rowRanges = subRanges,
                    metadata = STM@metadata
                )

                Anno <- chromVAR::getAnnotations(motifList,
                                                    rowRanges = subRanges)

                Obj1 <- chromVAR::addGCBias(Obj1, genome = S4Vectors::metadata(TSAM_Object)$Genome)
                if(any(is.na(SummarizedExperiment::rowData(Obj1)$bias))){

                    naList <- is.na(SummarizedExperiment::rowData(Obj1)$bias)
                    message(paste(sum(naList), "NaNs found within GC Bias", sep =" "))

                    SummarizedExperiment::rowData(Obj1)$bias[which(naList)] = mean(rowData(Obj1)$bias, na.rm = TRUE)

                }
                                            
                dev <- chromVAR::computeDeviations(object = Obj1, 
                                        annotations = Anno)

                devList <- append(devList, list(dev))

            }

        }

        names(devList) <- cellTypes
              
        return(devList)

    }


}



#' @title \code{subsetDev}
#'
#' @description \code{subsetDev} subsets a chromVAR output by a given metadata column (subsetBy) and the specific groups within (groupList)
#'
#'
#' @param Obj chromVAR output, from a MOCHA object. Should come from runChromVAR
#' @param subsetBy the variable to subset by. Can either be 'celltype', or a
#'   column from the sample metadata (see colData(Object))
#' @param groupList the list of cell type names or sample-associated data that
#'   should be used to subset the Object
#' @param na.rm removes groups that are NA if set to true. If set to false, then
#'   you filter for everything in the groupList and also NA values.
#'
#' @return Object the input Object, filtered down to either the cell type or
#'   samples desired.
#'
#' @details This is a wrapper for basic SummarizedExperiment-based subsetting. 
#'
#'
#' @export
#' 

subsetDev <- function(Object,
                              subsetBy,
                              groupList){

    sampleData <- SummarizedExperiment::colData(Object)

    #Check if subsetting by cell type
    if (tolower(subsetBy) == "celltypes") {
      if (!all(groupList %in% names(SummarizedExperiment::assays(Object)))) {
        stop("Error: groupList includes celltypes not found within Object.")
      }

      keep <- which(names(SummarizedExperiment::assays(Object)) %in% groupList)

      SummarizedExperiment::assays(Object) <- SummarizedExperiment::assays(Object)[keep]

      return(Object)
    }
    
    #check if subsetting by sample metadata
    if (subsetBy %in% colnames(sampleData)) {
        if (!all(groupList %in% unique(sampleData[[subsetBy]]))) {
        stop("Error: groupList includes names not found within the object sample data. Please check groupList.")
        }
    }

    keep <- rownames(sampleData)[which(sampleData[[subsetBy]] %in% groupList)]
    return(Object[, keep])
    
}

#' @title \code{modelDeviations}
#'
#' @description \code{modelDeviations} Runs a generalized linear model on ChromVAR deviation outputs. 
#'
#'
#' @param Obj chromVAR output, from a MOCHA object. Should come from runChromVAR
#' @param formula formula for the linear model
#' @param type Type of chromVAR output to model. Should be either 'deviations' or 'z' score.
#' @param verbose Set TRUE to display additional messages. Default is FALSE.
#' @param numCores 
#' 
#' @return lmem_res A list of linear model results for each chromVAR deviations. 
#'
#' @details relies on lmerTest to run a linear model on each chromVAR result. 
#'
#' @export
#' 

                  
modelDeviations <- function(Obj, formula, type = 'z', numCores = 1){
    

    meta1 <- as.data.frame(SummarizedExperiment::colData(Obj))

    if(type != 'z' & type != 'deviations'){
            stop('type must be with z or deviations')

    }
    mat1 <- SummarizedExperiment::assays(Obj)[[type]] 

    meta <- meta1[meta1$Sample %in% colnames(mat1),]

    cl <- parallel::makeCluster(numCores)

    parallel::clusterExport(cl=cl, varlist=c("meta", "formula", "mat1"), envir=environment())
    
    suppressMessages(lmem_res <- pbapply::pblapply(c(1:nrow(mat1)),
        function(x) {
            df <- data.frame(exp = as.numeric(mat1[x, ]), 
                meta, stringsAsFactors = FALSE)
            lmerTest::lmer(formula = formula, data = df)
        }, cl = cl), classes = "message")

    return(lmem_res)

    
}
              