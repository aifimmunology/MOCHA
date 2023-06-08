#' @title \code{makeChromVAR}
#'
#' @description \code{makeChromVAR} pseodubulks a Seurat object by sample and cell type into a SummarizedExperiment, similar to MOCHA
#'
#' @param TSAM_Object A SummarizedExperiment object generated from
#'   getSampleTileMatrix or combineSampleTileMatrix
#' @param cellPopulation Names of cell types to analyze. Must match assay names in the SummarizedExperiment object. 
#'  Alternative, if you want to run ChromVAR across all celltypes, you can provide the output of combineSampleTileMatrix, and set this parameter to 'counts'.
#' @param motifName Name of metadata slot that has motif information for analysis.
#' @param exportList Boolean. Default is true, and will export a seperate ChromVAR object for each cell type in cellPopulation. If set to false,
#'          It will combine this list into one output SummarizedExperiment. 
#' @param 
#' @param verbose Boolean.
#' @return A named list of ChromVAR objects. 
#'
#'
#' @export
makeChromVAR <- function(TSAM_Object,
                      cellPopulation,
                      motifName,
                      exportList = TRUE,
                      verbose = TRUE) {

    if (
        any(!cellPopulation %in% names(SummarizedExperiment::assays(TSAM_Object)))
    ) {

        stop("cellPopulation was not found within TSAM_Object.")

    } else if (length(cellPopulation) > 1) {
   
        #Generate a list of combined objects
        newObj <- lapply(cellPopulation, function(x){
            combineSampleTileMatrix(subsetMOCHAObject(TSAM_Object, subsetBy = 'celltype', groupList = x, subsetPeaks = TRUE))
        })
        names(newObj) <- cellPopulation

    } else if(cellPopulation == 'counts'){
        newObj <- TSAM_Object
    }else{
        newObj <- combineSampleTileMatrix(subsetMOCHAObject(TSAM_Object, subsetBy = 'celltype', groupList = cellPopulation, subsetPeaks = TRUE))
    }
    genome <- S4Vectors::metadata(TSAM_Object)$Genome
    genome <- BSgenome::getBSgenome(genome)

    BiocParallel::register(BiocParallel::SerialParam())
    #Either iterate ovewr a list of newObj, or just directly on newObj to generate ChromVAR
    if(any(tolower(class(newObj)) %in% 'list')){

        chromVAROut <- lapply(cellPopulation, function(XX){
            if(verbose){ message('Analyzing ', XX)}
            newObj[[XX]] <- addGCBias_ChAI(newObj[[XX]], genome)
            anno_ix <- chromVAR::getAnnotations(newObj[[XX]]@metadata[[motifName]], 
                        rowRanges = SummarizedExperiment::rowRanges(newObj[[XX]]))
            chromVAR::computeDeviations(object =newObj[[XX]], 
                annotations = anno_ix)
        })

    }else{

        if(verbose){ message('Analyzing ', cellPopulation)}
        newObj[[XX]] <- addGCBias_ChAI(newObj[[XX]], genome)
        anno_ix <- chromVAR::getAnnotations(newObj@metadata[[motifName]], 
                        rowRanges = SummarizedExperiment::rowRanges(newObj))
        chromVAROut <- list(chromVAR::computeDeviations(object = newObj, 
                    annotations = anno_ix))

    }

    names(chromVAROut) <- cellPopulation
    if(exportList){
        return(chromVAROut)
    }

    newOut_Dev <- reformatChromVARList(chromVAROut, selectDev =TRUE)
    newOut_Z <- reformatChromVARList(chromVAROut, selectDev =FALSE)

    newOut <- list('Z_Score' = newOut_Z, 'Deviations' = newOut_Dev)
    return( newOut)

}

addGCBias_ChAI <- function(obj1, genome){

    obj1 <- addGCBias(obj1, genome = genome)
    if (any(is.na(SummarizedExperiment::rowData(obj1)$bias))) {
        naList <- is.na(SummarizedExperiment::rowData(obj1)$bias)
        
        if (verbose) {
        warning(paste(sum(naList), "NaNs found within GC Bias", sep = " "))
        }
        
        SummarizedExperiment::rowData(obj1)$bias[which(naList)] <- mean(rowData(obj1)$bias, na.rm = TRUE)
    }
    return(obj1)
}

#' @title \code{reformatChromVARList}
#'
#' @description \code{reformatChromVARList} pseodubulks a Seurat object by sample and cell type into a SummarizedExperiment, similar to MOCHA
#'
#' @param TSAM_Object A SummarizedExperiment object generated from
#'   getSampleTileMatrix or combineSampleTileMatrix
#' @param cellPopulation Names of cell types to analyze. Must match assay names in the SummarizedExperiment object. 
#'  Alternative, if you want to run ChromVAR across all celltypes, you can provide the output of combineSampleTileMatrix, and set this parameter to 'counts'.
#' @param motifName Name of metadata slot that has motif information for analysis.
#' @return A named list of ChromVAR objects. 
#'
#'
#' @export

reformatChromVARList <- function(chromVARList, selectDev = TRUE){

    if(!class(chromVARList) == 'list'){
        stop('chromVARList is not a list.')
    }else if(!all(unlist(lapply(chromVARList, class)) == 'chromVARDeviations')){
        stop('Some or all indices of chromVARList are not chromVarDeviations objects. Check your list.')
    }
    #browser()
    allMeta <- do.call('rbind', lapply(chromVARList, SummarizedExperiment::colData))
    #Identify the end of the original metadata, after which CellType, and various other CellType-Sample metadata was tacked on. 
    limitation1 <-which(colnames(allMeta) == 'CellType')
    #Transform metadata for cell type, and find the sample-level data
    newCellType <- paste(gsub(" |//.","_", unique(allMeta$CellType)), collapse ="|")
    newSample <- sub("__","",gsub(newCellType,"", allMeta$Sample))
    originMeta <- allMeta[, 1:(limitation1-1)]
    originMeta$Sample = newSample
    sampleData <- dplyr::distinct(as.data.frame(originMeta))
    rownames(sampleData) <- sampleData$Sample


    ## Transform the CellType-Sample metadata
    sumData <- allMeta[,(limitation1+1):length(colnames(allMeta))] %>% as.data.frame() 
    sumDatalist <- lapply(colnames(sumData), function(x){
            tmpMat <- sumData[,x, drop = FALSE]
            tmpMat$CellType = allMeta$CellType
            tmpMat$Sample = newSample
            newTmp <- tidyr::pivot_wider(tmpMat,id_cols = 'CellType', names_from = 'Sample',
                            values_from = {{x}}) %>% as.data.frame()
            rownames(newTmp) <- newTmp$CellType
            dplyr::select(newTmp, !CellType)
    })

    names(sumDatalist) <-   colnames(sumData)     

    summarizedData = SummarizedExperiment::SummarizedExperiment(sumDatalist, colData = sampleData)

    # Pull out the assays needed and rename them after the celltypes. 
    if(selectDev){
        assayType = 'deviations'
    }else{
        assayType = 'z'
    }
    assayList <- lapply(chromVARList, function(x){
        tmpMat <- SummarizedExperiment::assays(x)[[assayType]]
        colnames(tmpMat) <- sub("__","",gsub(newCellType,"", colnames(tmpMat)))
        tmpMat
        })

    outputSE <- SummarizedExperiment::SummarizedExperiment(assayList, colData = sampleData, 
                    metadata = list('SummarizedData' = summarizedData))

    return(outputSE)
}

