#' @title \code{getGeneMotifLinks}
#'
#' @description \code{getGeneMotifLinks} takes aset of genes, pseudobulk expression and accessibility, and finds region accessibility to gene expression correlations.
#'
#' @param ChromVARObj The SummarizedExperiment object output from getSampleTileMatrix containing your sample-tile matrices
#' @param SampleTileObj The SummarizedExperiment object output from getSampleTileMatrix containing your sample-tile matrices
#' @param cellPopulation A string denoting the cell population of interest, which must be present in SampleTileObj
#' @param DEGList A list of gene names to test. Must aligned with gene names from the SampleGeneObj. Should be differential.
#' @param DUMList a list of differential ChromVAR Z-scores to test. 
#' @param windowSize the size of the window, in basepairs, around each input region to search for co-accessible links
#' @param backgroundSize Integer. the size of background tileset to use. If NULL, the background tileset will be sized matched to the foreground set.
#' @param returnBackground Boolean. Set TRUE to return a list of both foreground and background correlations. 
#' @param numCores Optional, the number of cores to use with multiprocessing. Default is 1.
#' @param verbose Set TRUE to display additional messages. Default is FALSE.

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
getGeneMotifLinks <- function(ChromVARObj,
                          SampleGeneObj,
                          DUMList = NULL,
                          DEGList = NULL,
                          returnBackground = FALSE,
                          numCores = 20,
                          verbose = TRUE) {

    
    . <- NULL


    if(length(SummarizedExperiment::assays(SampleGeneObj)) > 1){
        stop("SampleGeneObj needs to be subsetted to only one cell type.")
    }

    if(!all(colnames(ChromVARObj) %in% colnames(SampleGeneObj))){
        stop("Samples are not aligned between ChromVARObj and SampleGeneObj.")
    }

    if(!all(DEGList %in% rownames(SampleGeneObj))){
        stop("Not all gene names are found in SampleGeneObj.")
    }

    if(!all(DUMList %in% rownames(ChromVARObj))){
        stop("Motif names not found in ChromVARObj.")
    }

    if(length(DUMList)*2 > length(rownames(ChromVARObj))){
        stop("More than half of motifs are listed as Differentially Used. This is unlikely to be accurate. Provide a shorter DUMList.")
    }

    if(length(DEGList)*2 > length(rownames(SampleGeneObj))){
        stop("More than half of genes are listed as Differentially Expressed. This is unlikely to be accurate. Provide a shorter DEGList.")
    }

    ChromVARObj

    linkOutput <- getGeneralLinks(Obj1 = ChromVARObj, 
                    Obj2 = SampleGeneObj
                    Obj1List = DUMList,
                    Obj2List = DEGList,
                    returnBackground = returnBackground,
                    ZI_inflated = FALSE,
                    numCores = numCores,
                    verbose = verbose)

    if(returnBackground){

        colnames(linkOutput$Foreground) <- c('Motifs','Genes','Correlations','pValues', 'FDR')
        colnames(linkOutput$Background) <- c('Motifs','Genes','Correlations')
        return(linkOutput)
        
    }else{

        colnames(linkOutput$Background) <- c('Motifs','Genes','Correlations','pValues', 'FDR')
        return(linkOutput)

    }
 
}
