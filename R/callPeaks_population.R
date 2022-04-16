#' @title \code{callPeaks}
#'
#' @description \code{callPeaks} is the main peak-calling function in scMACS
#'              that serves as a wrapper function to call peaks provided a set 
#'              of fragment files and an ArchR Project for meta-data purposes
#'
#'
#' @param ArchRProj an ArchR Project 
#' @param cellSubsets vector of strings. Cell subsets for which to call peaks. Optional, if cellSubsets='ALL', then peak calling is done on all cell populations in the ArchR project metadata
#' @param cellCol_label_name string indicating which column in the meta data file contains 
#'        the cell population label
#' 
#' @param returnAllPeaks boolean. Indicates whether scMACS should return object containing all genomic regions or just the positive (+) called peaks. Default to the latter, only positive peaks. 
#' 
#' @param numCores integer. Number of cores to parallelize peak-calling across
#'                 multiple cell populations 
#'
#' @param totalFrags # of fragments in that sample for that cell population
#' @return scMACs_PeakList an list containing peak calls for each cell population passed on in the 
#'         cell subsets argument. Each peak call is returned as as Genomic Ranges object.
#' 
#' @details The technical details of the algorithm are found in XX.
#'
#' @references XX
#'
#' @export

callPeaks_by_population <- function(ArchRProj, 
                      cellSubsets=NULL,
                      cellCol_label_name=NULL,
                      returnAllPeaks=FALSE,
                      numCores=10,
                      totalFrags,
                      fragsList=NULL
                     
                     ){
    
                               
    ########################################################################
    ########################################################################
    
    ### User-input/Parameter Checks
    
      if(class(cellSubsets)!='character'){
        stop('cellSubsets must be a vector of strings indicating cell populations')
      }
    
      if(class(cellCol_label_name)!='character'){
        stop('cellCol_label_name must be a vector of strings indicating cell populations')      
      }   

      if(class(ArchRProj)!='ArchRProject'){
        stop('ArchRProject must be an ArchR Project')
      } 
      if(is.null(fragsList)){
          stop('Load fragments prior to running scMACS')
          
      }
    

    
                               
    ########################################################################
    ########################################################################
    ## coefficients trained on ~ 3600 frags per cell 
    ## and future datasets need to be calibrated to
    ## these coefficients 
    finalModelObject = scMACS::finalModelObject
    fragsList <- S4Vectors::SimpleList(fragsList)
    thresholdModel = scMACS::thresholdModel
    
    ### obtain meta data from ArchR Project
    meta = getCellColData(ArchRProj)
    
    if(cellSubsets=='ALL'){
       cellPopulations= unique(meta[,cellCol_label_name])   
      
    } else{
       cellPopulations=cellSubsets
    }
    
    ### get barcodes by cell pop for 
    ### peak-calling by different 
    ### cell population 
    
    barcodes_by_cell_pop <- lapply(cellPopulations, function(x)
        row.names(meta)[which(meta[,cellCol_label_name]==x)]
           )
    
    ### create named list 
    names(barcodes_by_cell_pop) <- cellPopulations  

    cellsPerPop <- sapply(barcodes_by_cell_pop , function(x) length(x)
           )
    
    df <- cbind(cellPopulations, cellsPerPop)
        
    cat('Running peak calls for the following cell populations:')
    print(df)
    
    
    ### call peaks by cell-subsets

    callPeaks_by_cell_population <- function(fragsList, cellNames,ArchRProj, totalFrags,
                                            normScale=10^9){

        print(length(cellNames))

        ### subset ArchR Project
        cellSubsetArchR <- subsetCells(ArchRProj, cellNames=cellNames)
        
        subset_Frag <- function(cellNames, tmp){
            tmp_df <- GenomicRanges::as.data.frame(tmp)
            idx <- which(tmp_df$RG %in% cellNames)
            tmp[ idx]
           
        }
        fragsList_by_cell <- mclapply(fragsList, 
                                function(x) subset_Frag(cellNames, x)
                                )
        #print(fragsList_by_cell)
                                         

        fragsList_by_cell = SimpleList(fragsList_by_cell)
        fragMat <- unlist(fragsList_by_cell)    

        FinalBins <-  determine_dynamic_range(fragsList_by_cell,
                                         cellSubsetArchR, 
                                         binSize=500, 
                                         doBin=FALSE)

        countsMatrix <- calculate_intensities(fragMat, 
                                              FinalBins,
                                              totalFrags=totalFrags
                                            )
        
        countsMatrix = countsMatrix[countsMatrix$TotalIntensity > 0,]
        scMACS_peaks <- make_prediction(countsMatrix, finalModelObject,thresholdModel
                                       )

        if(returnAllPeaks){
            return(scMACS_peaks)

        } else {
            scMACS_peaks = scMACS_peaks[scMACS_peaks$peak==T]
            return(scMACS_peaks)

        }


    }
    
        
    scMACs_PeakList <- mclapply(1:length(barcodes_by_cell_pop), 
                                function(ZZ) callPeaks_by_cell_population(fragsList,                                                                    cellNames=barcodes_by_cell_pop[[ZZ]], 
                                                       ArchRProj,
                                                       totalFrags,
                                            normScale=10^9),
                                mc.cores =numCores,
                                mc.preschedule=TRUE,
                                mc.allow.recursive=FALSE
    ) 
    
    names(scMACs_PeakList) <- cellPopulations
                                              
                           
    ########################################################################
    ########################################################################


    return(scMACs_PeakList)

}
