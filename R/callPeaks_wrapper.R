#' @title \code{callPeaks}
#'
#' @description \code{callPeaks} is the main peak-calling function in scMACS
#'              that serves as a wrapper function to call peaks provided a set 
#'              of fragment files and an ArchR Project for meta-data purposes
#'
#'
#' @param X an intensityMatrix output from \code{calculate_intensities}
#' @param finalModel is a matrix with coefficients and an index indicating
#'        the number of cell used to train that model
#'
#' @return the original intensityMatrix with the two intensity parameters required
#' to calculate the probability of a (+) peak, with an additional two columns
#' that include the prediction probability
#'
#' @details The technical details of the algorithm are found in XX.
#'
#' @references XX
#'
#' @export

callPeaks <- function(ArchRProj, 
                      cellSubsets='CD14 Mono', 
                      cellCol_label_name='predictedGroup_Col2.5'
                     
                     ){
    
                               
    ########################################################################
    ########################################################################
    
    
      if(class(cellSubsets)!='character'){
        stop('cellSubsets must be a vector of strings indicating cell populations')
      }
    
      if(class(cellCol_label_name)!='character'){
        stop('cellCol_label_name must be a vector of strings indicating cell populations')      
      }   

      if(class(ArchRProject)!='ArchRProject'){
        stop('ArchRProject must be an ArchR Project')
      }
    


    
    
                               
    ########################################################################
    ########################################################################
    medianFrags_training = 3618
    ### load fragment files from ArchR Project 
    fragsList<-  getFragmentsFromProject(ArchRProj)
    
    ### obtain meta data from ArchR Project
    meta = getCellColData(ArchRProj)
    
    ### obtain median # of fragments
    ### per cell to calibrate features
    ### to pre-trained model 
    medianFrags_current = median(ArchRProj@cellColData$nFrags)
    
    ### identify scaling factor 
    scaleFactor= medianFrags_training/ medianFrags_current
    
    ### get barcodes by cell pop for 
    ### peak-calling by different 
    ### cell population 
    
    barcodes_by_cell_pop <- lapply(cellSubsets, function(x)
        row.names(meta)[which(meta[,cellCol_label_name]==x)]
           )
    
    ### create named list 
    names(barcodes_by_cell_pop) <- cellSubsets  
    
    ### call peaks by cell-subsets
    
    callPeaks_by_population <- function(fragsList, cellSubset){
        
        ### subset ArchR Project
        cellSubsetArchR <- subsetCells(ArchRProj,barcodes_by_cell_pop[[cellSubset]])
        metaSub= getCellColData(cellSubsetArchR)
         
        fragsList_by_cell <- SimpleList(lapply(fragsList, 
                                function(x) x[x$RG %in% barcodes_by_cell_pop[[cellSubset]] ]
                                )
                                         )
        
        fragMat <- unlist(fragsList_by_cell)    
        
        FinalBins <-  determine_dynamic_range(fragsList_by_cell,
                                         ArchRProj, 
                                         binSize, 
                                         doBin=FALSE)
        
        countsMatrix <- calculate_intensities(fragMat, 
                                             FinalBins,
                                             NULL,
                                             normalizeBins=FALSE,
                                             theta=0.001,
                                             width=20
                                            )
        
        countsMatrix$lambda1 <- countsMatrix$lambda1 * scaleFactor
        countsMatrix$lambda2 <- countsMatrix$lambda2 * scaleFactor        
        
        scMACS_peaks <- scMACS::make_prediction(countsMatrix, finalModelObject )
        
        return(scMACS_peaks)
        
    }
    
    scMACs_PeakList <- mclapply(cellSubsets, function(x)      
        callPeaks_by_population(fragsList, x),
                                mc.cores =10
                                )
    
    names(scMACs_PeakList) <- cellSubsets
                                              
                           
    ########################################################################
    ########################################################################


    return(scMACs_PeakList)

}
