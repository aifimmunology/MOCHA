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
#'                 multiple cell populationsdevto
#' 
#'
#' @return scMACs_PeakList an list containing peak calls for each cell population passed on in the 
#'         cell subsets argument. Each peak call is returned as as Genomic Ranges object.
#' 
#' @details The technical details of the algorithm are found in XX.
#'
#' @references XX
#'
#' @export

callPeaks <- function(ArchRProj, 
                      cellSubsets=NULL,
                      cellCol_label_name=NULL,
                      returnAllPeaks=FALSE,
                      numCores=10
                     
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
    

    
                               
    ########################################################################
    ########################################################################
    ## coefficients trained on 3600 frags per cell 
    ## and future datasets need to be calibrated to
    ## these coefficients 
    finalModelObject = scMACS::finalModelObject
    
    medianFrags_training = 3618
    
    ### load fragment files from ArchR Project 
    fragsList<-  getFragmentsFromProject(ArchRProj)
    
    ### obtain meta data from ArchR Project
    meta = getCellColData(ArchRProj)
    
    if(cellSubsets=='ALL'){
       cellPopulations= unique(meta[,cellCol_label_name])   
      
    } else{
       cellPopulations=cellSubsets
    }

    
    ### obtain median # of fragments
    ### per cell to calibrate features
    ### to pre-trained model 
    medianFrags_current = median(ArchRProj@cellColData$nFrags)
    
    ### identify scaling factor 
    scaleFactor= medianFrags_training/ medianFrags_current
    
    cat(paste('\nScale factor is ', round(scaleFactor,2),'\n\n'))

    
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

    callPeaks_by_population <- function(fragsList, cellNames,ArchRProj, scaleFactor){

        print(length(cellNames))

        ### subset ArchR Project
        cellSubsetArchR <- subsetCells(ArchRProj, cellNames=cellNames)

        print(cellSubsetArchR)
        
        subset_Frag <- function(cellNames, tmp){
            tmp_df <- GenomicRanges::as.data.frame(tmp)
            idx <- which(tmp_df$RG %in% cellNames)
            tmp[ idx]
           
        }
        fragsList_by_cell <- mclapply(fragsList, 
                                function(x) subset_Frag(cellNames, x)
                                )
        print(fragsList_by_cell)
                                         

        fragsList_by_cell = SimpleList(fragsList_by_cell)
        fragMat <- unlist(fragsList_by_cell)    

        FinalBins <-  determine_dynamic_range(fragsList_by_cell,
                                         cellSubsetArchR, 
                                         binSize=500, 
                                         doBin=FALSE)

        countsMatrix <- calculate_intensities(fragMat, 
                                             FinalBins,
                                             NULL,
                                             normalizeBins=FALSE,
                                             theta=0.001,
                                             width=20
                                            )
        
        countsMatrix = countsMatrix[countsMatrix$lambda1 > 0,]

        countsMatrix$lambda1 <- countsMatrix$lambda1 * scaleFactor
        countsMatrix$lambda2 <- countsMatrix$lambda2 * scaleFactor        

        scMACS_peaks <- make_prediction(countsMatrix, finalModelObject )

        if(returnAllPeaks){
            return(scMACS_peaks)

        } else {
            scMACS_peaks = scMACS_peaks[scMACS_peaks$Peak==T]
            return(scMACS_peaks)

        }


    }
    
        
    scMACs_PeakList <- mclapply(1:length(barcodes_by_cell_pop), 
                                function(ZZ) callPeaks_by_population(fragsList,                                                                           cellNames=barcodes_by_cell_pop[[ZZ]], 
                                                                            ArchRProj, 
                                                                            scaleFactor
                                                                            ),
                                mc.cores =numCores,
                                mc.preschedule=TRUE,
                                mc.allow.recursive=FALSE
    ) 
    
    names(scMACs_PeakList) <- cellPopulations
                                              
                           
    ########################################################################
    ########################################################################


    return(scMACs_PeakList)

}
