#' @title \code{callPeaks_by_sample}
#'
#' @description \code{callPeaks_by_sample} is the main peak-calling function in scMACS
#'              that serves as a wrapper function to call peaks provided a set 
#'              of fragment files and an ArchR Project for meta-data purposes
#'
#'
#' @param ArchRProj an ArchR Project 
#' @param cellSubsets vector of strings. Cell subsets for which to call peaks. Optional, if cellSubsets='ALL', then peak calling is done on all cell populations in the ArchR project metadata
#' @param cellCol_label_name string indicating which column in the meta data file contains 
#'        the cell population label
#' @param sample_label_name string indicating which column in the meta data file contains 
#'        the samples
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

callPeaks_by_sample <- function(ArchRProj, 
                                meta,
                                cellSubsets=NULL,
                                cellCol_label_name=NULL,
                                sampleCol_label_name=NULL,         
                                batchCol_label_name=NULL,
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
    
    ## calculate internal scaling factors 
    meta$BatchColumn = meta[,batchCol_label_name]
    meta_dt <- as.data.table(meta)
    scaling_factors <- meta_dt[, list(ScalingFactor=median(nFrags)), by=BatchColumn]

    
    ## override in case you need to down-sample
    ## per MP 
    meta <- getCellColData(ArchRProj)
                               
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
    
    if(cellSubsets=='ALL'){
       cellPopulations= unique(meta[,cellCol_label_name])   
      
    } else{
       cellPopulations=cellSubsets
    }

    unique_samples <- unique(meta[,sampleCol_label_name])
    
    ### obtain median # of fragments
    ### per cell to calibrate features
    ### to pre-trained model 
#     medianFrags_current = median(ArchRProj@cellColData$nFrags)
    ### identify scaling factor 
#     scaleFactor= medianFrags_training/ medianFrags_current


    
        
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
    
    df_cell <- cbind(cellPopulations, cellsPerPop)
        
    cat('Running peak calls for the following cell populations:')
    print(df_cell)
    
    cat('Running peak calls for the following samples:')
    df_sample <- as.data.frame(table(meta[,sampleCol_label_name]))
    colnames(df_sample) <- c('Sample', '# Cells')
    print(df_sample)
    
    ### call peaks by cell-subsets

    callPeaks_by_population_sample <- function(fragsList, 
                                               cellNames,
                                               ArchRProj, 
                                               scaling_factors,
                                              batchCol_label_name){

        #print(length(cellNames))

        ### subset ArchR Project
        cellSubsetArchR <- subsetCells(ArchRProj, cellNames=cellNames)
        metaSubset <- getCellColData(cellSubsetArchR)
    

        print(cellSubsetArchR)
        
        
        ### Subset fragments 
        ### by cell population 
        subset_Frag <- function(cellNames, tmp){
            tmp_df <- GenomicRanges::as.data.frame(tmp)
            idx <- which(tmp_df$RG %in% cellNames)
            tmp[ idx]
           
        }
        
        ### obtain 
        fragsList_by_cell <- mclapply(fragsList, 
                                function(x) subset_Frag(cellNames, x)
                                )                                       

        fragsList_by_cell = SimpleList(fragsList_by_cell)
        
        
        fragMat <- unlist(fragsList_by_cell)    

        FinalBins <-  determine_dynamic_range(fragsList_by_cell,
                                         cellSubsetArchR, 
                                         binSize=500, 
                                         doBin=FALSE)
        
        ### Subset fragments 
        ### by samples 
        
        fragMat_df <- GenomicRanges::as.data.frame(fragMat, row.names=1:length(fragMat))
        
        subset_Frag_sample <- function(cell_sample, fragMat_df, fragMat){          
           idx <- which(fragMat_df$RG %in% cell_sample)
           fragMat[idx] 
        }
           
        
        
        barcodes_by_sample <- lapply(unique_samples,
                                     function(x) 
                            row.names(metaSubset)[which(metaSubset[,sampleCol_label_name] %in% x)]
                                     )
        

        
        
        frags_by_sample <- lapply(barcodes_by_sample,
                                  function(x) subset_Frag_sample(x,fragMat_df, fragMat)
                                  )
        
        

        countsMatrix_List <- mclapply(1:length(frags_by_sample), 
                                      function(x) calculate_intensities(frags_by_sample[[x]], 
                                                         FinalBins,
                                                         NULL,
                                                         normalizeBins=FALSE,
                                                         theta=0.001,
                                                         width=20
                                                        ),
                                      mc.cores=3
                            )
        
        countsMatrix_List = lapply(countsMatrix_List, 
                                   function(x) x[x$lambda1 > 0,]
                                   )
        
        scale_matrices <- function(tmp, meta, sampleID){
            
            sampleBatch <- meta[which(meta[,sampleCol_label_name]==sampleID)[1], batchCol_label_name]
            
            scaling_factors=as.data.frame(scaling_factors)
            medianFrags_current <- scaling_factors[scaling_factors$BatchColumn == sampleBatch,2]
            
            scaleFactor= medianFrags_training/ medianFrags_current
                
            tmp$lambda1 = tmp$lambda1 * scaleFactor
            tmp$lambda2 = tmp$lambda2 * scaleFactor

            return(tmp)
        }
        
        ## identify which sample 
        ## falls in which batch
        ## 

        countsMatrix_List = lapply(1:length(countsMatrix_List),
                                   function(x) scale_matrices(tmp=countsMatrix_List[[x]],
                                                              meta=meta,
                                                              sampleID = unique_samples[x]
                                                             )
                                   )
        
        idx <-  sapply(countsMatrix_List, function(x) x$numCells[1]
                       )
        nonzero_samples <- which(idx> 5)
        
        scMACS_peaks <- lapply(countsMatrix_List[nonzero_samples],
                               function(x) make_prediction(x, finalModelObject )
                               )
        
        names(scMACS_peaks) <- unique_samples[nonzero_samples]

                                         
        if(returnAllPeaks){
          return(list(scMACS_peaks=scMACS_peaks,
                        SampleFragments=frags_by_sample
                       ))

        } else {
            scMACS_peaks <- lapply( scMACS_peaks,
                                function(x) x[x$Peak==T]
                                )
            
            return(list(scMACS_peaks=scMACS_peaks,
                        SampleFragments=frags_by_sample
                       ))

        }


    }
    
        
    scMACs_PeakList <- mclapply(1:length(barcodes_by_cell_pop), 
                                function(ZZ) callPeaks_by_population_sample(fragsList,                                                                           cellNames=barcodes_by_cell_pop[[ZZ]], 
                                                                            ArchRProj, 
                                                                            scaling_factors,
                                                                            batchCol_label_name
                                                                            ),
                                mc.cores =numCores-1,
                                mc.preschedule=TRUE,
                                mc.allow.recursive=FALSE
    ) 
    
    names(scMACs_PeakList) <- cellPopulations
           
    
    
    #### take all the individual samples
    #### and take the union to identify
    #### all possible peak combinations
    #### across individual samples
    
    create_union_peakset <- function(cellType_peaks_by_sample){
    
        NSamples = length(cellType_peaks_by_sample)
        
        if(length(NSamples)>0){
            ## Take Union of Peakset 
            combinedPeakset <- as.data.table(do.call(rbind, cellType_peaks_by_sample)) 

            summarized_union <- combinedPeakset[, list(AverageLambda1=mean(lambda1), 
                                   AverageLambda2=mean(lambda2),
                                   SamplesWithPeak = .N,
                                   PropSamplesWithPeak = .N/NSamples), 
                            by=c('PeakID','seqnames','start','end','strand')]


            cellType_peaks_by_sample[['Union']] <- summarized_union
            return(cellType_peaks_by_sample)
        } else {
            print('No peaks were called')
            return(NULL)
            
        }
    }
    
    scMACs_PeakList2<- mclapply(1:length(scMACs_PeakList),
                                    function(x) try(create_union_peakset(scMACs_PeakList[[x]][[1]])),
                                      mc.cores=numCores,
                                      mc.preschedule=TRUE,
                                      mc.allow.recursive=FALSE
                                    )
    
    names(scMACs_PeakList2) <- cellPopulations
    
    for(i in 1:length(scMACs_PeakList)){
        
         scMACs_PeakList[[i]][['scMACS_peaks']] <- scMACs_PeakList2[[i]]

        
    }
                               
                           
    ########################################################################
    ########################################################################


    return(scMACs_PeakList)

}
