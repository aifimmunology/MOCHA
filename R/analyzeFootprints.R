#' @title plotMotifs - plots motif footprints, exported from motifFootprint
#'
#' @description Generates plot 
#'
#' @param motifSE A SummarizedExperiment with motif footprinting information, from motifFootprint
#' @param footprint Optional string, to describe which footprint within the motifSE to analyze. Default is NULL, at which point it will pull from all footprint present. 
#' @param groupColumn An optional string, that will contain the group-level labels for either calling footprints as significant, or for comparing motif footprintings stats betweem groups. Default is null, at which point no wilcoxon tests will be conducted. groupColumn can be a sample-level group from the MOCHA's object colData slot, 'CellType' (if you want to compare across cell types), or 'Motif' if you want to compare across motifs. 
#' @param returnDF Boolean, default is FALSE, determines whether or not to return a data.frame, rather than plotting. 
#' @param returnPlotList A boolean, default is false, determines whether to return the full plot, or a list of subplots (ggplot2-based) for custom arrangements. 
#' @param topPercentages. A number 1 or less, that describes the top percent of motif locations to use for plotting. If 0.9, then the top 90% of regions, by total insertions, will be used in the plot. If 0.1, then the top 10% of regions will be used. Default of 1 uses all regions.
#' @param plotIndividualRegions A boolean, default is TRUE, that determines whether to plot the individual motif regions. If FALSE, only the motif footprint will be returned. 
#' @param relHeights A vector of two numbers, describing the relative space in the plot given to the motif summary and individual location plot. 
#'
#' @return a data.frame, containing motif footprint stats. 
#'
#' @examples
#' \dontrun{
#' p1 <- MOCHA::plotMotifs(
#'   motifSE,
#'   footprint = "CD4_Naive_ARID5A", groupColumn = 'COVID_status', topPercentages = 0.1)
#' }
#'
#' @export
#' @keywords exporting

plotMotifs <- function(motifSE, 
                       footprint = NULL, 
                       groupColumn = NULL, 
                       returnDF = FALSE, 
                       returnPlotList = FALSE,
                       topPercentage = 1,
                       plotIndividualRegions = TRUE,
                          relHeights = c(0.3, 0.7)){

    if(is.null(footprint)){
    
        footprint = names(SummarizedExperiment::assays(motifSE))
        
    }
    if(!any(names(SummarizedExperiment::assays(motifSE)) %in% footprint)){
    
        stop('footprint name was not found within the motifSE object')
        
    }
    
    colData1 = SummarizedExperiment::colData(motifSE)

    if(!methods::is(topPercentage, 'numeric')){

        stop('topPercentage must be a number')
        
    }else if(topPercentage <= 0 | topPercentage > 1){


        stop('topPercentage must be greater than 0 and less then or equal to 1')
        
    }
    


    ## Pull out footprints across regions before running stats and add the direction feature. 
    if(length(footprint) > 1){
        
        motif_spec <- data.table::rbindlist(lapply(footprint, function(XX){
                    if(topPercentage < 1){
                        regionDistribution = Matrix::rowSums(motifSE@ExperimentList[[XX]])
                        motifSE@ExperimentList[[XX]] = motifSE@ExperimentList[[XX]][regionDistribution > 
                                                                    quantile(regionDistribution, probs = 1-topPercentage),]
                    }
                    tmpMat = Matrix::colMeans(SummarizedExperiment::assays(motifSE)[[XX]])
                    Sample= gsub('__-[0-9].*|__[0-9].*', '' ,names(tmpMat))
                    Position = gsub('__', '', gsub(paste0(unique(Sample), collapse='|'), '', names(tmpMat)))
                    data.table::data.table(Insertions = tmpMat, Sample = Sample,
                                           Position = as.numeric(Position), Footprint = XX)
            }))
        ## Process single location data as well
        if(plotIndividualRegions){
            dt1 = data.table::rbindlist(lapply(footprint, function(XX){
                        tmpMat = SummarizedExperiment::assays(motifSE)[[XX]]
                        tmpMat2 = as.data.table(Matrix::as.matrix(tmpMat))
                        tmpMat2$Location = rownames(tmpMat)
                        dt2 = data.table::melt(tmpMat2,
                                            id.vars = 'Location',
                                            measureVar = colnames(tmpMat),
                                               variable.name = 'Index',
                                               value.name = 'Insertions')
                        dt2$Footprint = XX
                        
                        rm(tmpMat2)
                        dt2                    
                }))
        }
        
    }else{
        if(topPercentage < 1){
            regionDistribution = Matrix::rowSums(motifSE@ExperimentList[[footprint]])
            motifSE@ExperimentList[[footprint]] = motifSE@ExperimentList[[footprint]][regionDistribution > 
                                                       quantile(regionDistribution, probs = 1-topPercentage),]
        }
        tmpMat = Matrix::colMeans(SummarizedExperiment::assays(motifSE)[[footprint]])
        Sample= gsub('__-[0-9].*|__[0-9].*', '' ,names(tmpMat))
        Position = as.numeric(gsub('__', '', gsub(paste0(unique(Sample), collapse='|'), '', names(tmpMat))))
        motif_spec <- data.table::data.table(Insertions = tmpMat, Sample = Sample, 
                                             Position = Position, Footprint = footprint)
        ## Process single location data as well        
        if(plotIndividualRegions){
            tmpMat = SummarizedExperiment::assays(motifSE)[[footprint]]
            tmpMat2 = as.data.table(Matrix::as.matrix(tmpMat))
            tmpMat2$Location = rownames(tmpMat)
            dt1 = data.table::melt(tmpMat2,
                                id.vars = 'Location',
                                measureVar = colnames(tmpMat),
                                   variable.name = 'Index',
                                   value.name = 'Insertions')
            dt1$Footprint = footprint
            rm(tmpMat2)
        }
    }    
    
    ## Split out Cell type and motif info
    motif_spec <- motif_spec[, ':='(CellType = gsub('__.*','', Footprint), 
                                    Motif = gsub('.*__','', Footprint))]

    if(plotIndividualRegions){
        dt1 <- dt1[, ':='(CellType = gsub('__.*','', Footprint), 
                                        Motif = gsub('.*__','', Footprint))]
        ## Now pull out Sample vs Position info
        dt1 <- dt1[, Sample := gsub('__-[0-9].*|__[0-9].*', '' , Index)]
        allSamples = paste0(unique(dt1$Sample), collapse='|')
        dt1 <- dt1[, Position := as.numeric(gsub('__', '', 
                                gsub(allSamples, '', Index)))]

    }else{

        dt1 = NULL
    }
    
     ### Add in sample metadata, if motif footprinting done on the sample level. 
    if(any(names(motifSE@metadata) == 'metadata')){
    
        metadf = as.data.table(motifSE@metadata[['metadata']])
        
        if(any(metadf$Sample %in% motif_spec$Sample)){
            motifDT <- motif_spec[metadf, on = .(Sample)]
            if(! is.null(dt1)){ singleDT <- dt1[metadf, on = .(Sample)] }
            
        }else{
            motifDT = motif_spec
            if(! is.null(dt1)){ 
                singleDT <- dt1 
            }else{
               singleDT = NULL
            }
        }
    }else{

        motifDT = motif_spec
       if(! is.null(dt1)){  
           singleDT <- dt1 
       }else{
           singleDT = NULL
        }
    }
    rm(dt1)
    
    if(!is.null(groupColumn)){

        ## Check that groupColumn only has one value and that it is present in the motifDT metadata
        if(length(groupColumn) != 1){

            stop('More than one value provided to groupColumn.')

        }

        if(!any(colnames(motifDT) %in% groupColumn)){

            stop('Group column not found in sample metadata')
        }

        if(!any(colnames(motifDT) %in% groupColumn)){

            stop('Group column not found in sample metadata')
        }

        #Setting grouping variable
        group_by_cols = unique(c('CellType', 'Motif', groupColumn, 'Position', 'Footprint'))
    
        ## Now conduct grouping
        motifDT = motifDT[, .(Insertions = mean(Insertions, na.rm = TRUE)),
                                     by =  group_by_cols]
        motifDT$PlotGroup = paste(groupColumn, motifDT[,get(groupColumn)], motifDT$Footprint, sep = '_')
        if(! is.null(singleDT)){ 
            singleDT = singleDT[, .(Insertions = mean(Insertions, na.rm = TRUE)),
                                         by =  c('Location', group_by_cols)]
            singleDT$PlotGroup = paste(groupColumn, singleDT[,get(groupColumn)], singleDT$Footprint, sep = '_')
        }
        
    }else{
        motifDT$PlotGroup = paste(motifDT[,Sample], motifDT$Footprint, sep = '_')
        
        if(! is.null(singleDT)){  singleDT$PlotGroup = paste(singleDT[,Sample], singleDT$Footprint, sep = '_') }
    }

    if(! is.null(singleDT)){ 
        ##Order the singleDT by average strength across conditions
        avgStrength = singleDT[, .(totalInsertions = sum(Insertions, na.rm = TRUE)), by = .(Location, PlotGroup)]
        avgStrength = avgStrength[order(totalInsertions,decreasing=TRUE )]
        singleDT$PlotGroup = factor(singleDT$PlotGroup, levels = rev(unique(avgStrength$PlotGroup)))
        singleDT$Location = factor(singleDT$Location, levels = rev(unique(avgStrength$Location)))
    }
    
    if(returnDF){
        
      returnList = list('MotifAverage' = motifDT, 'LocationSpecific' = singleDT)
        return(returnList)
    }
    
    maxPos = max(motifDT$Position)
    
    p1 = ggplot2::ggplot(motifDT, ggplot2::aes(x = Position, y = Insertions,
                                               group = PlotGroup, color= PlotGroup)) + 
            ggplot2::geom_line() + ggplot2::xlab(NULL) + 
             ggplot2::theme_bw() + ggplot2::theme(legend.position = 'top')+
            ggplot2::theme(plot.margin = ggplot2::unit(c(0, 0.1, 0, 0.1), "cm")) +
            ggplot2::xlim(-maxPos, maxPos) + 
            ggplot2::scale_x_continuous(expand=c(0,0)) 

    if(plotIndividualRegions){
        if (!requireNamespace("ggrastr", quietly = TRUE)) {
          stop(
          "Package 'ggrastr' is required for plotting individual locations.",
          "Please install ggrastr to proceed, or set plotIndividualRegions to FALSE."
          )
        }
    
        p2 = ggplot2::ggplot(singleDT, ggplot2::aes(x = Position, fill = Insertions, y = Location)) + 
                ggrastr::geom_tile_rast() +
                ggplot2::xlim( -maxPos, maxPos) +
                ggplot2::theme_bw() + 
                ggplot2::scale_x_continuous(expand=c(0,0)) + 
                ggplot2::facet_wrap(~ PlotGroup, ncol = 1 ) + #, strip.position = 'left') + 
                ggplot2::scale_fill_gradient(transform = 'log1p', low = 'white', high = 'red') + 
                ggplot2::theme(legend.position = 'right', 
                               axis.text.y = ggplot2::element_blank()) +
                ggplot2::theme(strip.background = ggplot2::element_rect(fill="white"),
                              strip.placement = "outside",
                             strip.text.y.left = ggplot2::element_text(angle = 90, vjust = 1),
                              plot.background = ggplot2::element_rect(fill = 'white',
                                                                      colour = 'white'))+
                ggplot2::guides(r = ggplot2::guide_axis(angle = 45)) +
                ggplot2::theme(plot.margin = ggplot2::unit(c(0, 0.1, 0.1, 0.1), "cm"))
        
       
        
        if(returnPlotList){
             returnList = list('MotifAverage' =  p1, 'LocationSpecific' = p2)
            return(returnList)
        }   
         returnList = list('MotifAverage' =  p1 + 
                    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                             axis.ticks.x = ggplot2::element_blank()), 
                'LocationSpecific' = p2)
       

        cp1 = cowplot::plot_grid(plotlist = returnList,
                                 ncol=1, align = 'v',
                                 axis = "brl",
                                rel_heights = relHeights)
        
    }else{

        if(returnPlotList){
            return(list('MotifAverage' =  p1, 'LocationSpecific' = NULL))
        }   
        return(p1)
        
    }
    
         
    return(cp1)  
    
}