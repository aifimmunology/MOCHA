#' @title Generate motif footprint metrics and/or statistics across samples, cell types, and motifs
#'
#' @description Identifies the width and depth of a given motif for a given group & cell type. Requires the output of motifFootprint. If a groupColumn and background is provided, the function will also conduct wilcoxon tests on the motif baseline, depth, and width between groups. 
#'
#' @param motifSE A SummarizedExperiment with motif footprinting information, from motifFootprint
#' @param footprint Optional string, to describe which footprint within the motifSE to analyze. Default is NULL, at which point it will pull from all footprint present. 
#' @param callFootprint A boolean, to determine whether to call footprints as significant across samples. Requires sample-specific motif footprinting.
#' @param motifSize A number describing the approximate motif size you want to use, when looking at the insertion rate over the motif core. Default is 6 (3 bp up and 3 bp downstream)
#' @param searchWidth A number, describing the maximum motif footprint width. The function will look for the maximum insertions within +/- this number of the motif center. Default is 50.  
#' @param groupColumn An optional string, that will contain the group-level labels for either calling footprints as significant, or for comparing motif footprintings stats betweem groups. Default is null, at which point no wilcoxon tests will be conducted. groupColumn can be a sample-level group from the MOCHA's object colData slot, 'CellType' (if you want to compare across cell types), or 'Motif' if you want to compare across motifs. 
#' @param background An optional string if you want to compare two different conditions. This string is the background group that you want to use for running wilcoxon tests. 
#'
#' @return a data.frame, containing motif footprint stats. 
#'
#' @examples
#' \dontrun{
#' summaryStats <- MOCHA::motifStats(
#'   motifSE,
#'   footprint = "CD4 Naive_ARID5A", groupColumn = 'COVID_status', background = 'Negative')
#' }
#'
#' @export
#' @keywords exporting

motifStats <- function(motifSE, footprint = NULL, 
                       callFootprint = FALSE, groupColumn = NULL,
                       background = NULL, motifSize = 6, searchWidth = 50){
    
        browser()
    
    if(is.null(footprint)){
    
        footprint = names(SummarizedExperiment::assays(motifSE))
        
    }
    if(!any(names(SummarizedExperiment::assays(motifSE)) %in% footprint)){
    
        stop('footprint name was not found within the motifSE object')
        
    }
    if(!methods::is(motifSize, 'numeric') | !methods::is(searchWidth, 'numeric')){
    
        stop('motifSize and searchWidth must be numeric')
        
    }
    
    colData1 = SummarizedExperiment::colData(motifSE)

    ## Pull out footprints across regions before running stats and add the direction feature. 
    if(length(footprint) > 1){
        
        motif_spec <- do.call('rbind', lapply(footprint, function(XX){
                    tmpMat = colMeans(SummarizedExperiment::assays(motifSE)[[XX]])
                    Sample= gsub('__-[0-9].*|__[0-9].*', '' ,names(tmpMat))
                    Position = gsub('__', '', gsub(paste0(unique(Sample), collapse='|'), '', names(tmpMat)))
                    data.table::data.table(Insertions = tmpMat, Sample = Sample, Position = as.numeric(Position), Footprint = XX)
            }))
        
    }else{
        
        tmpMat = colMeans(SummarizedExperiment::assays(motifSE)[[footprint]])
        Sample= gsub('__-[0-9].*|__[0-9].*', '' ,names(tmpMat))
        Position = gsub('__', '', gsub(paste0(unique(Sample), collapse='|'), '', names(tmpMat)))
        motif_spec <- data.table::data.table(Insertions = tmpMat, Sample = Sample, Position = Position, Footprint = footprint)
    }    
    
    motif_spec <- dplyr::mutate( motif_spec,
                             Direction = ifelse(Position > 0, 'Downstream', 'Upstream'))
    
    ## Focus on the 50 basepairs up and downstream to identify the footprint size and depth
    ## Most footprints are less than 20 bp wide anyways
    local_motif_spec <- dplyr::mutate(motif_spec[abs(motif_spec$Position) <= searchWidth,], 
                                       pos = abs(Position), pos_squared = abs(pos)^2)
    ## Now for each sample/group and Footprint (CellType + Motif combination), split it into a seperate data.frame for 
    ## the left and right hand side of the motif footprint. 
    ## Then fit a quadratic model to it, and use that quadratic to fit to identify whether there is evidence for a footprint, 
    ## as well as what the height and radius of the motif is on that side.  
    halfSplits <- dplyr::group_split(local_motif_spec, Footprint, Sample, Direction)
    
    peakSum <- do.call('rbind',lapply(halfSplits, function(hs){
            tryCatch({
                initFit = summary(stats::lm(data = hs, formula = Insertions ~ pos + pos_squared))
                quadFit = as.data.frame(initFit$coefficients)
                Height = quadFit$Estimate[2]*peakWidth + 
                            quadFit$Estimate[3]*peakWidth^2 +
                            quadFit$Estimate[1]
                Radius = -quadFit$Estimate[2]/(2*quadFit$Estimate[3])
                R_squared = initFit$r.squared
                Footprint_pvalue = pf(initFit$fstatistic[1], initFit$fstatistic[2], 
                                      initFit$fstatistic[3], lower.tail = FALSE)
                ## If the y-intercept of the model is below the height at the vertex, 
                ## then the shape is downwards, not upwards.
                ## In that case, it's not a footprint, so set p-value to 1 and R squared statistic to 0
                if(quadFit$Estimate[1] > Height){
                    Footprint_pvalue = 1
                    R_squared = 0
                }
                data.table(
                    Sample = unique(hs$Sample),
                    Footprint = unique(hs$Footprint),
                    Direction = unique(hs$Direction),
                    Height = Height,
                    Radius = Radius,
                    R_squared = R_squared,
                    PValue = Footprint_pvalue
                )
            }, error = function(e){
                data.table(
                    Sample = unique(hs$Sample),
                    Footprint = unique(hs$Footprint),
                    Direction = unique(hs$Direction),
                    Radius = NA,
                    Height = NA,
                    R_squared = NA,
                    PValue = NA
                )
            })
        }))

    
    ### Remaining tests:
    ##### Test for normality for heigh and Radius (left and Right)
    ##### Join baseline to existing plot
    ##### Run statistical tests (left and right) for baseline vs Height. Get joint p-value
    
    ## Identify the baseline insertions around the motif center 
    averageMin =  as.data.table(motif_spec[abs(motif_spec$Position) < motifSize/2,])[ ,
                        .(Baseline = mean(Insertions, na.rm = TRUE)), by=.(Sample, Footprint)]
    ## Join the baseline with the other peak metrics. 
    peakSum <- peakSum[averageMin, on = .(Sample, Footprint),]
    
    ## Generate the estimated depth
    peakSum <- peakSum[, Depth := Height - Baseline]
    summaryStats <- peakSum[, ':='(CellType = gsub('__.*','', Footprint), 
                                    Motif = gsub('.*__','', Footprint))]
    
    ### Add in sample metadata, if motif footprinting done on the sample level. 
    if(any(names(motifSE@metadata) == 'metadata')){
    
        metadf = as.data.table(motifSE@metadata[['metadata']])
        
        if(any(metadf$Sample %in% summaryStats$Samples)){
            summaryStats <- summaryStats[metadf, on = .(Sample)]
        }
    }
    
    
    ## Check if the user wants to run a statistic tests on the motif footprints across samples. 
    if(!is.null(groupColumn)){
        
        ## Check that groupColumn only has one value and that it is present in the summaryStats metadata
        if(length(groupColumn) != 1){
            
            stop('More than one value provided to groupColumn.')
            
        }
        
        if(!any(colnames(summaryStats) %in% groupColumn)){
        
            stop('Group column not found in sample metadata')
        }
        
        if(!any(colnames(summaryStats) %in% groupColumn)){
        
            stop('Group column not found in sample metadata')
        }
        
        ## First check if the user wants a statistics between motifs or between cell types
        group_by_cols = unique(c('CellType', 'Motif', 'Direction', groupColumn))
        if(any(groupColumn %in% c('CellType', 'Motif'))){

                group_by_cols = group_by_cols[!group_by_cols %in% groupColumn ]

        }
        
        if(is.null(background)){
        
            
            ### If no background is defined, just generate group averages
            ## Width, Footprint Pvalue, and R-squared need to combine both directions, 
            ## so we remove it from the footprinting
            eachGrouping = group_by_cols[group_by_cols != 'Direction']

            ### Summary left and right heigh, turn Radius into Width, and take the min R-squared and max P-value
            HL = summaryStats[Direction == 'Upstream',][, 
                            .(Height_Left = mean(Height, na.rm = TRUE)), by = eachGrouping]
            HR = summaryStats[Direction == 'Downstream',][,
                            .(Height_Right = mean(Height, na.rm = TRUE)), by = eachGrouping]
            ## Take calculate Width by adding the left and right radius
            ## Find the overall footprint p-value using the max PValue from left and right
            ## Generate the overall footprint R-squared by taking the minimum R_squared for left and right model
            Width = summaryStats[, .(Width = sum(Radius),
                                    Footprint_PValue = max(PValue), 
                                     R_squared= min(R_squared)), 
                                     by = c('Sample',eachGrouping)]
            ## Now take the average of the width, the max of the footprint P-value and mean of the R_squared within each group
            Width_avg = Width[,.(Width = mean(Width),
                                    Footprint_PValue = max(Footprint_PValue, 'fdr'), 
                                     R_squared= mean(R_squared)), 
                                     by = eachGrouping]
            ## Now test the radius of left and right within each group as compared to random distribution
            ## between 0 and the searchWidth. This is a way of verifying that we're getting a true height, 
            ## and not a random value
            widthTest = summaryStats[,.(Radius_PValue = 
                                    wilcox.test(Radius, sample(0, searchWidth, length(Radius)), 
                                                alternative = 'greater')$p.value), 
                                        by = group_by_cols]
            ## Now combine p-values for width of left and right by taking the max p_value
            widthTest = widthTest[, .(Width_PValue = max(Radius_PValue)), 
                                     by = eachGrouping]
            #Now combine the results
            finalStats = HL[HR[
                            Width_avg[widthTest, on = eachGrouping], 
                               on = eachGrouping],
                                on = eachGrouping]
            
        
        
          }else{
            
            ## Not implemented yet
            #### Now run additional statistics if background is define
            groups = unique(unlist(summaryStats[,groupColumn]))
            
            if(!any(groups == background)){

                stop('Background group not found within groups column')
            }
            
              ## First check if the user wants a differential between motifs or between cell types
        group_by_cols = c('CellType', 'Motif', 'Direction')
        if(any(groupColumn %in% c('CellType', 'Motif'))){
           
            group_by_cols = group_by_cols[!group_by_cols %in% groupColumn ]
            
        }
        sumValues = c('Height', 'Radius', 'R_squared', 'PValue', 'Baseline', 'Depth')
        
        ### Pivot two heights to left and right, and generate a Width parameter for the sum of Radius
        
        ## Test whether the radius is significantly different from random. 
        widthTest = summaryStats[,.(Width_Pvalue = 
                                wilcox.test(Radius, sample(0, searchWidth, length(Radius)), 
                                            alternative = 'greater')$p.value), 
                                    by = group_by_cols]
        summaryStats = summaryStats[widthTest, on = c('CellType','Motif','Direction')]
        
        meanSumStats = summaryStats[,lapply(.SD, mean, na.rm = TRUE), by = c(group_by_cols), .SDcols = sumValues]
        
        evalStats = meanSumStats[normSumStats, on = group_by_cols]
        
        
        summaryStats2 = suppressMessages(do.call('rbind', 
            lapply(groups[groups != background], function(GG){

                Baseline = dplyr::summarise(dplyr::group_by_at(summaryStats, dplyr::vars(dplyr::one_of(group_by_cols))),
                         Metric = 'Baseline',
                         Foreground = GG,
                         Background = background,
                         MeanBackground = mean(Baseline[!!rlang::sym(groupColumn) == background]),
                         MeanDiff = mean(Baseline[!!rlang::sym(groupColumn) == GG]) - 
                                        mean(Baseline[!!rlang::sym(groupColumn) == background]),
                         P_Value = stats::wilcox.test(Baseline[!!rlang::sym(groupColumn) == GG], 
                                                Baseline[!!rlang::sym(groupColumn) == background], paired=FALSE)$p.value) 

                Width = dplyr::summarise(dplyr::group_by_at(summaryStats, dplyr::vars(dplyr::one_of(group_by_cols))),
                         Metric = 'Width',
                         Foreground = GG,
                         Background = background,
                         MeanBackground = mean(Width[!!rlang::sym(groupColumn) == background]),
                         MeanDiff = mean(Width[!!rlang::sym(groupColumn) == GG]) - 
                                        mean(Width[!!rlang::sym(groupColumn) == background]),
                         P_Value = stats::wilcox.test(Width[!!rlang::sym(groupColumn) == GG], 
                                                Width[!!rlang::sym(groupColumn) == background], paired=FALSE)$p.value) 

                 Depth = dplyr::summarise(dplyr::group_by_at(summaryStats, dplyr::vars(dplyr::one_of(group_by_cols))),
                         Metric = 'Depth',
                         Foreground = GG,
                         Background = background,
                         MeanBackground = mean(Depth[!!rlang::sym(groupColumn) == background]),
                         MeanDiff = mean(Depth[!!rlang::sym(groupColumn) == GG]) - 
                                        mean(Depth[!!rlang::sym(groupColumn) == background]),
                         P_Value = stats::wilcox.test(Depth[!!rlang::sym(groupColumn) == GG], 
                                                Depth[!!rlang::sym(groupColumn) == background], paired=FALSE)$p.value)

                do.call('rbind', list(Baseline, Width, Depth))
            })))
            
            
          }
        
    }else if(is.null(groupColumn)){
        
        # If there's no group column, just return the basic footprinting values for each 'sample' (could be a group-average, or an actual sample. 
        
        ### Summary left and right heigh, turn Radius into Width, and take the min R-squared and max P-value
        HL = summaryStats[Direction == 'Upstream',][, .(Height_Left = Height), by = .(Sample, CellType, Motif)]
        HR = summaryStats[Direction == 'Downstream',][, .(Height_Right = Height), by = .(Sample, CellType, Motif)]
        Width = summaryStats[, .(Width = sum(Radius), Footprint_PValue = max(PValue), R_squared= min(R_squared)), 
                                 by = .(Sample, CellType, Motif)]
        finalStats = HL[HR[Width, on = .(Sample,CellType,Motif)], on = .(Sample, CellType,Motif)]
        
    }
    
    return(finalStats)
}
           
#' @title Generate motif footprint metrics and/or statistics across samples, cell types, and motifs
#'
#' @description Identifies the width and depth of a given motif for a given group & cell type. Requires the output of motifFootprint. If a groupColumn and background is provided, the function will also conduct wilcoxon tests on the motif baseline, depth, and width between groups. 
#'
#' @param motifSE A SummarizedExperiment with motif footprinting information, from motifFootprint
#' @param footprint Optional string, to describe which footprint within the motifSE to analyze. Default is NULL, at which point it will pull from all footprint present. 
#' @param callFootprint A boolean, to determine whether to call footprints as significant across samples. Requires sample-specific motif footprinting.
#' @param motifSize A number describing the approximate motif size you want to use, when looking at the insertion rate over the motif core. Default is 6 (3 bp up and 3 bp downstream)
#' @param searchWidth A number, describing the maximum motif footprint width. The function will look for the maximum insertions within +/- this number of the motif center. Default is 50.  
#' @param groupColumn An optional string, that will contain the group-level labels for either calling footprints as significant, or for comparing motif footprintings stats betweem groups. Default is null, at which point no wilcoxon tests will be conducted. groupColumn can be a sample-level group from the MOCHA's object colData slot, 'CellType' (if you want to compare across cell types), or 'Motif' if you want to compare across motifs. 
#' @param background An optional string if you want to compare two different conditions. This string is the background group that you want to use for running wilcoxon tests. 
#'
#' @return a data.frame, containing motif footprint stats. 
#'
#' @examples
#' \dontrun{
#' summaryStats <- MOCHA::motifStats(
#'   motifSE,
#'   footprint = "CD4 Naive_ARID5A", groupColumn = 'COVID_status', background = 'Negative')
#' }
#'
#' @export
#' @keywords exporting

plotMotifs <- function(motifSE, 
                       footprint = NULL, 
                       groupColumn = NULL, 
                       returnDF = FALSE, 
                       returnPlotList = FALSE,
                          relHeights = c(0.3, 0.7)){

    if(is.null(footprint)){
    
        footprint = names(SummarizedExperiment::assays(motifSE))
        
    }
    if(!any(names(SummarizedExperiment::assays(motifSE)) %in% footprint)){
    
        stop('footprint name was not found within the motifSE object')
        
    }
    
    colData1 = SummarizedExperiment::colData(motifSE)
    


    ## Pull out footprints across regions before running stats and add the direction feature. 
    if(length(footprint) > 1){
        
        motif_spec <- data.table::rbindlist(lapply(footprint, function(XX){
                    tmpMat = colMeans(SummarizedExperiment::assays(motifSE)[[XX]])
                    Sample= gsub('__-[0-9].*|__[0-9].*', '' ,names(tmpMat))
                    Position = gsub('__', '', gsub(paste0(unique(Sample), collapse='|'), '', names(tmpMat)))
                    data.table::data.table(Insertions = tmpMat, Sample = Sample, Position = as.numeric(Position), Footprint = XX)
            }))
        ## Process single location data as well
        fullDT = data.table::rbindlist(lapply(footprint, function(XX){
                    tmpMat = SummarizedExperiment::assays(motifSE)[[XX]]
                    tmpMat2 = as.data.table(tmpMat)
                    tmpMat2$Location = rownames(tmpMat)
                    dt1 = data.table::melt(tmpMat2,
                                        id.vars = 'Location',
                                        measureVar = colnames(tmpMat),
                                           variable.name = 'Index',
                                           value.name = 'Insertions')
            }))
        
    }else{
        
        tmpMat = colMeans(SummarizedExperiment::assays(motifSE)[[footprint]])
        Sample= gsub('__-[0-9].*|__[0-9].*', '' ,names(tmpMat))
        Position = as.numeric(gsub('__', '', gsub(paste0(unique(Sample), collapse='|'), '', names(tmpMat))))
        motif_spec <- data.table::data.table(Insertions = tmpMat, Sample = Sample, Position = Position, Footprint = footprint)
        ## Process single location data as well        
        tmpMat = SummarizedExperiment::assays(motifSE)[[footprint]]
        tmpMat2 = as.data.table(tmpMat)
        tmpMat2$Location = rownames(tmpMat)
        dt1 = data.table::melt(tmpMat2,
                            id.vars = 'Location',
                            measureVar = colnames(tmpMat),
                               variable.name = 'Index',
                               value.name = 'Insertions')
        dt1$Footprint = footprint
    }    
    
    ## Split out Cell type and motif info
    motif_spec <- motif_spec[, ':='(CellType = gsub('__.*','', Footprint), 
                                    Motif = gsub('.*__','', Footprint))]
    dt1 <- dt1[, ':='(CellType = gsub('__.*','', Footprint), 
                                    Motif = gsub('.*__','', Footprint))]
    ## Now pull out Sample vs Position info
    dt1 <- dt1[, Sample := gsub('__-[0-9].*|__[0-9].*', '' , Index)]
    allSamples = paste0(unique(dt1$Sample), collapse='|')
    dt1 <- dt1[, Position := as.numeric(gsub('__', '', 
                            gsub(allSamples, '', Index)))]
    
     ### Add in sample metadata, if motif footprinting done on the sample level. 
    if(any(names(motifSE@metadata) == 'metadata')){
    
        metadf = as.data.table(motifSE@metadata[['metadata']])
        
        if(any(metadf$Sample %in% motif_spec$Sample)){
            motifDT <- motif_spec[metadf, on = .(Sample)]
            singleDT <- dt1[metadf, on = .(Sample)]
        }
    }
    
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
        singleDT = singleDT[, .(Insertions = mean(Insertions, na.rm = TRUE)),
                                     by =  c('Location', group_by_cols)]
        singleDT$PlotGroup = paste(groupColumn, singleDT[,get(groupColumn)], singleDT$Footprint, sep = '_')
        
    }else{
        motifDT$PlotGroup = paste(motifDT[,Sample], motifDT$Footprint, sep = '_')
        singleDT$PlotGroup = paste(groupColumn, singleDT[,Sample], singleDT$Footprint, sep = '_')
    }
    
    ##Order the singleDT by average strength across conditions
    avgStrength = singleDT[, .(totalInsertions = sum(Insertions, na.rm = TRUE)), by = .(Location, PlotGroup)]
    avgStrength = avgStrength[order(totalInsertions,decreasing=TRUE )]
    singleDT$PlotGroup = factor(singleDT$PlotGroup, levels = rev(unique(avgStrength$PlotGroup)))
    singleDT$Location = factor(singleDT$Location, levels = rev(unique(avgStrength$Location)))
    
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
            ggplot2::scale_x_continuous(expand=c(0,0)) + 
            ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                 axis.ticks.x = ggplot2::element_blank())
    p2 = ggplot2::ggplot(singleDT, ggplot2::aes(x = Position, fill = Insertions, y = Location)) + 
            ggplot2::geom_tile() +
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
                          plot.background = ggplot2::element_rect(fill = 'white', colour = 'white'))+
            ggplot2::guides(r = ggplot2::guide_axis(angle = 45)) +
            ggplot2::theme(plot.margin = ggplot2::unit(c(0, 0.1, 0.1, 0.1), "cm"))
    
    returnList = list('MotifAverage' =  p1 , 'LocationSpecific' = p2)
    if(returnPlotList){
        
        return(returnList)
        
    }
    
    cp1 = cowplot::plot_grid(plotlist = returnList,
                                 ncol=1, align = 'v',
                                 axis = "brl",
                                rel_heights = relHeights)
         
    return(cp1)  
    
}
                                
callFootprint <- function(summaryStats, groupColumn = NULL){
    

    group_by_cols = c('CellType', 'Motif', groupColumn)

    Depth = dplyr::summarise(dplyr::group_by_at(summaryStats, dplyr::vars(dplyr::one_of(group_by_cols))),
                         Metric = 'Depth',
                         Foreground = 'Height',
                         Background = 'Baseline',
                         MeanBackground = mean(Baseline),
                         MeanDiff = mean(Depth),
                         P_Value = stats::wilcox.test(Height, Baseline, paired=TRUE, alternative = 'greater')$p.value)
    
    
}

