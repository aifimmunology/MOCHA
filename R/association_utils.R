

#' @title \code{visualizeAssociations}
#'
#' @description \code{visualizeAssociations} Takes the output of multiModalModeling (or other) and visualizes the top interaction features
#'  for one of the modalities. 
#' @param interAct A dataframe of coefficients related to a modality, extracted from multiModalModeling or other. This should be prefiltered for significance (FDR < 0.1)
#' @param dataType1 Boolean flag. This determines whether it should look at the Modality1 or Modality2 column (i.e. the first or second assay given to multiModalIntegration)
#' @param dataTypeName A string to label the modality in the plots generated. 
#' @param threshold A number, describing how many points to label. By default, it will label the top 20 points. 
#' @param max.overlaps The maximum number of overlapping labels, a parameter which is passed to ggrepel when labeling points on the plot. 
#' @param verbose Boolean flag to determine verbosity. 
#' 
#' @return A ggplot object, showing number of interactions on the y-axis and the rank order of measurements on the x-axis 
#'    (Measurements rank-ordered by the number of interactions they have with the other modality.)
#'
#' @details
#'
#'
#' @export

visualizeAssociations <- function(interAct, dataType1 = TRUE, dataTypeName = 'TF', threshold = 20,  max.overlaps = 10,
                                 verbose = FALSE){
    if(dataType1){
        factor1 = 'Modality1'
    }else{
        factor1 = 'Modality2'
    }
    df1 <- as.data.frame(table(interAct[,factor1]))
    colnames(df1) = c(dataTypeName,'Interactions')
    df1 <- df1[order(df1$Interactions,decreasing =T),]
    df1$Rank <- as.numeric(c(1:length(df1$Interactions)))
    df1$Factor <- df1[,dataTypeName]
    df1$Factor[df1$Rank > threshold] = NA
    if(verbose){
       head(df1)
    }
    ggplot(df1, aes(x = Rank, y = Interactions, label = Factor)) + ggrastr::geom_point_rast() + theme_minimal() +
    ggtitle(factor1) + 
    xlab('Rank') + ylab('Number of Interactions')  + ggrepel::geom_text_repel(max.overlaps = max.overlaps)
}


#' @title \code{volcanoPairs}
#'
#' @description \code{volcanoPairs} Takes the output of multiModalModeling (or other) and generates a volcano plot of Estimate vs -log10(FDR) for that pair
#' @param interAct A dataframe of coefficients related to a modality, extracted from multiModalModeling or other. This should be not be filtered. 
#' @param factor1 Name for the 
#' @param dataTypeName A string to label the modality in the plots generated. 
#' @param FDR_threshold A number, which thresholds which points are considered significant. It will also add a red line to the plot to mark significance.
#' @param topLimit  A number, limiting how many points to label. By default, it will label the top 20 points by FDR, if there are more than 20 points.
#' @param max.overlaps The maximum number of overlapping labels, a parameter which is passed to ggrepel when labeling points on the plot. 
#' @param addSignificanceLine Boolean flag to determine whether or not to add a line at the FDR significance threshold. 
#' 
#' @return A ggplot object, showing number of interactions on the y-axis and the rank order of measurements on the x-axis 
#'    (Measurements rank-ordered by the number of interactions they have with the other modality.)
#'
#' @details
#'
#'
#' @export

volcanoPairs <- function(model, variable, FDR_threshold = 0.05,  topLimit = 20, max.overlaps = 10, returnDF = FALSE, addSignificanceLine = TRUE){


    interAct <- as.data.frame(SummarizedExperiment::assays(model)[[variable]])
    #Generate range for plotting
    xRange <- range(interAct$Estimate)*0.10
    yRange <- range(-log10(interAct$FDR))*0.10
    #Generate labels for each plot
    interAct$Pair = rownames(interAct)
    #Filter out which pairs to label. First by FDR threshold, then by the topLimit parameter. 
    interAct$Pair[interAct$FDR >= FDR_threshold] = NA
    if(sum(interAct$FDR <= FDR_threshold)){
        topPairs <-  interAct$Pair[order(interAct$FDR)[c(1:topLimit)]]
        interAct$Pair[!interAct$Pair %in% topPairs] = NA
    }

    if(returnDF){
      return(interAct)
    }

    #Generate ggplot.  
    p1 <- ggplot(interAct, aes(x = Estimate, y = -log10(FDR), label = Pair)) + 
           geom_point() +   coord_cartesian(clip = 'off') + 
                theme(plot.margin = margin(1,1,1.5,1.2, "cm")) +
            ggrepel::geom_text_repel(max.overlaps = max.overlaps, force = 4,
                            max.iter = 100000) +
    theme_minimal() 
    if(addSignificanceLine){
      p1 <- p1+ 
            geom_hline(yintercept = -log10(FDR_threshold), color = 'red', alpha = 0.25) 
    }
    return(p1)

}