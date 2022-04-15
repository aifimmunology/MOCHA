#' @title \code{plot_coaccessibility}
#'
#' @description \code{plot_coaccessibility} allows you to visually determine 
#'              whether two regions are correlated. 
#'
#'
#' @param mat: sample-peak matrix with regions to analyze
#' @param numCores: integer to determine # of paralell cores
#'
#' @return a 3-column data.frame containing 
#'         - zi_mat= Melted Zero-inflated Spearman Correlation Matrix
#'         - index = numeric integer indicating which combination to plot
#'         - peak_mat= peak-sample matrix
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
#' @references XX
#'
#' @export


plot_coaccessibility <- function(zi_mat, index, peak_mat){
    
    df = data.frame(
        Values_X = as.numeric(peak_mat[zi_mat$Peak1[index], 2:ncol(peak_mat)]),
        Values_Y = as.numeric(peak_mat[zi_mat$Peak2[index], 2:ncol(peak_mat)])
        )

    
    
    ## Calculate 3-metrics of correlation
    spear_cor = cor(df[,1], df[,2], method='spearman')
    zi_spear_cor = weightedZISpearman(df[,1], df[,2])

    ## Filter out nonzeros 
    df=as.data.table(df)
    
    ## calculate spearman on nonzeros 
    nonzero_df = df[Values_X != 0 & Values_Y != 0]
    spear_nonzero = cor(nonzero_df$Values_X, nonzero_df$Values_Y, method='spearman')
    
    res=data.frame(Spearman=spear_cor,
               ZI_Spearman=zi_spear_cor,
               Spearman_NonZero=spear_nonzero)

    fname =paste(regions$region[index],'.png',sep='')
    png(fname)
    p <- ggplot(df,
                aes(x=log2(Values_X+1),
                    y=log2(Values_Y+1)
                   ))+    ThemeMain+geom_point()+
        ggtitle(paste("Co-Accessibility"))+xlim(0,17)+ylim(0,17)+
            xlab(regions$Peak1[index])+ylab(regions$Peak2[index])+
    geom_text(  x = c(10), y = c(2),
              label=paste('ZI-Spearman=',
                           round(zi_spear_cor,2)),
                                   size=8)+
        geom_text(  x = c(10), y = c(4),
              label=paste('Spearman=',
                           round(spear_cor,2)),
                                   size=8)+
            geom_text(  x = c(10), y = c(6),
              label=paste('NonZero Spearman=',
                           round(spear_nonzero,2)),
                                   size=8)
    
    print(p)
    dev.off()
    
}
