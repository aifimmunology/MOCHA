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
#' @references Pimentel, Ronald Silva, "Kendall's Tau and Spearman's Rho 
#'               for Zero-Inflated Data" (2009). Dissertations. 
#'
#' @export


plot_coaccessibility <- function(zi_mat, index, peak_mat,fname, addLabel=F){
    
    df = data.frame(
        Values_X = as.numeric(peak_mat[zi_mat$Peak1[index], 2:ncol(peak_mat)]),
        Values_Y = as.numeric(peak_mat[zi_mat$Peak2[index], 2:ncol(peak_mat)])
        )

    
        fname =paste(fname,'.png',sep='')

    ## Calculate 3-metrics of correlation
    spear_cor = cor(df[,1], df[,2], method='spearman')
    zi_spear_cor = weightedZISpearman(df[,1], df[,2])

    ## Filter out nonzeros 
    df=as.data.table(df)
    df$Values_Y2 <- df$Values_Y
    df$Values_Y2[df$Values_Y2==0] <- NA
    ## calculate spearman on nonzeros 
    nonzero_df = df[Values_X != 0 & Values_Y != 0]
    spear_nonzero = cor(nonzero_df$Values_X, nonzero_df$Values_Y, method='spearman')
    
    nonzero_fit <- lm(log2(Values_Y+1) ~ log2(Values_X+1), nonzero_df)
    int = nonzero_fit$coefficients[1]
    slp = nonzero_fit$coefficients[2]
    
    res=data.frame(Spearman=spear_cor,
               ZI_Spearman=zi_spear_cor,
               Spearman_NonZero=spear_nonzero)
    
    x1 = log2(min(nonzero_df$Values_X))-1
    x2 = log2(max(nonzero_df$Values_X))+1
    y1 = int + slp* x1
    y2 = int + slp* x2
    
    png(fname)
    p <- ggplot(df,
                aes(x=log2(Values_X+1),
                    y=log2(Values_Y+1)
                   ))+    geom_point()+
        ggtitle(paste("Co-Accessibility"))+xlim(0,17)+ylim(0,17)+
            xlab(zi_mat$Peak1[index])+ylab(zi_mat$Peak2[index])+
      theme(legend.position = 'none', 
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        strip.text.x = element_text(size=18),
        axis.text.x = element_text(size=14, angle = 90),
        axis.text.y = element_text(size=14),
          )
    
    if(addLabel){
            p= p+geom_text(  x = c(3), y = c(8),
                  label=paste('ZI-S=',
                               round(zi_spear_cor,2)),
                                       size=10)+
            geom_text(  x = c(3), y = c(6),
                  label=paste('S=',
                               round(spear_cor,2)),
                                       size=10,
                             col='blue')+
            geom_smooth(method='lm', se=FALSE,size=2)+
            geom_segment(aes(x = x1, xend =x2 , 
                                 y = y1,  yend = y2),size=2
                             )
        }

    
    print(p)
    dev.off()
    
}
