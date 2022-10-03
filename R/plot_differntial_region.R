#' @title \code{plot_differential_region}
#'
#' @description \code{plot_differential_region} allows you to graph a particular
#'              region across two different conditions
#'
#' @param tmp_mat: sample-peak matrix 
#' @param tileID: chromosome region to analyze 
#' @param group: vector of group indices (0,1)
#'
#' 
#' The output is a ggplot2 object showing the histograms
#' of the intensities at a region, by each condition.
#'
#' @references XX
#'
#' @export   

plot_differential_region<- function(tmp_mat,tileID,group){

    ## Filter out 
    ## 'tileID' column
    idx = which(tmp_mat$tileID == tileID)
    test_vec = as.numeric(tmp_mat[idx,1:ncol(tmp_mat)-1, with=F])

    ## log-transform the values 
    data=test_vec
   
  
    df = data.frame(
        Vals = data,
        Group= group,
        Group_label=ifelse(group==1, 'case','Control')
    )
    

        
        fname =paste(tileID,'.png',sep='')
    
        p <- ggplot(df,
                   aes(x=Vals, col=Group_label, fill=Group_label))+#geom_density()+
                   geom_histogram(aes(y=0.2*..density..),
                 alpha=0.5,position='identity',binwidth=0.2)+
    facet_wrap(~df$Group_label, ncol=1, scales='free')+
                ggtitle(tileID)+
        xlab('Log2 Intensity')+ylab('Frequency')+
            xlim(c(min(df$Vals)-1, max(df$Vals)+1))+
    theme(legend.position='none',
         axis.text.x=element_text(size=18),
         axis.title.x=element_text(size=18),          
         axis.title.y=element_text(size=18),         
         plot.title=element_text(size=18),                   
         axis.text.y=element_text(size=18)
          )
        return(p)
        
}