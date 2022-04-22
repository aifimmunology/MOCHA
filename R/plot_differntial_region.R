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
    test_vec = as.numeric(tmp_mat[idx,2:ncol(tmp_mat)])

    ## log-transform the values 
    data=log2(test_vec+1)
   
    ## conduct two part test
    two_part_results <- TwoPart(test_vec, group=group, test='wilcoxon', point.mass=0)

    ## permutation test
    n_a = sum(group)
    n_b = length(group) - n_a

    sampled_test_vec <- sample(test_vec, size=n_a+n_b, replace=F)
    permuted_pvalue <- TwoPart(sampled_test_vec, group=group, test='wilcoxon', 
                               point.mass=0)$pvalue

    
    ## filter non-zero values for group1
    nonzero_dx <- data[group==1]
    nonzero_dx=nonzero_dx[nonzero_dx!=0]
    
    ## filter non-zero values for group0
    nonzero_control <- data[group==0]
    nonzero_control=nonzero_control[nonzero_control!=0]
    
    df = data.frame(
        Vals = data,
        Group= group,
        Group_label=ifelse(group==1, 'case','Control')
    )
    

        
        fname =paste(res$Peak,'.png',sep='')
        png(fname)
        p <- ggplot(df,
                   aes(x=Vals, col=Group_label, fill=Group_label))+
               geom_histogram(aes(y=0.1*..density..),
                 alpha=0.5,position='identity',binwidth=0.1)+facet_wrap(~df$Group_label, ncol=1, scales='free')+ThemeMain+
                ggtitle(paste("Differential Accessibility",res$Peak,sep='\n'))+
        xlab('Log2 Intensity')+ylab('Frequency')+
            xlim(c(min(df$Vals)-1, max(df$Vals)+1))
        print(p)
        dev.off()
        print(paste('file saved to ',fname))
        
   }