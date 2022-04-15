#' @title \code{determine_dynamic_range}
#'
#' @description \code{make_prediction} is an R helper function, part of the single-cell peak calling
#' algorithm scMACS by (XXX et al, 2021) that determines which genomic regions, or bins,
#' will be used for de-novo peak calling. The function "make_prediction" applies the
#' logistic regression model predictions and then calls peaks that exceed a given
#' threshold.
#'
#'
#' @param X an intensityMatrix output from \code{calculate_intensities}
#' @param finalModel is a matrix with coefficients and an index indicating
#'        the number of cell used to train that model
#' @param thresholdModel is an internal threhsod
#'
#' @return the original intensityMatrix with the two intensity parameters required
#' to calculate the probability of a (+) peak, with an additional two columns
#' that include the prediction probability
#'
#' @details The technical details of the algorithm are found in XX.
#'
#' @references XX
#'
#' 

make_prediction <- function(X, finalModelObject, thresholdModel){

    ### Model was trained on varying
    ### cell abundances. Identify 
    ### Number of cells in the sample
    ### to apply appropriate model
    cell_model = X$numCells[1]
    
    ### If the number of cells < 5
    ### We do not make predictions 
    
    if(cell_model  < 5){
        print('Cannot make peak calls with < 5 cells, returning NULL object')
        return(NULL)
    } else {
    
       ### Apply model fit based on the 
       ### cell abundance. First Part is a loess
       ### Fit, final part is a linear fit,
       ### Middle part is an average between 
       ### Loess & Linear fits. 
       if(cell_model<=100000){

           ## Loess Fit 
           Intercept <- predict(finalModelObject$Loess$Intercept, data.frame(NumCells=cell_model))
           Total <- predict(finalModelObject$Loess$Total, data.frame(NumCells=cell_model))
           Max <- predict(finalModelObject$Loess$Max, data.frame(NumCells=cell_model))
           tmpModel = c(Intercept, Total, Max)
           tmpModel = as.matrix(tmpModel)


       } else if(cell_model>100000 & cell_model <= 140000){ 

           ## Average of Loess & Linear Fit           
           Intercept <- mean(predict(finalModelObject$Loess$Intercept, 
                                     data.frame(NumCells=cell_model)),
                           predict(finalModelObject$Linear$Intercept, 
                                   data.frame(NumCells=cell_model))
                           )
           
           Total <- mean(predict(finalModelObject$Loess$Total, data.frame(NumCells=cell_model)),
                           predict(finalModelObject$Linear$Total, data.frame(NumCells=cell_model))
                           )
           Max <- mean(predict(finalModelObject$Loess$Max, data.frame(NumCells=cell_model)),
                           predict(finalModelObject$Linear$Max, data.frame(NumCells=cell_model))
                           )       

           tmpModel = c(Intercept, Total, Max)
           tmpModel = as.matrix(tmpModel)


       } else {
           ## Linear Fit           
           Intercept <- predict(finalModelObject$Linear$Intercept, data.frame(NumCells=cell_model))
           Total <- predict(finalModelObject$Linear$Total, data.frame(NumCells=cell_model))
           Max <- predict(finalModelObject$Linear$Max, data.frame(NumCells=cell_model))
           tmpModel = c(Intercept, Total, Max)
           tmpModel = as.matrix(tmpModel)
        }
           

   
        }
    
    ### Create Design Matrix 
    designX <- X[,c('TotalIntensity','maxIntensity')]
    designX$Intercept =1 
    designX <- designX[,c('Intercept','TotalIntensity','maxIntensity')]   
    colnames(designX) <- c('Intercept','Total','Max')

    ### Apply Prediction Model to 
    ### Obtain Probability of Accessibility
    z = as.matrix(designX) %*% tmpModel
    preds = 1/(1+exp(-z))

    ### Round Predictions & 
    ### Create a Prediction Strength Feature 
    X$Prediction = round(preds ,4)
    X$PredictionStrength = X$TotalIntensity

    ### Smoothened Thresholding 
    ### To create stable model 
    ### Across varying cell abundances
    ### Using the Youden Index 
    ### to create the smoothened threshold
    ### model. Loess & Linear fits were 
    ### used to create a smoothened 
    ### threshold to be robust against
    ### sparsity & noise. A global offset
    ### is included to account for the fact
    ### that the nFrags normalizing factor 
    ### will shift probabilities of accessiblities
    ### slightly. 
    
    if(cell_model <= 2000){ 
           newdata = data.frame(Ncells = cell_model)
           adaptiveThreshold = predict(thresholdModel$loess_low, 
                                       newdata=newdata)+.02
    
    } else if(cell_model > 2000 & cell_model < 5000){

           newdata = data.frame(Ncells = cell_model)
           adaptiveThreshold = mean(predict(thresholdModel$loess_low, 
                                       newdata=newdata),
                                    predict(thresholdModel$loess_mid, 
                                       newdata=newdata)   
                                    )+.03
    } else if(cell_model >= 5000 & cell_model <= 60000){
        
           newdata = data.frame(Ncells = cell_model)
           adaptiveThreshold = predict(thresholdModel$loess_mid, 
                                       newdata=newdata)+.02
        
    } else if(cell_model > 60000 & cell_model <= 150000){
           newdata = data.frame(Ncells = cell_model)
           adaptiveThreshold = predict(thresholdModel$linear, 
                                       newdata=newdata)+.02                
    } else if(cell_model > 150000){
           ### Fix threshold based on
           ### Youden Index Convergence 
           adaptiveThreshold = 0.20
    }
           
   ## enforce threshold to be numeric object
   adaptiveThreshold=as.numeric(adaptiveThreshold)
   
   ### a global tolerance offset 
   ### included to account for numerical rounding
   ### and call peaks based on probabilities > threshold  
   tolerance=0.001
   X$peak = X$Prediction > round(adaptiveThreshold,4) + tolerance

    return(X)


        
    
}
