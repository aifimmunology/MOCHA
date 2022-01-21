#' @title \code{determine_dynamic_range}
#'
#' @description \code{make_prediction} is an R helper function, part of the single-cell peak calling
#' algorithm scMACS by (XXX et al, 2021) that determines which genomic regions, or bins,
#' will be used for de-novo peak calling. The function "make_prediction" applies the
#' logistic regression model predictions
#'
#'
#' @param X an intensityMatrix output from \code{calculate_intensities}
#' @param finalModel is a matrix with coefficients and an index indicating
#'        the number of cell used to train that model
#'
#' @return the original intensityMatrix with the two intensity parameters required
#' to calculate the probability of a (+) peak, with an additional two columns
#' that include the prediction probability
#'
#' @details The technical details of the algorithm are found in XX.
#'
#' @references XX
#'
#' @export

make_prediction <- function(X, finalModelObject, tolerance=0.001, thresholdCutoff=6){

    cell_model = X$numCells[1]
    
    if(cell_model  < 5){
        print('Cannot make peak calls with < 5 cells, returning NULL object')
        return(NULL)
    } else {
    

       if(cell_model<=80000){

           lambda0 <- predict(finalModelObject$LogFits[[1]], data.frame(numCells=cell_model))
           lambda1 <- predict(finalModelObject$LogFits[[2]], data.frame(numCells=cell_model))
           lambda2 <- predict(finalModelObject$LogFits[[3]], data.frame(numCells=cell_model))
           tmpModel = c(lambda0, lambda1, lambda2    )
           tmpModel = as.matrix(tmpModel)


       } else if(cell_model>80000 & cell_model <= 120000){ 


           lambda0 <- mean(predict(finalModelObject$LogFits[[1]], data.frame(numCells=cell_model)),
                           predict(finalModelObject$LinearFits[[1]], data.frame(numCells=cell_model))
                           )
           
           lambda1 <- mean(predict(finalModelObject$LogFits[[2]], data.frame(numCells=cell_model)),
                           finalModelObject$LinearFits[[2]]$coefficients[1]
                           )
           lambda2 <- mean(predict(finalModelObject$LogFits[[3]], data.frame(numCells=cell_model)),
                           predict(finalModelObject$LinearFits[[3]], data.frame(numCells=cell_model))
                           )       

           tmpModel = c(lambda0, lambda1, lambda2    )
           tmpModel = as.matrix(tmpModel)



       } else {

           lambda0 <- predict(finalModelObject$LinearFits[[1]], data.frame(numCells=cell_model))
           lambda1 <- finalModelObject$LinearFits[[2]]$coefficients[1]
           lambda2 <- predict(finalModelObject$LinearFits[[3]], data.frame(numCells=cell_model))

           tmpModel = c(lambda0, lambda1, lambda2    )
           tmpModel = as.matrix(tmpModel)    

        }
           

   
        }
    
    designX <- X[,c('lambda1','lambda2')]
    designX$Intercept =1 
    designX <- designX[,c('Intercept','lambda1','lambda2')]       

    z = as.matrix(designX) %*% tmpModel
    preds = 1/(1+exp(-z))

    X$Prediction = preds 
    X$PredictionStrength = X$lambda1

    if(cell_model <= 1000){
           ## controls against noise
           noiseZ <- c(1,min(X$lambda1), min(X$lambda2)) %*% tmpModel
           adaptiveThreshold <- max(0.5, 1/(1+exp(-noiseZ)))
           
    } else if(cell_model > 1000 & cell_model < 135000){
        
           newdata = data.frame(NCells = cell_model)
           adaptiveThreshold = predict(scMACS::thresh_model_youden, newdata=newdata)     
            
    } else if(cell_model > 135000){
            adaptiveThreshold = 0.029 # max value of loess model 
        
    }
           

   ### tolerance
   ### included to account for 
   ### precision errors in low 
   ### values 
   X$Peak = X$Prediction > adaptiveThreshold + tolerance

    return(X)


        
    
}
