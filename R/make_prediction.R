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

make_prediction <- function(X, finalModel){

    cell_model = X$numCells[1]
    
    if(cell_model  < 20){
        stop('Cannot make peak calls with < 20 cells')
    }
    
    designX = cbind(rep(1, nrow(X)),
              X[,c('lambda1','lambda2')]
    )

    model.idx <- which(cell_model == finalModel$numCells)
    
    if(length(model.idx)==0){
        distances <- (abs(cell_model - finalModel$numCells))
        sorted.distances <- sort(distances, decreasing=F)      
        
        model.idx <- which(distances %in% sorted.distances[1])
        
        if(length(model.idx)==2){
            
           interpolatedModel <- (finalModel[model.idx[1],] + finalModel[model.idx[2],])/2

        } else if(finalModel$numCells[model.idx] > cell_model){
            
           interpolatedModel <- (finalModel[model.idx,] + finalModel[model.idx-1,])/2

        } else {
            interpolatedModel <- (finalModel[model.idx,] + finalModel[model.idx+1,])/2
            
        }    
        
        
        tmpModel <- interpolatedModel[,c('Intercept','lambda1','lambda2')]
        tmpModel = as.matrix(tmpModel)
        tmpModel = t(tmpModel)
        designX = as.matrix(designX)
        preds = designX%*% tmpModel
        preds = 1/(1+exp(-preds))

        X$Prediction = preds 
        X =makeGRangesFromDataFrame(X, keep.extra.columns=T)
        X$PredictionStrength = X$lambda1
        X$Peak = X$Prediction > 0.5
        
          return(X)
        
    } else {
    
        tmpModel <- finalModel[model.idx,c('Intercept','lambda1','lambda2')]
        tmpModel = as.matrix(tmpModel)
        tmpModel = t(tmpModel)
        designX = as.matrix(designX)
        preds = designX%*% tmpModel
        preds = 1/(1+exp(-preds))

        X$Prediction = preds 
        X =makeGRangesFromDataFrame(X, keep.extra.columns=T)
        X$PredictionStrength = X$lambda1
        X$Peak = X$Prediction > 0.5
        
        return(X)
        }

  }
