#' finalModelObject
#'
#' Trained MOCHA models - LOESS and linear regression
#'
#' @format A list of lists containing 2 items: "Loess" and "Linear" each with "Total" "Max" and "Intercept"  
#' \describe{
#' \item{Loess}{LOESS model}
#' \item{Linear}{Linear model}
#' }
"finalModelObject"

#' youden_threshold
#'
#' Trained regression model for predicting a cutoff threshold 
#' for peak calling.
#' Call:
#' loess(formula = OptimalCutpoint ~ Ncells, data = thresh_df)
#' 
#' Number of Observations: 27 
#' Equivalent Number of Parameters: 5.98 
#' Residual Standard Error: 0.02121 
#'
#' @format A list of 18 regression variables
#' 
"youden_threshold"