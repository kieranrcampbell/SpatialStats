
##############################################
## Data processing utilities for SpatialPRo ##
## kieranrcampbell@gmail.com                ##
##############################################


###############################
## Variance Inflation Factor ##
###############################

#' Most Variance inflation factor (VIF) calculations involve supplying
#' a linear model, when in fact only the matrix of predictors is required.
#' VIF returns the variance inflation factors for predictors using only the
#' matrix of predictors
#'
#' @param X Predictor matrix
#'
#' @export
VIF <- function(X) {
    ## returns the vif coefficients without need for the
    ## entire model
    k <- dim(X)[2]
    sapply(1:k, function(i) {
        y <- X[,i] ; x <- X[,-i]
        fit <- lm(y ~ x)
        r <- summary(fit)$r.squared
        1/(1-r)
    })
}

#' Takes a matrix of predictors and sequentially removes them until the largest
#' VIF is below a threshold (default=10)
#'
#' @param X Matrix of predictors
#' @param vif.threshold Maximum VIF threshold allowed
#'
#' @export
VIFRemove <- function(X, vif.threshold=10) {
    x <- X
    while(TRUE) {
        v <- VIF(x)
        if(max(v) > vif.threshold) {
            rem <- which.max(v)
            x <- x[,-rem]
        } else {
            return(x)
        }
    }
}
