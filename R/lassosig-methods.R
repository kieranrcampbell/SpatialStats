
#####################################################################
## High dimensional lasso p-vals                                   ##
## Meinshausen, Nicolai, Lukas Meier, and Peter Buhlmann.          ##
## "P-values for high-dimensional regression."                     ##
## Journal of the American Statistical Association 104.488 (2009). ##
##                                                                 ##
## kieran.campbell@dpag.ox.ac.uk                                    ##
#####################################################################

#' High dimensional significance test for LASSO using the multi-sample split method
#' 
#' @description
#' Provides p-values for lasso regression. This method implements the multi-sample splitting method for significance testing in a high-dimensional
#' regression context. The basic idea is to split a sample in two, perform variable selection using LASSO
#' on one half and derive p-values using ordinary least squares (OLS) on the other half. 
#' 
#' Note that the penalisation parameter \eqn{\lambda} and \emph{s} are identical - they are named
#' this way for consistency with the package \link{glmnet}.
#' 
#' This method is only implemented for a single response variable, since general lasso regression
#' requires the same set of parameters to be selected for every response variable, which is
#' overly restrictive in some cases.
#' 
#' @details
#' The method works by partitioning the dataset randomly in two halves. Lasso regression is
#' performed on one half, and using a particular value of the penalisation parameter lambda
#' then a subset of the predictor variables are chosen. Ordinary least squares regression is
#' then performed on the other half of the data. If \eqn{S} variables are chosen for a given
#' split, then the p-values are Bonferroni moderated to \eqn{S.p}. The p-values of all variables
#' not selected for a given split are then set to 1. This process is repeated \code{B} times, and
#' subsequently \code{B} sets of p-values are generated. These p-values are then aggregated across
#' splits to provide a given p-value for each predictor variable. For full details see the 
#' original paper. The aggregation requires an extra parameter \eqn{\gamma_min}, which is recommended
#' to be 0.05 (and set by the parameter \code{gamma.min}). 
#' 
#' Care must be taken with regards to the number of measurements (\code{length(y)}) and 
#' the number of folds for cross-validation (\code{nfolds}). The package \code{glmnet}
#' requres at least 3 samples in a cross-validation split in finding the optimal
#' \eqn{\lambda} (\emph{not the same as a multi-sample split}). Therefore, if we start
#' with \eqn{N} samples, \code{glmnet} receives at least \eqn{floor(N/2)} which it then splits
#' in to \code{nfolds} for cross validation. As such we necessarily need
#' \deqn{
#' floor((floor(N/2))/nfolds) > 3
#' }
#' which is safely satisfied provided \eqn{N > 6*nfolds + 3}
#'
#' When choosing \code{B}, there is a trade-off between bias and efficiency. A larger \code{B}
#' will lead to a less biased result (i.e. less sensitive to the random sampling of folds) but
#' can require significantly more computation time. A heuristically 'good' value is \code{B=50}.
#' 
#' The choice of \eqn{\lambda} = \code{s} is detailed in the \code{glmnet} package. For standard
#' problems the best choice may be \code{lambda.min}, though if you are specifically trying to
#' minimise the number of parameters necessary, \code{lambda.1se} (a one, not an L) is a good choice.
#' Alternatively it may be advantageous to select a fixed number of parameters on every split. This
#' can be performed by setting \code{s="usefixed"} and \code{fixedP} to the desired number of 
#' parameters.
#' 
#' Occasionally it is necessary to force the inclusion of predictors into the OLS significance testing.
#' These can be included by setting \code{include} to the numeric indices  (i.e. the column numbers) 
#' of the predictors to force-include. 
#' 
#' Note that force exclusion of an intercept in OLS (by setting \code{intercept = FALSE}) can seriously
#' bias results - only do this if you are sure at \eqn{x_i = 0} for all \eqn{i} that \eqn{y = 0} and that
#' all relationships are perfectly linear.
#' 
#' @references  Meinshausen, Nicolai, Lukas Meier, and Peter Buhlmann.
#' "P-values for high-dimensional regression."
#' Journal of the American Statistical Association 104.488 (2009).
#' 
#' @param y Response vector
#' @param X Design matrix
#' @param B Number of times to partition the sample
#' @param s The value of lambda to use in lasso. Can be:
#' \itemize{
#' \item{An actual value of lambda you have determined}
#' \item{\code{lambda.min} The lambda that minimises the mean squared error from n-fold cross-validation (recommended)}
#' \item{\code{lambda.1se} The lambda that accounts for the smallest number of parameters and is within 1 standard error from \code{lambda.min}}
#' \item{\code{halfway} The lambda that is halfway between \code{lambda.min} and \code{lambda.1se}}
#' \item{\code{usefixed} Use the smallest lambda that accounts for a given number of parameters (set by \code{fixedP})}
#' }
#' @param include Set of predictors to be force-included in OLS analysis
#' @param gamma.min Lower bound for gamma in the adaptive search for the best p-value (default 0.05)
#' @param fixedP The fixed number of parameters to use (if \code{s == "usefixed"})
#' @param nfolds Number of folds of cross-validation in the glmnet n-fold crossvalidation
#' @param intercept Whether to include an intercept in the OLS regression (default = \code{TRUE})
#' 
#' @return A vector of p-values, where the \eqn{i}th entry corresponds to the p-value for the predictor
#' defined by the \eqn{i}th column of \code{X}.
#' 
#' 
#' @export
#' @author Kieran Campbell \email{kieran.campbell@@dpag.ox.ac.uk}
#' @examples 
#' \dontrun{
#' set.seed(123)
#' X <- matrix(rnorm(300),ncol=3)
#' colnames(X) <- paste("predictor",1:3,sep="")
#' y <- rnorm(100, X[,1] - 2*X[,3],1)
#' pvals <- LassoSig(y, X, B=50, s="lambda.min")
#' 
#' ## Output:
#' ## predictor1   predictor2   predictor3 
#' ## 4.567217e-07 1.000000e+00 1.187136e-17 
#' }
#' 
LassoSig <- function(y,X, B=100, s=c("lambda.min","lambda.1se","halfway","usefixed"),
                     gamma.min=0.05, fixedP=NULL, include=NULL, nfolds=10, intercept = TRUE) {
  
  if(nrow(X) != length(y)) stop("Number of samples doesn't match")
  if(length(y) <= 6 * nfolds + 3) stop("Number of samples too small for required cross validation folds. See help('LassoSig') for more details")
  
  p.mat <- replicate(B, doSingleSplit(y, X, s, fixedP, include, nfolds, intercept) )
  ##mode(p.mat) <- "numeric" # weird
  
  adjusted.pvals <- apply(p.mat, 1, adaptiveP, gamma.min)
  
  names(adjusted.pvals) <- colnames(X)
  return(adjusted.pvals)
}

## Performs a single split of the data into floor(N/2) and N - floor(N/2)
## groups and performs lasso & LS estimation
doSingleSplit <- function(y, X, s, fixedP, include, nfolds, intercept) {
  n.samp <- length(y)
  sample.in <- sample(1:n.samp, size = floor(n.samp / 2))
  sample.out <- setdiff(1:n.samp, sample.in)
  
  y.in <- y[sample.in] ; y.out <- y[sample.out]
  X.in <- X[sample.in,] ; X.out <- X[sample.out,]
  
  predictors <- doLasso(y.in, X.in, s, fixedP, nfolds)
  
  p.vals <- doLSReg(y.out, X.out, predictors, include, intercept)
  return(p.vals)
}

## Performs lasso on y & X, with custom option of s to be midway point between lambdamin
## and lambda within 1se
doLasso <- function(y,X,s, fixedP = NULL, nfolds) {
  
  cv.fit <- cv.glmnet(X,y, standardize=FALSE)
  if(s == "halfway") {
    s <- (cv.fit$lambda.min + cv.fit$lambda.1se) / 2
  }
  
  if(s == "usefixed") {
    if(is.null(fixedP)) stop("Minimum predictors selected but fixedP is NULL")
    fit <- cv.fit$glmnet.fit
    df <- fit$df
    lambda.loc <- max(which(df < fixedP))
    s  <-  fit$lambda[ lambda.loc ]
  }
  m <- coef(cv.fit, s=s)[-1,1]
  
  ## get names of predictors for which lasso returns non-zero
  predictors <- which(m != 0)
  names(predictors) <- NULL
  
  return(predictors)
}

# Performs standard least squares and returns
# p-values adjusted for multiple testing (bonferroni)
doLSReg <- function(y, X, predictors, include, intercept) {
  ## magnitude of the prediction subset
  nPredict <- ncol(X)
  
  if(length(predictors) == 0) return( rep(1, nPredict )) # lasso gives no predictors -> return p=1
  
  full.pred <- c(predictors, include) # full set of predictors
  full.pred <- unique(full.pred) # remove any duplicates (unlikely)
  s.mag <- length(full.pred)
  
  fit <- NULL
  
  if(!intercept) {
    fit <- lm(y ~ 0 + X[,full.pred] )
  } else {
    fit <- lm(y ~ X[,full.pred] )
  }  
  ## vector of p values for all predictors to return
  pVec <- rep(1, nPredict)
  
  pvals <- coef(summary(fit))[-1,4]
  
  ## adjust p-values for multiple testing
  pvals <- pvals * s.mag
  pvals <- sapply(pvals, min, 1) # scale > 1 to 1
  
  pVec[predictors] <- pvals
  names(pVec) <- colnames(X)
  return(pVec)
}

## Q(gamma) as defined in eq 2.2
Q.gamma <- function(gam, P) {
  q <- quantile(P / gam, gam)
  names(q) <- NULL
  return(min(q,1))
}

## Implements adaptive gamma search and returns the final P values with
## family-wise error correction (eq. 2.3)
adaptiveP <- function(P, gamma.min) {
  if(is.list(P)) {
    msg <- "P is List! \n"
    msg <- cat(msg, length(P))
    msg <- cat(msg,I)
    stop(msg)
  }
  if(gamma.min == 0) {
    stop("Gamma.min must be greater than 0: about to log it!")
  }
  
  grange <- seq(from=gamma.min, to=0.99, length.out=100)
  Qs <- sapply(grange, Q.gamma, P)
  inf.Q <- min(Qs)
  
  pj <- (1 - log(gamma.min)) * inf.Q
  pj[pj > 1] <- 1 # resize p-values down to 1
  return( pj )
}



