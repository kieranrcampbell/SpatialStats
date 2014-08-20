
#####################################################################
## High dimensional lasso p-vals                                   ##
## Meinshausen, Nicolai, Lukas Meier, and Peter Bühlmann.          ##
## "P-values for high-dimensional regression."                     ##
## Journal of the American Statistical Association 104.488 (2009). ##
##                                                                 ##
## kieran.campbell@sjc.ox.ac.uk                                    ##
#####################################################################

#' High dimensional significance test for LASSO using the multi-sample split method
#' 
#' @references  Meinshausen, Nicolai, Lukas Meier, and Peter Bühlmann.
#' "P-values for high-dimensional regression."
#' Journal of the American Statistical Association 104.488 (2009).
#' 
#' @param y Response vector
#' @param X Design matrix
#' @param B Number of times to partition the sample
#' @param s The value of lambda to use in lasso
#' @param include Set of predictors to be force-included in OLS analysis
#' @param gamma.min
#' @param alpha
#' @param minP
#' 
#' @return A matrix of p-values
#' 
#' @details
#' \itemize{
#'  \item{\code{y} }{Text about y}
#'  \item{\code{X} }{Text about X}
#'  \item{\code{B} }{Text about B}
#'  \item{\code{s} }{Lots of text about s}
#' }
#' 
#' @export
#' @author Kieran Campbell \email{kieran.campbell@@dpag.ox.ac.uk}
#' @examples 
#' \dontrun{
#' pvals <- hdimLasso(y,X)
#' 
#' }
lassoSig <- function(y,X, B=100, s=c("lambda.min","lambda.1se","halfway","usefixed"),
                                      gamma.min=0.05, alpha=0.05, minP=NULL, include) {
    require(glmnet)

    if(nrow(X) != length(y)) stop("Number of samples doesn't match")

    p.mat <- replicate(B, doSingleSplit(y,X,s,alpha, minP, include=NULL) )
    mode(p.mat) <- "numeric" # weird

    adjusted.pvals <- apply(p.mat, 1, adaptiveP, gamma.min)

    return(adjusted.pvals)
}

#' Performs a single split of the data into floor(N/2) and N - floor(N/2)
#' groups and performs lasso & LS estimation
doSingleSplit <- function(y,X,s,alpha, minP, include) {
    n.samp <- length(y)
    sample.in <- sample(1:n.samp, size = floor(n.samp / 2))
    sample.out <- setdiff(1:n.samp, sample.in)

    y.in <- y[sample.in] ; y.out <- y[sample.out]
    X.in <- X[sample.in,] ; X.out <- X[sample.out,]

    predictors <- doLasso(y.in, X.in,s, minP)

    p.vals <- doLSReg(y.out, X.out, predictors, alpha, include)
    return(p.vals)
}

#' Performs lasso on y & X, with custom option of s to be midway point between lambdamin
#' and lambda within 1se
doLasso <- function(y,X,s, minP = NULL) {
    cv.fit <- cv.glmnet(X,y, standardize=FALSE)
    if(s == "halfway") {
        s <- (cv.fit$lambda.min + cv.fit$lambda.1se) / 2
    }

    if(s == "usefixed") {
        if(is.null(minP)) stop("Minimum predictors selected but minP is NULL")
        fit <- cv.fit$glmnet.fit
        df <- fit$df
        lambda.loc <- max(which(df < minP))
        s = fit$lambda[ lambda.loc ]
    }
    m <- coef(cv.fit, s=s)

    ## get names of predictors for which lasso returns non-zero
    predictors <- setdiff(which(m != 0), 1) - 1

    return(predictors)
}

#' Performs standard least squares and returns
#' p-values adjusted for multiple testing (bonferroni)
#' keep.last: we want to force it to take the samples into account
#' so make sure the last four predictors are always in the model
doLSReg <- function(y, X, predictors, alpha, include) {
    ## magnitude of the prediction subset
    s.mag <- length(predictors) + length(include)
    nPredict <- ncol(X)

    if(length(predictors) == 0) return( rep(1, nPredict )) # lasso gives no predictors -> return p=1

    full.pred <- c(predictors, include) # full set of predictors
    full.pred <- unique(full.pred) # remove any duplicates (unlikely)

    fit <- lm(y ~ X[,full.pred] )

    ## vector of p values for all predictors to return
    pVec <- rep(1, nPredict)

    co <- coef(summary(fit))

    ## find significant p-values
    which.signif <- which(co[-1,4] < alpha)
    p.signif <- co[-1,4][which.signif]

    ## adjust p-values for multiple testing
    p.signif <- p.signif * s.mag
    p.signif <- sapply(p.signif, min, 1) # scale > 1 to 1

    pVec[ full.pred[which.signif] ] <- p.signif
    names(pVec) <- colnames(X)
    pVec
}

#' Q(gamma) as defined in eq 2.2
Q.gamma <- function(gam, P) {
    ##print(gam)
    ##print(class(P))
    q <- quantile(P / gam, gam)
    names(q) <- NULL
    return(min(q,1))
}

#' Implements adaptive gamma search and returns the final P values with
#' family-wise error correction (eq. 2.3)
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


