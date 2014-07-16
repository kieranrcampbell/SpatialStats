
#####################################################################
## High dimensional lasso p-vals                                   ##
## Meinshausen, Nicolai, Lukas Meier, and Peter Bühlmann.          ##
## "P-values for high-dimensional regression."                     ##
## Journal of the American Statistical Association 104.488 (2009). ##
##                                                                 ##
## kieran.campbell@sjc.ox.ac.uk                                    ##
#####################################################################

#' @param y Vector of response variables
#' @param X Design matrix
#' @param B Number of times to partition the matrix
#' @param s The value of lambda to use in lasso (normally use
#' either "lambda.min" or "lambda.1se" as per glmnet package)
hdimLasso <- function(y,X, B=100, s="lambda.1se", gamma.min=0.05) {
    require(glmnet)

    p.mat <- replicate(B, doSingleSplit(y,X,s) )

    adjusted.pvals <- apply(p.mat, 1, adaptiveP, gamma.min)

    return(adjusted.pvals)
}

doSingleSplit <- function(y,X,s,alpha=0.05) {
    n.samp <- length(y)
    sample.in <- sample(1:n.samp, size = floor(n.samp / 2))
    sample.out <- setdiff(1:n.samp, sample.in)

    y.in <- y[sample.in] ; y.out <- y[sample.out]
    X.in <- X[sample.in,] ; X.out <- X[sample.out,]

    predictors <- doLasso(y.in, X.in,s)

    print(paste("Using", length(predictors), "out of total", ncol(X), sep=" "))

    p.vals <- doLSReg(y.out, X.out, predictors, alpha)
    return(p.vals)
}

doLasso <- function(y,X,s) {
    cv.fit <- cv.glmnet(X,y)
    if(s == "halfway") {
        s = (cv.fit$lambda.min + cv.fit$lambda.1se) / 2
    }

    m <- coef(cv.fit, s=s)

    ## get names of predictors for which lasso returns non-zero
    predictors <- setdiff(which(m != 0), 1) - 1

    return(predictors)
}

doLSReg <- function(y, X, predictors, alpha) {
    ## magnitude of the prediction subset
    s.mag <- length(predictors)
    nPredict <- ncol(X)

    if(length(predictors) == 0) return( rep(1, nPredict )) # lasso gives no predictors -> return p=1

    fit <- lm(y ~ X[,predictors] )

    ## vector of p values for all predictors to return
    pVec <- rep(1, nPredict)

    co <- coef(summary(fit))

    ## find significant p-values
    which.signif <- which(co[-1,4] < alpha)
    p.signif <- co[-1,4][which.signif]

    ## adjust p-values for multiple testing
    p.signif <- p.signif * nPredict
    p.signif <- sapply(p.signif, min, 1) # scale > 1 to 1

    pVec[ predictors[which.signif] ] <- p.signif
    names(pVec) <- colnames(X)
    pVec
}

Q.gamma <- function(gamma, P) {
    q <- quantile(P / gamma, gamma)
    names(q) <- NULL
    q
}


adaptiveP <- function(P, gamma.min) {
    if(gamma.min == 0) {
        stop("Gamma.min must be greater than 0: about to log it!")
    }

    grange <- seq(from=gamma.min, to=0.99, length.out=100)
    Qs <- sapply(grange, Q.gamma, P)
    inf.Q <- min(Qs)

    return( (1 - log(gamma.min) ) * inf.Q )
}
