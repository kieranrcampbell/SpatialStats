## regression models for nearest neighbour data
## kieranrcampbell@gmail.com

#' Basic nearest neighbour regression
#' @export
setMethod("nnReg", "SPData", function(object) {
    Y <- cells(object)

    X <- t(sapply(NN(object), function(x) {
        if(!is.matrix(x)) x
        else colMeans(x)
    }))

    fit <- lm(Y ~ X)
    fit
})

#' Get adjacency matrix from fit
#'
#' This function constructs a nProt by nProt directed
#' adjacency matrix from the fit
#'
#' @param spdata The SPdata object the fit was made from
#' @param fit The fit returned by calling nnReg
#' @param alpha The significance level for the tests.
#'
#' @export
getAdj <- function(spdata, fit, alpha=0.05) {
    ## returns an adjacency matrix from a general linear model fit
    s <- summary(fit)

    n.proteins <- nProt(spdata)

    A <- matrix(0,nrow=n.proteins,ncol=n.proteins)

    for(i in 1:n.proteins) {
        ## looking at which proteins influence i
        coeff <- s[[i]]$coefficients
        p.vals <- as.numeric(coeff[,4])
        p.vals <- p.vals[-1] ## remove intercept consideration
        A[,i] <- as.numeric(p.vals < alpha)
    }

    rownames(A) <- colnames(A)  <- pNames(spdata)
    return ( A )
}

reweightReg <- function(spdata, A) {
    fits <- list()
    length(fits) <- nProt(spdata)

    for(p in 1:nProt(spdata)) {
        Y <- as.numeric(cells(spdata[,p]))

        ## which independent variables are we regressing against?
        regressors <- which(A[,p] == 1)

        ## sometimes no significant interactions are found, in which case
        ## we just return NULL
        if(length(regressors) > 0) {

            ## reduced NN matrix - only use significant regression coefficients
            red.nn <- spdata[,as.logical(A[,p])]

            X <- NULL
            ## need to treat n = 1 and n > 1 proteins differently

            if(nProt(red.nn) > 1) {
                X <- t(sapply(NN(red.nn), colMeans))
            } else {
                X <- as.matrix(sapply(NN(red.nn), mean))
                                        #X <- t(X)
            }

            if(dim(X)[1] != length(Y)) {
                stop("dims don't match!!")
            }
            fit <- lm(Y ~ X)
            fits[[p]] <- fit
        }

    }
    return(fits)
}

