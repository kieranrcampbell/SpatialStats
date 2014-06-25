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


#' Get p-value matrix from fit
#'
#' This function constructs a nProt by nProt directed
#' matrix of p values
#'
#' @param spdata The SPdata object the fit was made from
#' @param fit The fit returned by calling nnReg
#' @param setZero Sets all entries to zero that don't pass the threshold
#' @param alpha The significance level for the tests.
#'
#' @export
getPMat <- function(spdata, fit, setOne=FALSE, alpha=0.05) {
    ## returns an adjacency matrix from a general linear model fit
    s <- summary(fit)

    n.proteins <- nProt(spdata)

    A <- matrix(0,nrow=n.proteins,ncol=n.proteins)

    for(i in 1:n.proteins) {
        ## looking at which proteins influence i
        coeff <- s[[i]]$coefficients
        p.vals <- as.numeric(coeff[,4])
        p.vals <- p.vals[-1] ## remove intercept consideration
        pcol <- p.vals
        if(setOne) {
            p.vals <- p.vals < alpha
            pcol[!p.vals] <- 1
        }
        A[,i] <- pcol
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

#' Performs weighted (cell edge weights) regression along the tumour-stromal
#' boundary for a given set of response channels and a class (tumour or stromal)
#' to focus on in the reponses
#' @param sp The SPData class for the experiment
#' @param cellClasses A vector of length nCells(sp) where each entry is either 1 or 2
#' @param responseClass A numeric 1 or 2 indicating which cells should be included in the response variable
#' @param boundary A vector with the indices of cells that lie along the beoundary
#' @param weights A list of length nCells where each entry gives the length of the boundary between
#' it and all nearest neighbour cells
#' @param responseSubset A vector containing the indices of interesting proteins to use in the response variables
#' @param dependentSubset A vector containing the indices of proteins to include in regression coefficients
#'
#' @export
weightedSubsetBoundaryRegression <- function(sp, cellClasses, responseClass, boundary,
                                             weights, responseSubset, dependentSubset) {
    cl1 <- which(cellClasses == 1) ; cl2 <- which(cellClasses == 2)

    cellBoundary <- list()
    cellBoundary[[1]] <- intersect(cl1, boundary)
    cellBoundary[[2]] <- intersect(cl2, boundary)
    print(length(cellBoundary[[1]]))
    print(length(cellBoundary[[2]]))

    X <- weightNN(NN(sp), weights)
    ## list involving only nn of class 1
    X.class1 <- separateNN(X, nnID(sp), cl1)

    ## list involving only nn of class 2
    X.class2 <- separateNN(X, nnID(sp), cl2)

    nProtein <- length(pNames(sp))

    X1 <- sumNN(X.class1, nProtein)
    X2 <- sumNN(X.class2, nProtein)

    #X1 <- t(X1) ; X2 <- t(X2)
    ## select out response variables of interest
    Y <- cells(sp)

    ## ## convert to standard normal
    ## Y <- Y - rowMeans(Y)
    ## sds <- apply(Y, 1, sd)
    ## Y <- Y / sds

    Y <- as.matrix(Y[,responseSubset])

    ## select out only cells along the boundary
    Y <- Y[cellBoundary[[ responseClass ]], ]
    X1 <- X1[cellBoundary[[ responseClass ]], dependentSubset]
    X2 <- X2[cellBoundary[[ responseClass ]], dependentSubset]

    X <- NULL
    if(responseClass == 1) X <- X2
    if(responseClass == 2) X <- X1

    fit <- lm(Y ~ log(X))

    return(list(fit=fit,Y=Y,X=X))
}

#' Separates out the nearest neighbours into only those
#' of a desired class, given the list X, list of nearest
#' neighbour ids nn.ids and the list of indices in desired
#' class (classList)
separateNN <- function(X, nnIds, classList) {
    X.new <- lapply(1:length(X), function(i) {
        x <- X[[i]]
        nnid <- nnIds[[i]]
        if(is.matrix(x)) {
            x <- x[nnid %in% classList, ]
            if(any(nnid %in% classList)) x else  numeric(0)
        } else {
            if(nnid %in% classList) x else  numeric(0)
        }
    })
    return(X.new)
}


#' Takes the nearest neighbour list and returns the weighted version across the cells
weightNN <- function(X, weights) {
    X.w <- lapply(1:length(X), function(i) {
        x <- X[[i]] ; weight <- weights[[i]]

        if(!is.matrix(x)) {
            x
        } else {
            totalBoundary <- sum(weight)
            return(  x * weight / totalBoundary )
        }
    })
    return(X.w)
}

#' Takes a nearest neighbour list (X) and sums the columns
#' (this assumes they're already properly weighted)
sumNN <- function(X, nprotein) {
    t(sapply(X, function(x) {
        if(!is.matrix(x)) {
            if(length(x) == 0) rep(0,nprotein) else x
        } else {
            cs <- colSums(x)
            cs
        }
    }))
}
