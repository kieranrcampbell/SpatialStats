
######################################
## Utility functions for SpatialPRo ##
######################################

#' This takes the entire set of cells and returns sets whose
#' members aren't neighbours of each other
findDisjointCellSet <- function(sp) {
    nc <- nCells(sp)
    nn <- neighbourIDs(sp)

    x <- allCells <- 1:nc
    listofSets <- list()
    counter <- 1

    while(length(x) > 0) {
        disSet <- removeNeighbours(nn, x)
        listofSets[[ counter ]] <- disSet
        counter <- counter + 1 ;
        x <- setdiff(allCells, unlist(listofSets))
        print(x)
    }

    return(listofSets)
}

removeNeighbours <- function(nn, x) {
    xcopy <- x
    for(i in xcopy) {
        if(i %in% x) {
            x <- setdiff(x, nn[[i]])
        }
    }
    x
}


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
    vifs <- sapply(1:k, function(i) {
        y <- X[,i] ; x <- X[,-i]
        fit <- lm(y ~ x)
        r <- summary(fit)$r.squared
        1/(1-r)
    })
    names(vifs) <- colnames(X)
    vifs
}

#' Takes a matrix of predictors and sequentially removes them until the largest
#' VIF is below a threshold (default=10)
#'
#' @param X Matrix of predictors
#' @param vif.threshold Maximum VIF threshold allowed
#'
#' @export
VIFRemove <- function(X, vif.threshold=10) {
    vif.class <- list()
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


#' Remove predictors based on correlation, and keep a list
#' of groupings
#'
#' @export
corrRemove <- function(X, threshold = 0.5) {
    xcor <- cor(X)
    channel.names <- colnames(X)

    diag(xcor) <- 0 # don't want to remove self interactions

    toRemove <- NULL

    pathwayList <- list()

    while(TRUE) {
        if(max(xcor) > threshold) {
            cn <- colnames(xcor)
            maxloc <- which(xcor == max(xcor), arr.ind=TRUE)[1,] ## find position of maximum
            meanCorr1 <- mean(xcor[maxloc[1],])
            meanCorr2 <- mean(xcor[maxloc[2],])
            r <- maxloc[which.max(c(meanCorr1,meanCorr2))]

            ## pathway stuff
            q <- maxloc[which.min(c(meanCorr1,meanCorr2))]
            pathwayList[[ cn[q] ]] <- cn[r]

            names(r) <- NULL
            toRemove <- c(toRemove,which(channel.names %in% rownames(xcor)[r] ))
            xcor <- xcor[-r,-r]
        } else {
            return(list(toRemove=toRemove,pathways=pathwayList))
        }
    }
}

