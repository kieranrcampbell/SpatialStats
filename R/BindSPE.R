
#' Bind multiple sample data together
#'
#' Given an SPExp object, separate out each cell type,
#' calculate neighbour means and bind together, providing
#' more statistical power. Returns a list with cell readouts (Y),
#' neighbour means (X) and the number of cells from each sample (sizes)
#'
#' @param SPE The SPExp object to use
#' @param pickTumour Logical indicating whether or not to separate out cell types
#' @param useWeights Passed to neighbourMeans
#' @param normalise Whether to centre the columns
#'
#' @export
BindSPE <- function(SPE, pickTumour=TRUE, useWeights=TRUE, normalise=TRUE) {
    ## variable selection
    XY <- lapply(SPlist(SPE), function(sp) {
        Y <- cells(sp)
        X <- neighbourMean(sp, useWeights, normalise)

        if(pickTumour) {
            tumourID <- findTumourID(sp)
            tumourCells <- which(cellClass(sp) == tumourID)
            Y <- Y[tumourCells, ]
            X <- X[tumourCells, ]
        }

        list(X=X,Y=Y)
    })
    sizes <- sapply(XY, function(xy) nrow(xy$Y))

    all.x <- lapply(XY, function(xy) xy$X)
    all.y <- lapply(XY, function(xy) xy$Y)

    X <- do.call(rbind, all.x)
    Y <- do.call(rbind, all.y)

    ## normalise Y & X to have mean 0 and unit variance
    if(normalise) {
        m0uv <- function(x) (x - mean(x)) / sd(x)
        Y <- apply(Y, 2, m0uv)
        X <- apply(X, 2, m0uv)
    }

    return( list( X=X, Y=Y, sizes=sizes))
}

#' Construct model matrix for integrating multiple samples
#'
#' Given the output of the function BindSPE we would like to control
#' for biases introduced by binning all the data together. ConstructSampleFactors
#' takes the output of BindSPE and returns a model matrix of factors highlighting
#' what sample a given cell is from
#'
#' @param XY Output of BindSPE
#'
#' @export
ConstructSampleFactors <- function(XY, sample.ids) {
    cell.sizes <- XY$sizes
    Nfactors <- length(cell.sizes) - 1
    factors <- NULL
    for(i in 1:Nfactors) {
        tcol <- rep(0, sum(cell.sizes))
        cs <- cumsum(cell.sizes)
        range <- (cs[i] + 1):cs[i+1]
        tcol[ range ] <- 1
        factors <- cbind(factors, tcol)
    }
    colnames(factors) <- paste("sample", sample.ids[1:Nfactors], sep="")
    factors
}


#' SPExp object with 5 samples used in report.
#'
#' A dataset of 5 samples of type SPData held in an SPExp object.
#'
#' @docType data
#' @keywords SPExp
#' @name SPE_report
#' @usage data(SPE_report)
#' @format An SPExp object with 5 samples
NULL
