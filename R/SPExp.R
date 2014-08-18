## SPExp - Spatial Proteomics Experiment

#' Class SPExp
#'
#' Class \code{SPExp} defines an entire spatial proteomics experiment.
#' It stores a list of the samples (of type \code{SPData}) as well as
#' the original file names and directories.
#'
#' @name SPExp-class
#' @rdname SPExp-class
#' @aliases SPExp
#'
#' @exportClass SPExp
SPExp <- setClass("SPExp",
                  representation = list(dir = "character",
                      files = "character",
                      spdata = "list"))


#' @export
setMethod("show", "SPExp", function(object) {
    cat("An object of class ", class(object), "\n",sep="")
    cat(" Location: ", getDir(object), "\n", sep= "")
    cat( " ", "With ", length(SPlist(object)), " samples \n", sep="")
    cat(" ", IDs(object), "\n", sep=" \n")
    invisible(NULL)
})


#' @rdname getdir-methods
#' @aliases getDir,SPExp-methods
setMethod(f = "getDir",
          signature = "SPExp",
          definition = function(object) object@dir)

#' @rdname files-methods
#' @aliases files,SPExp-methods
setMethod(f = "files",
          signature = "SPExp",
          definition = function(object) object@files)

#' @rdname splist-methods
#' @aliases SPlist,SPExp-methods
setMethod(f = "SPlist",
          signature = "SPExp",
          definition = function(object) object@spdata)

#' @rdname ids-methods
#' @aliases IDs,SPExp-methods
setMethod(f = "IDs",
          signature = "SPExp",
          definition = function(object) sapply(SPlist(object), ID))


#' @export
setMethod("initialize", "SPExp",
          function(.Object, dir, files, spdata) {
              .Object@dir <- dir
              .Object@files <- files
              .Object@spdata <- spdata
              return(.Object)
          })


#' @rdname loadexp-methods
#' @aliases loadExp,SPExp-methods
setMethod("loadExp", "SPExp",
          function(object) {
              object@spdata <- list()
              N <- length(files(object))
              length(object@spdata) <- N
              for(i in 1:N) {
                  object@spdata[[ i ]] <- loadCells(paste(getDir(object),
                                                          files(object)[i],
                                                          sep=""), id=IDs(object)[i])
              }

              return(object)
          })

#' Subset a \code{SPExp} object.
#'
#' Returns a new \code{SPExp} object subsetted using the first \code{i} instances
#'
#' @param x The \code{SPExp} instance to subset
#' @param i The samples to retain
#' @name [
#' @aliases [,SPExp-methods
#' @return A subsetted \code{SPExp} object
#' @rdname extract-methods
#' @exportMethod [
setMethod("[", "SPExp",
          function(x, i) {
              x@files <- files(x)[i]
              x@spdata <- x@spdata[i]
              x@ids <- x@ids[i]
              return(x)
          })

#' Access a sample
#'
#' Extracts a single sample from a \code{SPExp} object
#'
#' @param x The SPExp instance to use
#' @param i The index of the \code{SPData} object within the \code{SPExp} to extract
#'
#' @return The \code{SPData} object in slot \code{i} of the \code{SPExp} instance
#' @name [[
#' @aliases [[,SPExp-methods
#' @rdname dextract-methods
#' @exportMethod [[
setMethod(f = "[[",
          signature = "SPExp",
          function(x,i) {
              return( x@spdata[[i]])
          })


#' Sets a single sample in a \code{SPExp} object
#'
#' @param x The SPExp instance to use
#' @param i The index of the \code{SPData} object within the \code{SPExp} to set
#' @param value An object of class \code{SPData} to replace slot \code{i}
#'
#' @name [[<-
#' @aliases [[,SPExp-methods
#' @rdname dextract-methods
#' @exportMethod [[
setReplaceMethod(f = "[[",
                 signature = "SPExp",
                 definition = function(x,i, value) {
                     x@spdata[[i]] <- value
                     x
                 })

#' Number of samples
#'
#' @name length
#' @param x \code{SPExp instance to use}
#' @return Number of samples in the \code{SPExp} instance
#' @export
setMethod(f = "length",
          signature = "SPExp",
          def = function(x) length(x@spdata))

#' @export
SPExperiment <- function(dir, files, spdata) {
  return ( new("SPExp", dir, files, spdata) )
}

#' Load an experiment from cytobank in matlab format
#'
#' @param directory The directory containing the experiment files
#' @param files The filenames to load. Default is NULL, in which case all files
#' in the directory are used
#' @return A \code{SPExp} object created from the experiment directory.
#' @export
SPExperimentfromDir <- function(directory, files=NULL) {
    if(is.null(files)) files <- dir(directory)

    filesToLoad <- paste(directory,files,sep="/")
    ids <- sapply(filesToLoad, getIDfromTMAname)
    names(ids) <- NULL
    sps <- lapply(1:length(ids), function(i) { loadCells(filesToLoad[i], ids[i])})
    return(SPExperiment(directory, files, sps, ids))
}

getIDfromTMAname <- function(str) {
    splt1 <- strsplit(str, "_ID")[[1]]
    splt2 <- strsplit(splt1[2],"_")[[1]]
    id <- as.numeric(splt2[1])
    id
}



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
