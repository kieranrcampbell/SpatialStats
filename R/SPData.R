################################################
## SpatialPRo R class                         ##
## kieran.campbell@dpag.ox.ac.uk              ##
## Data container for spatial proteomics data ##
################################################

#' Class SPdata
#'
#' Class \code{SPData} defines a spatial proteomics sample. It
#' stores a sample (e.g. a tissue microarray) of a given
#' tissue as single cell proteomics data. It contains
#' cell measurements, size and classification.
#'
#' @name SPData-class
#' @rdname SPData-class
#' @aliases SPData
#'
#' @exportClass SPData
SPData <- setClass(Class = "SPData",
                   representation = list(channelNames = "character",
                       readouts = "matrix", ## +min
                       raw = "matrix", ## + min + log
                       cellNeighbours = "list",
                       nn.ids = "list",
                       size = "numeric",
                       id = "numeric",
                       weights = "list",
                       pos = "matrix", # nCell by 2 matrix of cell locations
                       cellClass = "numeric"))

#' @rdname cells-methods
#' @aliases cells,SPData-methods
setMethod(f = "cells",
          signature = "SPData",
          definition = function(object) object@readouts )

#' @rdname cells-methods
#' @name cells<-
#' @aliases cells<-,SPData-methods
setReplaceMethod("cells", signature = "SPData",
                 function(object, value) {
                     object@readouts <- value
                     validObject(object)
                     return(object)
                 })

#' @rdname raw-methods
#' @aliases rawData,SPData-methods
setMethod(f = "rawData",
          signature = "SPData",
          definition = function(object) object@raw)

#' @rdname nCells-methods
#' @aliases nCells,SPData-methods
setMethod(f = "nCells",
          signature = "SPData",
          definition = function(object) dim(object@raw)[1] )

#' @rdname nChannel-methods
#' @aliases nChannel,SPData-methods
setMethod(f = "nChannel",
          signature = "SPData",
          definition = function(object) dim(object@readouts)[2] )

#' @rdname channels-methods
#' @aliases channels,SPData-methods
setMethod(f = "channels",
          signature = "SPData",
          definition = function(object) object@channelNames )

#' @rdname neighbours-methods
#' @aliases neighbours,SPData-methods
setMethod(f = "neighbours",
          signature = "SPData",
          definition = function(object) {
    object@cellNeighbours
})

#' @rdname neighbours-methods
#' @name neighbours<-
#' @aliases neighbours<-,SPData-methods
setReplaceMethod("neighbours", signature = "SPData",
                 function(object, value) {
                     object@cellNeighbours <- value
                     validObject(object)
                     return(object)
                 })


#' @rdname size-methods
#' @aliases size,SPData-methods
setMethod(f = "size",
          signature = "SPData",
          definition = function(object) object@size)

#' @rdname weight-methods
#' @aliases weight,SPData-methods
setMethod(f = "weight",
          signature = "SPData",
          def = function(object) object@weights)

#' @rdname weight-methods
#' @name weight<-
#' @aliases weight<-,SPData-methods
setReplaceMethod(f = "weight",
                 signature="SPData",
                 function(object, value) {
                     object@weights <- value
                     return(object)
                 })


#' Dimension of underlying matrix representation
#'
#' Vector of length 2 that represents the dimensions of the
#' underlying cell matrix. The first entry is the number of cells
#' and the second is the number of channels. Equivalent to
#' dim(cells(x))
#'
#' @param x The SPData object to use
#' @name dim
#' @rdname dim-methods
#' @exportMethod dim
setMethod(f = "dim",
          signature = "SPData",
          def  = function(x) c(nCells(x), nChannel(x)))

#' @rdname id-methods
#' @aliases ID,SPData-methods
setMethod(f = "ID",
          signature = "SPData",
          def  = function(object) object@id)

#' @name ID<-
#' @rdname id-methods
#' @aliases ID<-,SPData-methods
setReplaceMethod(f = "ID",
                 signature = "SPData",
                 function(object, value) {
                     object@id <- value
                     return(object)
                 })

#' @rdname neighbourid-methods
#' @aliases neighbourIDs,SPData-methods
setMethod(f = "neighbourIDs",
          signature = "SPData",
          def = function(object) object@nn.ids)

#' @rdname xy-methods
#' @aliases xy,SPData-methods
setMethod(f = "xy",
          signature = "SPData",
          def = function(object) object@pos)

#' @rdname xy-methods
#' @aliases xy<-,SPData-methods
#' @name xy<-
setReplaceMethod(f = "xy",
                 signature="SPData",
                 function(object, value) object@pos <- xy)



#' Default show call
#' @export
setMethod("show", "SPData", function(object) {
    cat("An object of class ", class(object), "\n",sep="")
    if(ID(object) > -1) cat(" Sample ID: ", ID(object), "\n", sep="")
    cat(" ", nCells(object), " cells with ",
        nChannel(object), " channel(s)\n", sep="")
    invisible(NULL)

})


setValidity("SPData", function(object) {
    msg <- NULL
    valid <- TRUE
    if(length(neighbours(object)) > 1) {
        if(nCells(object) != length(neighbours(object))) {
            valid <- FALSE
            msg <- c(msg, "Nearest neighbour data not available for all cells")
        }
    }

    if(!is.null(xy(object)) && length(xy(object)) > 1) {
        if(nCells(object) != nrow(xy(object))) {
            valid <- FALSE
            msg <- c(msg, "Length mismatch between location info and number of cells")
        }
    }

    if(length(weight(object)) > 1 ) { ## okay for object not to have weights
        if(length(weight(object)) != nCells(object)) {
            valid <- FALSE
            msg <- c(msg, "Number of weights and number of cells differ")
        }
    }

    if(nCells(object) != nrow(rawData(object))) {
        valid <- FALSE
        ##print(nCells(object))
        ##print(nrow(cells(object)))
        msg <- c(msg, "Number of cells must be equal to number of rows in cell by protein matrix")
    }

    if(nChannel(object) != dim(cells(object))[2]) {
        valid <- FALSE
        msg <- c(msg, "Number of proteins must be equal to number of columns in cell by protein matrix")
    }

    if(valid) TRUE else msg

})


#' Subset an SPData set
#'
#' Select SPData[i,j] for cells \code{i} and channels \code{j}.
#' Note that this does not subset out nearest neighbours also.
#'
#' @param i Cells to subset
#' @param j Channels to subset
#' @name [
#'
#' @return An SPData object reduced to cells \code{i} and channels \code{j}
#'
#' @aliases [,SPData-methods
#' @rdname extract-methods
#'
#' @export
setMethod("[", "SPData", function(x, i, j) {
    if(missing(j)) j <- 1:nChannel(x)
    if(missing(i)) i <- 1:nCells(x)

    .n.proteins <- length(j)
    .channelNames <- channels(x)[j]
    .id <- ID(x)
    .weight <- weight(x)[i]
    .nnid <- neighbourIDs(x)[i]
    .cell.class <- cellClass(x)[i]
    .pos <- xy(x)

    .Y <- as.matrix(cells(x)[i,j])
    .raw <- as.matrix(rawData(x)[i,j])

    .X <- lapply(neighbours(x), function(nn.cells) {
        if(is.matrix(nn.cells)) nn.cells[,j] else t(as.matrix(nn.cells[j]))
    })

    .X <- .X[i]

    .size <- size(x)[i]

    SPData(channelNames = .channelNames,
           readouts = .Y,
           cellNeighbours = .X,
           size = .size, id=.id,
           weights = .weight, nn.ids = .nnid,
           raw = .raw, cellClass = .cell.class, pos=.pos)
})

###########################################
## Methods for neighbours and cell class ##
###########################################

#' @rdname neighbourmean-methods
#' @aliases neighbourMean,SPData-methods
setMethod(f = "neighbourMean",
          signature = signature("SPData", "logical", "logical"),
          def = function(object, useWeights, normalise) {
              #if(missing(useWeights)) useWeights <- FALSE
              #if(missing(normalise)) normalise <- TRUE


              ## average over nearest neighbours then means
              X <- neighbours(object)
              weights <- NULL
              if(useWeights) weights <- weight(object)

              X <- lapply(1:length(X), function(i) {
                  nn <- X[[i]]

                  if(is.matrix(nn)) {
                      if(useWeights) {
                          w <- weights[[i]]
                          total.boundary <- sum(w)
                          nn <- nn * w / total.boundary ## IMPORTANT: matrix * vector multiplication is by column
                          return( colSums(nn) )
                      } else {
                          return( colMeans(nn) )
                      }
                  }
                  else {
                      return( nn )
                  }
              })

              X <- matrix(unlist(X), nrow=length(X), byrow=TRUE)

              if(normalise) X <- apply(X, 2, function(x) (x - mean(x))/sd(x))

              colnames(X) <- channels(object)
              X
          })


#' @rdname cellclass-methods
#' @aliases cellClass,SPData-methods
setMethod("cellClass", signature="SPData", function(object) object@cellClass)

#' @name cellClass<-
#' @rdname cellclass-methods
#' @aliases cellClass<-,SPData-methods
setReplaceMethod("cellClass", signature="SPData",
                 function(object, value) {
                     object@cellClass <- value
                     return(object)
                 })

#' @rdname neighbourclass-methods
#' @aliases neighbourClass,SPData-methods
setMethod("neighbourClass", signature("SPData","numeric"),
          function(object, cell.class) {
              X <- neighbours(sp)
              nn.ids <- neighbourIDs(sp)
              cell.select <- which(cellClass(sp) == cell.class)

              nn <- lapply(1:length(X), function(i) {
                  Xi <- X[[i]] ; ids <- nn.ids[[i]]

                  lvec <- ids %in% cell.select
                  if(is.matrix(Xi)) {
                      if(any(lvec)) return(Xi[lvec,]) else return(numeric(0))
                  } else {
                      if(lvec) return(Xi) else return(numeric(0))
                  }
              })
              nn
          })

#' Selects only particular channels to be returned as nearest neighbours
#'
#' @param NN Neighbour list (e.g. via neighbours(sp))
#' @param channel.ids A channel list to select out
#'
#' @export
neighbourChannel <- function(NN, channel.ids) {
    nn <- lapply(NN, function(Xi) {
        if(is.matrix(Xi)) {
            Xi[,channel.ids]
        } else {
            if(length(Xi) > 0) Xi[channel.ids] else numeric(0)
        }
    })
    nn
}

#' Cells on a class boundary
#'
#' Find the cells that lie on the boundary between two classes
#' (currently only implemented for 2 classes)
#'
#' @param sp The SPData object to use
#'
#' @return A vector of cell identifiers relating to those
#' that lie along the boundary
#'
#' @export
findBoundary <- function(sp) {
    classes <- cellClass(sp)
    cl1 <- which(classes == 1) ; cl2 <- which(classes == 2)
    nn.ids <- neighbourIDs(sp)

    cell1neighbours <- sapply(cl1, function(cellid) {
        nn.id <- nn.ids[[cellid]]
        any(nn.id %in% cl2)
    })

    cell2neighbours <- sapply(cl2, function(cellid) {
        nn.id <- nn.ids[[cellid]]
        any(nn.id %in% cl1)
    })
    boundary <- sort(c(cl1[cell1neighbours], cl2[cell2neighbours]))
}





