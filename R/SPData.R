## spatialpro R class
## kieranrcampbell@gmail.com

#' Class containing spatial proteomics data
#'
#' @export
SPData <- setClass("SPData",
                   representation = list(channelNames = "character",
                       readouts = "matrix", ## +min + loess & N(0,1) normalised
                       raw = "matrix", ## + min + log
                       cellNeighbours = "list",
                       nn.ids = "list",
                       size = "numeric",
                       id = "numeric",
                       weights = "list",
                       pos = "numeric", ## not yet implemented
                       cellClass = "numeric"))

#' Extracts the cell proteomics data
#' @export
setMethod("cells", "SPData", function(object) object@readouts )

#' Extract the raw (logged) data
#' @export
setMethod("rawData", "SPData", function(object) object@raw)

#' Returns the number of cells in the sample
#' @export
setMethod("nCells", "SPData", function(object) dim(object@readouts)[1] )

#' Returns the number of proteins measured (number of channels)
#' @export
setMethod("nChannel", "SPData", function(object) dim(object@readouts)[2] )

#' Returns the names of the proteins measured
#' @export
setMethod("channels", "SPData", function(object) object@channelNames )

#' Returns a list of nearest neighbour readouts
#'
#' The ith entry is an n by m matrix, for cell i having n neighbours
#' each of which have m channels
#' @export
setMethod("neighbours", "SPData", function(object) {
    #nn <- lapply(object@nn.ids, function(nn.id) { object@readouts[nn.id,] })
    #return(nn)
    object@cellNeighbours
})

#' Returns the cell sizes
#'
#' The ith entry is the size of the ith cell, as ordered by cells(X)
#' @export
setMethod("size", "SPData", function(object) object@size)

#' Returns a list of nearest neighbour boundary sizes (weights)
#'
#' @export
setMethod("weight", "SPData", function(object) object@weights)

#' Sets the boundary weights
#'
#' @name weight<-
#' @export
setReplaceMethod("weight", signature="SPData",
                 function(object, value) {
                     object@weights <- value
                     return(object)
                 })


#' Gives dimension of underlying matrix representation
#'
#' Returns number of cells and number of proteins as dimension of
#' underlying matrix
#' @export
setMethod("dim", "SPData", function(x) c(nCells(x), nChannel(x)))

#' Returns the sample id
#'
#' @export
setMethod("id", "SPData", function(object) object@id)

#' Sets the sample id
#'
#' @name id<-
#' @export
setReplaceMethod("id", signature = "SPData",
                 function(object, value) {
                     object@id <- value
                     return(object)
                 })

#' Returns the nearest neighbour ids
#'
#' @export
setMethod("neighbourIDs", "SPData", function(object) object@nn.ids)

#' @export
setMethod("show", "SPData", function(object) {
    cat("An object of class ", class(object), "\n",sep="")
    if(id(object) > -1) cat(" Sample ID: ", id(object), "\n", sep="")
    cat(" ", nCells(object), " cells with ", nChannel(object), " channel(s)\n", sep="")
    invisible(NULL)

})


setValidity("SPData", function(object) {
    msg <- NULL
    valid <- TRUE
    if(nCells(object) != length(neighbours(object))) {
        print(nCells(object))
        valid <- FALSE
        msg <- c(msg, "Nearest neighbour data not available for all cells")
    }

    if(length(weight(object)) > 0 ) { ## okay for object not to have weights
        if(length(weight(object)) != nCells(object)) {
            valid <- FALSE
            msg <- c(msg, "Number of weights and number of cells differ")
        }
    }

    if(nCells(object) != dim(cells(object))[1]) {
        valid <- FALSE
        msg <- c(msg, "Number of cells must be equal to number of rows in cell by protein matrix")
    }

    if(nChannel(object) != dim(cells(object))[2]) {
        valid <- FALSE
        msg <- c(msg, "Number of proteins must be equal to number of columns in cell by protein matrix")
    }

    if(valid) TRUE else msg

})

#' Sets the cell by cell measurements
#'
#' @name cells<-
#' @export
setReplaceMethod("cells", signature = "SPData",
                 function(object, value) {
                     object@readouts <- value
                     validObject(object)
                     return(object)
                 })

#' Sets the neighbour measurements
#'
#' @name neighbours<-
setReplaceMethod("neighbours", signature = "SPData",
                 function(object, value) {
                     object@cellNeighbours <- value
                     validObject(object)
                     return(object)
                 })


#' Subset an SPData set
#'
#' Select SPData[i,j] for cells i and channels j.
#' @export
setMethod("[", "SPData", function(x, i, j) {
    if(missing(j)) j <- 1:nChannel(x)
    if(missing(i)) i <- 1:nCells(x)

    .n.proteins <- length(j)
    .channelNames <- channels(x)[j]
    .id <- id(x)
    .weight <- weight(x)[i]
    .nnid <- neighbourIDs(x)[i]
    .cell.class <- cellClass(x)[i]

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
           raw = .raw, cellClass = .cell.class)
})


#############################################
## Methods for nearest neighbour averaging ##
#############################################

#' Averages over the nearest neighbour cells, going from a
#' list of length cell to a cell by channel matrix
#'
#' @param useWeights If TRUE then nearest neighbours are weighted by cell boundary size
#' @param normalise If TRUE then then each channel is normalised to mean 0 sd 1
#'
#' @export
setMethod("neighbourMean", signature("SPData", "logical", "logical"),
          function(object, useWeights, normalise) {
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
                          colSums(nn)
                      } else {
                          colMeans(nn)
                      }
                  }
                  else {
                      nn
                  }
              })

              X <- matrix(unlist(X), nrow=length(X), byrow=TRUE)

              if(normalise) X <- apply(X, 2, function(x) (x - mean(x))/sd(x))

              colnames(X) <- channels(object)
              X
          })


###########################################
## Methods for dealing with cell class & ##
##     boundaries start here             ##
###########################################


#' If each cell has a class (e.g. tumour or stromal) then it can be
#' assigned a numeric class which is retrieved through cellClass
#'
#' @export
setMethod("cellClass", signature="SPData", function(object) object@cellClass)

#' Sets the cell class
#'
#' @name cellClass<-
#' @export
setReplaceMethod("cellClass", signature="SPData",
                 function(object, value) {
                     object@cellClass <- value
                     return(object)
                 })

#' Filter out nearest neighbours by cell class. If cell i has no
#' nearest neighbours of class cell.class then numeric(0) is returned.
#'
#' @export
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

#' Filter out nearest neighbours by cell class. If cell i has no
#' nearest neighbours of class cell.class then numeric(0) is returned.
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

#' Find the cells that lie on the boundary between two classes
#' (currently only implemented for 2 classes)
#'
#' @export
setMethod("findBoundary", signature("SPData"),
          function(object) {
              classes <- cellClass(object)
              cl1 <- which(classes == 1) ; cl2 <- which(classes == 2)
              nn.ids <- neighbourIDs(object)

              cell1neighbours <- sapply(cl1, function(cellid) {
                  nn.id <- nn.ids[[cellid]]
                  any(nn.id %in% cl2)
              })

              cell2neighbours <- sapply(cl2, function(cellid) {
                  nn.id <- nn.ids[[cellid]]
                  any(nn.id %in% cl1)
              })

              boundary <- sort(c(cl1[cell1neighbours], cl2[cell2neighbours]))
          })

#################################################
## plotting & Visualisation methods start here ##
#################################################

#' Boxplots the distribution for each channel
#'
#' @param nrow Number of rows of the boxplots to plot
#' @param ncol Number of columns of boxplots to plot
#'
#' @export
setMethod("boxplots", signature("SPData", "numeric", "numeric"),
          function(object, nrow, ncol) {
              readouts <- cells(object)
              read.boundaries <- quantile(readouts, probs=c(0.01, 0.99))
              if(nrow * ncol > nChannel(object)) {
                  print("Not enough channels to boxplot")
                  return(NULL)
              } else {
                  par(mfrow=c(nrow,ncol), mar=c(1.8,1.8,1.8,1.8))
                  for(i in 1:(nrow*ncol)) {
                      boxplot(readouts[,i], ylim=read.boundaries,
                              main=channels(object)[i])
                  }
              }
          })

#' Boxplots the distributions of channels depending on their class
#'
#' @param channel.ids A numeric vector of channel indicies
#' @export
setMethod("channelPlot", signature("SPData", "numeric"),
          function(object, channel.ids) {
              require(ggplot2)
              require(reshape)
              readouts <- cells(object)
              classes <- as.factor(cellClass(object))
              readouts <- readouts[,channel.ids]
              colnames(readouts) <- channels(object)[channel.ids]
              readouts.melted <- melt(readouts)
              plotdf <- data.frame(exprs=readouts.melted$value,
                                   channel=readouts.melted$X2,
                                   cell.class = rep(classes,length(channel.ids)))
              ggplot(aes(x=channel,y=exprs,fill=cell.class), data=plotdf) + geom_boxplot() +
                  theme_bw() + theme(axis.text.x = element_text(angle=90, hjust=1))
          })

#' Correlation plot of channels within the sample
#'
#' @export
channelCorr <- function(sp, cell.class=NULL) {
    require(corrplot)
    Y <- cells(sp)
    if(!is.null(cell.class)) {
        Y <- Y[cellClass(sp) == cell.class,]
    }
    corr.y <- cor(Y)
    corrplot(corr.y)
}

#' Correlation plot of channels with neighbours
#'
#' @export
neighbourCorr <- function(sp, cell.class=NULL) {
    require(corrplot)
    Y <- cells(sp) ; X <- neighbourMean(sp, TRUE, TRUE)
    if(!is.null(cell.class)) {
        Y <- Y[cellClass(sp) == cell.class,]
        X <- X[cellClass(sp) == cell.class,]
    }
    xy.corr <- cor(X,Y)
    corrplot(xy.corr)
}


#######################################
## Methods for importing from matlab ##
## files start here                  ##
#######################################


#' Loads an Xell matlab file into the SPData format
#'
#' This function parses the matlab files, pulling out relevant proteins
#' in the 'D' channel and validates that the correct proteins are present. Matlab
#' files can be found at
#' https://s3.amazonaws.com/supplemental.cytobank.org/report_data/report_113/Figure_5/Figure_5_raw_image_files.zip
#'
#' @param filename The matlab file
#' @param control.isotopes The isotopes used for control to exclude from analysis
#' @param log.data Boolean of whether to log the data
#' @export
loadCells <- function(filename, id=-1, control.isotopes = c("Xe131","Cs133","Ir193")) {
    require(R.matlab)

    ## loads relevant data from matlab and parses into list

    m <- readMat(filename)

    n.cells <- dim(m$Xell)[1]
    n.channels <- -1


    ycolheads <- as.character(m$Xell.list.col)
    yp.id <- grep(")D", ycolheads)
    yp.id <- setdiff(yp.id, grep(paste(control.isotopes, collapse="|"),ycolheads))
    channelNames <- ycolheads[yp.id]

    Y <- m$Xell.list
    Y <- Y[,yp.id]

    colnames(Y) <- channelNames

    ## add minimum to make all values positive
    Y <- preprocess.addmin(Y)

    ## fork off 'raw' that this point
    raw <- Y # log(Y)

    ## want to LOESS normalise against cell size
    sizes <- as.numeric(m$Xell.size)
    ##Y <- loessNormalise(Y, sizes)
    ##Y <- totalProteinNormalise(Y)


    ## now on to constructing X, the nearest neighbour matrix
    ## nnids list of nearest neighbour IDs
    nnids <- lapply(m$Xell.nearest, function(xl) {
        if(is.matrix(xl)) {
            return ( xl[,1] )
        } else {
            xl[1]
        }
    })


    ## Y <- preprocess.centre(Y) # don't want to centre Y until after finding X

    sp <- SPData(channelNames=channelNames,
                 readouts=matrix(0), raw=raw, cellNeighbours=list(0),
                 size=sizes,id=id,
                 nn.ids=nnids)
    return( sp )
}

preprocess.addmin <- function(Y) {
    mu.bg <- -min(Y)
    Y <- Y + mu.bg + 1
}

## converts the measurements for each channel
## to N(0,1)
preprocess.centre <- function(Y) {
    Y <- apply(Y, 2, function(y) (y - mean(y)) / sd(y))
    Y
}


