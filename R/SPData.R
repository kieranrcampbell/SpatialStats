## spatialpro R class
## kieranrcampbell@gmail.com

#' Class containing spatial proteomics data
#'
#' @export
SPData <- setClass("SPData",
                   representation = list(names = "character",
                       readouts = "matrix",
                       cellNeighbours = "list",
                       nn.ids = "list",
                       size = "numeric",
                       id = "numeric",
                       weights = "list",
                       pos = "numeric"))

#' Extracts the cell proteomics data
#' @export
setMethod("cells", "SPData", function(object) object@readouts )

#' Returns the number of cells in the sample
#' @export
setMethod("nCells", "SPData", function(object) dim(object@readouts)[1] )

#' Returns the number of proteins measured (number of channels)
#' @export
setMethod("nChannel", "SPData", function(object) dim(object@readouts)[2] )

#' Returns the names of the proteins measured
#' @export
setMethod("channels", "SPData", function(object) object@names )

#' Returns a list of nearest neighbour readouts
#'
#' The ith entry is an n by m matrix, for cell i having n neighbours
#' each of which have m channels
#' @export
setMethod("neighbours", "SPData", function(object) {
    #nn <- lapply(object@nn.ids, function(nn.id) { object@Y[nn.id,] })
    #return(nn)
    object@cellNeighbours
})
    
#' Returns the cell sizes
#'
#' The ith entry is the size of the ith cell, as ordered by cells(X)
#' @export
setMethod("size", "SPData", function(object) object@size)

setMethod("weights", "SPData", function(object) object@weights)

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
                     object@Y <- value
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
    .protein.names <- channels(x)[j]
    .Y <- NULL

    .Y <- as.matrix(cells(x)[i,j])

    .X <- lapply(NN(x), function(nn.cells) {
        if(is.matrix(nn.cells)) nn.cells[,j] else t(as.matrix(nn.cells[j]))
    })

    .X <- .X[i]

    .size <- size(x)[i]

    SPData(protein.names = .protein.names,
           Y = .Y,
           X = .X,
           size = .size)
})



#' Loads an Xell matlab file into the SPData format
#'
#' This function parses the matlab files, pulling out relevant proteins
#' in the 'D' channel and validates that the correct proteins are present.
#'
#' @param filename The matlab file
#' @param control.isotopes The isotopes used for control to exclude from analysis
#' @param log.data Boolean of whether to log the data
#' @export
loadCells <- function(filename, id=-1, control.isotopes = c("Xe131","Cs133","Ir193"),
                        rescale.data=FALSE, scale.factor=10000) {
    library(R.matlab)


    ## loads relevant data from matlab and parses into list

    m <- readMat(filename)

    n.cells <- dim(m$Xell)[1]
    n.proteins <- -1

    ## xcolheads shows the data has the structure of a time cell, then three cells
    ## for each protein (here we use the D channel). The first two and last measurement
    ## are control isotopes, which we disregard.
    xcolheads <- as.character(m$Xell.nearest.col)

    xp.id <- grep(")D",xcolheads) ## x protein ids

    ## remove control isotopes
    xp.id <- setdiff(xp.id, grep(paste(control.isotopes, collapse="|"),xcolheads))
    xprotein.names <- xcolheads[xp.id]
    n.proteins <- length(xprotein.names)

    ## now get cell-by-cell protein data

    ycolheads <- as.character(m$Xell.list.col)
    yp.id <- grep(")D", ycolheads)
    yp.id <- setdiff(yp.id, grep(paste(control.isotopes, collapse="|"),ycolheads))
    yprotein.names <- ycolheads[yp.id]

    if(!all.equal(yprotein.names,xprotein.names)) {
        stop("Mismatch in X and Y protein names")
    }
    protein.names <- xprotein.names

    Y <- m$Xell.list
    Y <- Y[,yp.id]


    get.nn.count <- function(x,p.id) {
        ## returns the nearest neighbour count for proteins
        if(dim(x)[1] == 1) return(x[p.id])
        x <- x[,p.id]
        x
    }

    ## get nearest count data for ind
    nn.count <- m$Xell.nearest

    X <- lapply(nn.count, get.nn.count, xp.id)

    nnids <- lapply(m$Xell.nearest, function(xl) {
        if(is.matrix(xl)) {
            return ( xl[,1] )
        } else {
            xl[1]
        }
    })

    sp <- SPData(protein.names=protein.names,
                     Y=Y, X=X, size=as.numeric(m$Xell.size),id=id,
                 nn.ids=nnids)
    sp <- preprocess(sp, rescale.data, scale.factor)
    return( sp )
}

preprocess <- function(sp, scale.data=FALSE, scale.factor=10000) {
    Y <- cells(sp)
    X <- NN(sp)

    mu.bg <- -min(Y)

    print("Adding -background")
    Y <- log(Y + mu.bg + 1)

    X <- lapply(X, function(x) {
        x + mu.bg + 1
    })

    if(scale.data) {
        print("Rescaling data...")
        Y <- Y / scale.factor
        X <- lapply(X, function(x) x / scale.factor )
    }

    cells(sp) <- Y ; NN(sp) <- X

    return( sp )
}

#' Ridge regression
#' @export
ridgeReg <- function(sp) {
    Y <- cells(sp)

    X <- t(sapply(NN(sp), function(x) {
        if(!is.matrix(x)) x
        else colMeans(x)
    }))
}
