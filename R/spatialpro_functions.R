
#' Loads an Xell matlab file into the SpatialPRo formal
#'
#' This function parses the matlab files, pulling out relevant proteins
#' in the 'D' channel and validates that the correct proteins are present.
#'
#' @param filename The matlab file
#' @param control.isotopes The isotopes used for control to exclude from analysis
#' @param log.data Boolean of whether to log the data
#' @export
loadCells <- function(filename, control.isotopes = c("Xe131","Cs133","Ir193"),
                        log.data=TRUE) {
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

    if(log.data) Y <- log(Y)

    get.nn.count <- function(x,p.id) {
        ## returns the nearest neighbour count for proteins
        if(dim(x)[1] == 1) return(x[p.id])
        x <- x[,p.id]
        x
    }

    ## get nearest count data for ind
    nn.count <- m$Xell.nearest

    X <- lapply(nn.count, get.nn.count, xp.id)

    if(log.data) X <- lapply(X, log)

    sp <- SpatialPRo(n.cells=n.cells, n.proteins=n.proteins,protein.names=protein.names,
                     Y=Y, X=X)
    return( sp )
}
