## spatialpro R class
## kieranrcampbell@gmail.com

SpatialPRo <- setClass("SpatialPRo",
                       slots = c(n.cells = "integer",
                           n.proteins = "integer",
                           protein.names = "character",
                           Y = "matrix",
                           X = "list"))

#' Extracts the cell proteomics data
#' @export
setMethod("cells", "SpatialPRo", function(object) object@Y )

#' Returns the number of cells in the sample
#' @export
setMethod("nCells", "SpatialPRo", function(object) object@n.cells )

#' Returns the number of proteins measured (number of channels)
#' @export
setMethod("nProt", "SpatialPRo", function(object) object@n.proteins )

#' Returns the names of the proteins measured
#' @export
setMethod("pNames", "SpatialPRo", function(object) object@protein.names )

#' Returns a list of nearest neighbour readouts
#'
#' The ith entry is an n by m matrix, for cell i having n neighbours
#' each of which have m channels
#' @export
setMethod("NN", "SpatialPRo", function(object) object@X )

#' @export
setMethod("show", "SpatialPRo", function(object) {
    cat("An object of class ", class(object), "\n",sep="")
    cat(" ", nCells(object), " cells with ", nProt(object), " proteins\n", sep="")
    invisible(NULL)

})


setValidity("SpatialPRo", function(object) {
    msg <- NULL
    valid <- TRUE
    if(nCells(object) != length(NN(object))) {
        valid <- FALSE
        msg <- c(msg, "Nearest neighbour data not available for all cells")
    }

    if(nCells(object) != dim(cells(object))[1]) {
        valid <- FALSE
        msg <- c(msg, "Number of cells must be equal to number of rows in cell by protein matrix")
    }

    if(nProt(object) != dim(cells(object))[2]) {
        valid <- FALSE
        msg <- c(msg, "Number of proteins must be equal to number of columns in cell by protein matrix")
    }

    if(valid) TRUE else msg

})

#' @export
setMethod("[", "SpatialPRo", function(x,i,drop="missing") {
    .n.proteins <- length(i)
    .protein.names <- pNames(x)[i]
    .Y <- cells(x)[,i]
    .X <- lapply(NN(x), function(nn.cells) {
        if(is.matrix(nn.cells)) nn.cells[,i] else nn.cells[i]
    })

    SpatialPRo(n.cells = nCells(x),
               n.proteins = .n.proteins,
               protein.names = .protein.names,
               Y = .Y,
               X = .X)
})
