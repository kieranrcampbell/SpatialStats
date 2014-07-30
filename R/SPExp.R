## SPExp - Spatial Proteomics Experiment

SPExp <- setClass("SPExp",
                  representation = list(dir = "character",
                      files = "character",
                      spdata = "list",
                      ids = "numeric"))


#' @export
setMethod("show", "SPExp", function(object) {
    cat("An object of class ", class(object), "\n",sep="")
    cat(" Location: ", getDir(object), "\n", sep= "")
    cat( " ", "With ", length(SPlist(object)), " samples \n", sep="")
    cat(" ", IDs(object), "\n", sep=" \n")
    invisible(NULL)
})

#' Show files
#'
#' @export
setMethod("files", "SPExp", function(object) object@files)

#' Displays the directory of the experiment
#' @export
setMethod("getDir", "SPExp", function(object) object@dir)

#' @export
setMethod("initialize", "SPExp",
          function(.Object, dir, files, spdata, ids) {
              .Object@dir <- dir
              .Object@files <- files
              .Object@spdata <- spdata
              .Object@ids <- ids

              #s <- strsplit(files(.Object),"_ID")
              #.Object@ids <- as.numeric(sapply(s, function(x) strsplit(x[2],"_")[[1]][1]))

              return(.Object)
          })

#' Returns the sample ids.
#'
#' @export
setMethod("IDs", "SPExp",
          function(object) object@ids)

#' Loads in the experiment data
#'
#' @export
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

#' @export
setMethod("[", "SPExp",
          function(x, i) {
              x@files <- files(x)[i]
              x@spdata <- x@spdata[i]
              x@ids <- x@ids[i]
              return(x)
          })
#' @export
setMethod("[[", "SPExp",
          function(x,i) {
              return( x@spdata[[i]])
          })


#' Returns the list of SPData samples
#' @export
setMethod("SPlist", "SPExp",
          function(object) object@spdata)

#' Returns the number of samples
#' @export
setMethod("length", "SPExp", function(x) length(x@spdata))

#' @export
SPExperiment <- function(dir, files, spdata, ids) {
    return ( new("SPExp", dir, files, spdata, ids) )
}

#' Load an experiment from cytobank in matlab format
#'
#' @param directory The directory containing the experiment files
#' @param files The filenames to load. Default is NULL, in which case all files
#' in the directory are used
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
bindSPE <- function(SPE, pickTumour=TRUE, useWeights=TRUE, normalise=TRUE) {
    ## variable selection
    XY <- lapply(SPlist(SPE), function(sp) {
        if(pickTumour) {
            tumourID <- findTumourID(sp)
            tumourCells <- which(cellClass(sp) == tumourID)
        }

        Y <- cells(sp)
        X <- neighbourMean(sp, useWeights, normalise)

        Y <- Y[tumourCells, ]
        X <- X[tumourCells, ]

        list(X=X,Y=Y)
    })
    sizes <- sapply(XY, function(xy) nrow(xy$Y))

    all.x <- lapply(XY, function(xy) xy$X)
    all.y <- lapply(XY, function(xy) xy$Y)

    X <- do.call(rbind, all.x)
    Y <- do.call(rbind, all.y)
    return( list( X=X, Y=Y, sizes=sizes))
}
