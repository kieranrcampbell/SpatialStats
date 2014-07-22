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
