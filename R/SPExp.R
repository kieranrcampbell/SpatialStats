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
    cat( " ", "With ", length(files(object)), " samples \n", sep="")
    cat(" ", files(object), "\n", sep=" \n")
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
          function(.Object, edir) {
              .Object@dir <- edir
              .Object@files <- dir(edir)
              .Object@spdata <- list()

              s <- strsplit(files(.Object),"_ID")
              .Object@ids <- as.numeric(sapply(s, function(x) strsplit(x[2],"_")[[1]][1]))

              return(.Object)
          })

#' Returns the sample ids.
#'
#' @export
setMethod("ids", "SPExp",
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
                                                          sep=""), id=ids(object)[i])
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


#' @export
setMethod("data", "SPExp",
          function(object) object@spdata)

#' @export
SPExperiment <- function(directory) {
    return ( new("SPExp", edir=directory) )
}
