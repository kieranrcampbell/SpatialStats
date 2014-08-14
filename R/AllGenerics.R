
## For SPData

#######################
## Set / Get methods ##
#######################


#' Cell measurement matrix
#'
#' Returns the n by p matrix of normalised cell measurements for
#' n cells and p channels.
#'
#' @param object The SPData object to use.
#'
#' @name cells
#' @rdname cells-methods
#' @exportMethod cells
setGeneric(name = "cells",
           def = function(object) standardGeneric("cells"))

#' Raw data matrix
#'
#' Returns the n by p matrix of raw measurements for n cells
#' and p channels.
#'
#' @param object The SPData object to use
#'
#' @name rawData
#' @rdname rawData-methods
#' @exportMethod rawData
setGeneric(name = "rawData",
           def = function(object) standardGeneric("rawData"))

#' Number of cells in the sample
#'
#' @param object The SPData object to use
#'
#' @name nCells
#' @rdname nCells-methods
#' @exportMethod nCells
setGeneric(name = "nCells",
           def = function(object) standardGeneric("nCells"))

#' Number of channels in the sample
#'
#' Returns the number of channels, which can be proteins and
#' protein modifications in the given sample.
#'
#' @param object The SPData object to use
#'
#' @name nChannel
#' @rdname nChannel-methods
#' @exportMethod nChannel
setGeneric(name = "nChannel",
           definition = function(object) standardGeneric("nChannel"))

#' Channel names used in the sample
#'
#' Returns the names of the channels used in the experiment. This
#' can be the (abbreviated) names of proteins and protein modifications,
#' or whatever else was measured in a given experiment.
#'
#' @param object The SPData object to use
#'
#' @name channels
#' @rdname channels-methods
#' @exportMethod channels
setGeneric(name = "channels",
           definition = function(object) standardGeneric("channels"))

#' List of nearest neighbour readouts
#'
#' Returns a list of length nCells(object).
#' The ith entry is an n by m matrix, for cell i having n neighbours
#' each of which have nChannel(object) channels.
#'
#' @param object The SPData object to use
#'
#' @name neighbours
#' @rdname neighbours-methods
#' @exportMethod neighbours
setGeneric(name = "neighbours",
           definition = function(object) standardGeneric("neighbours"))

#' Cell sizes
#'
#' Returns a vector of cell sizes,
#'
#' @param object The SPData object to use
#'
#' @name size
#' @rdname size-methods
#' @exportMethod size
setGeneric(name = "size",
           definition = function(object) standardGeneric("size"))

#' Boundary weights
#'
#' Returns a list of boundary weights. The \emph{i}th item
#' will be a vector of length \emph{n} if cell \emph{i} has
#' \emph{n} nearest neighbours.
#'
#' @name weight
#' @rdname weight-methods
#' @exportMethod weight
setGeneric("weight", function(object) standardGeneric("weight"))

#' Set the boundary weights
#'
#' Set the boundary weights.
#'
#' @param object The SPData object in which to set the weights
#' @param value A list of length nCells(object) of neighbour weights.
#'
#' @name weight<-
#' @rdname w-methods
#' @exportMethod weight<-
setGeneric("weight<-", function(object, value) standardGeneric("weight<-"))

setGeneric("cells<-", function(object, value) standardGeneric("cells<-"))

setGeneric("neighbours<-", function(object, value) standardGeneric("neighbours<-"))

setGeneric("ID", function(object) standardGeneric("ID"))

setGeneric("ID<-", function(object, value) standardGeneric("ID<-"))

setGeneric("neighbourIDs", function(object) standardGeneric("neighbourIDs"))



setGeneric("xy", function(object) standardGeneric("xy"))

setGeneric("xy<-", function(object, value) standardGeneric("xy<-"))

########################
## Plotting functions ##
########################

setGeneric("boxplots", function(object, nrow, ncol) standardGeneric("boxplots"))

setGeneric("channelPlot", function(object, channel.ids) standardGeneric("channelPlot"))

#####################################
## Boundary / cell class functions ##
#####################################

setGeneric("cellClass", function(object) standardGeneric("cellClass"))

setGeneric("cellClass<-", function(object, value) standardGeneric("cellClass<-"))

setGeneric("neighbourClass", function(object, cell.class) standardGeneric("neighbourClass"))

setGeneric("findBoundary", function(object) standardGeneric("findBoundary"))

###########################################
## Nearest neighbour averaging functions ##
###########################################

setGeneric("neighbourMean", function(object, useWeights, normalise) standardGeneric("neighbourMean"))

## For SPExp

setGeneric("getDir", function(object) standardGeneric("getDir"))

setGeneric("files", function(object) standardGeneric("files"))

setGeneric("loadExp", function(object) standardGeneric("loadExp"))

setGeneric("SPlist", function(object) standardGeneric("SPlist"))

setGeneric("IDs", function(object) standardGeneric("IDs"))
