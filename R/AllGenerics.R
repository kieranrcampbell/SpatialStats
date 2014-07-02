
## For SPData

#######################
## Set / Get methods ##
#######################

setGeneric("cells", function(object) standardGeneric("cells"))

setGeneric("nCells", function(object) standardGeneric("nCells"))

setGeneric("nChannel", function(object) standardGeneric("nChannel"))

setGeneric("channels", function(object) standardGeneric("channels"))

setGeneric("neighbours", function(object) standardGeneric("neighbours"))

setGeneric("nnReg", function(object) standardGeneric("nnReg"))

setGeneric("size", function(object) standardGeneric("size"))

setGeneric("cells<-", function(object, value) standardGeneric("cells<-"))

setGeneric("neighbours<-", function(object, value) standardGeneric("neighbours<-"))

setGeneric("id", function(object) standardGeneric("id"))

setGeneric("id<-", function(object, value) standardGeneric("id<-"))

setGeneric("neighbourIDs", function(object) standardGeneric("neighbourIDs"))

setGeneric("weight", function(object) standardGeneric("weight"))

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

## For SPExp

setGeneric("getDir", function(object) standardGeneric("getDir"))

setGeneric("files", function(object) standardGeneric("files"))

setGeneric("loadExp", function(object) standardGeneric("loadExp"))

setGeneric("SPlist", function(object) standardGeneric("SPlist"))

setGeneric("ids", function(object) standardGeneric("ids"))
