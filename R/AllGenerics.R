
## For SPData
setGeneric("cells", function(object) standardGeneric("cells"))

setGeneric("nCells", function(object) standardGeneric("nCells"))

setGeneric("nChannel", function(object) standardGeneric("nChannel"))

setGeneric("channels", function(object) standardGeneric("channels"))

setGeneric("neighbours", function(object) standardGeneric("neighbours"))

setGeneric("nnReg", function(object) standardGeneric("nnReg"))

setGeneric("size", function(object) standardGeneric("size"))

setGeneric("cells<-", function(object, value) standardGeneric("cells<-"))

setGeneric("id", function(object) standardGeneric("id"))

setGeneric("id<-", function(object, value) standardGeneric("id<-"))

setGeneric("neighbourIDs", function(object) standardGeneric("neighbourIDs"))

setGeneric("weights", function(object) standardGeneric("weights"))

setGeneric("boxplots", function(object) standardGeneric("boxplots"))

## For SPExp

setGeneric("getDir", function(object) standardGeneric("getDir"))

setGeneric("files", function(object) standardGeneric("files"))

setGeneric("loadExp", function(object) standardGeneric("loadExp"))

setGeneric("data", function(object) standardGeneric("data"))

setGeneric("ids", function(object) standardGeneric("ids"))
