
## For SPData
setGeneric("cells", function(object) standardGeneric("cells"))

setGeneric("nCells", function(object) standardGeneric("nCells"))

setGeneric("nProt", function(object) standardGeneric("nProt"))

setGeneric("pNames", function(object) standardGeneric("pNames"))

setGeneric("NN", function(object) standardGeneric("NN"))

setGeneric("nnReg", function(object) standardGeneric("nnReg"))

setGeneric("size", function(object) standardGeneric("size"))

setGeneric("NN<-", function(object, value) standardGeneric("NN<-"))

setGeneric("cells<-", function(object, value) standardGeneric("cells<-"))

setGeneric("id", function(object) standardGeneric("id"))

setGeneric("id<-", function(object, value) standardGeneric("id<-"))

setGeneric("nnID", function(object) standardGeneric("nnID"))

## For SPExp

setGeneric("getDir", function(object) standardGeneric("getDir"))

setGeneric("files", function(object) standardGeneric("files"))

setGeneric("loadExp", function(object) standardGeneric("loadExp"))

setGeneric("data", function(object) standardGeneric("data"))

setGeneric("ids", function(object) standardGeneric("ids"))
