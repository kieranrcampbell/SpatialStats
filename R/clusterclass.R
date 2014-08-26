########################################
## Tools for clustering by cell class ##
########################################


#' Clusters cells into different classes using PCA-Gaussian Mixture model clustering
#' 
#' A method of using high-dimensional proteomics data to classify cellular phenotypes (e.g. 
#' epithelial and stromal cells may be present within a tumour biopsy). First principle component
#' analysis is performed on the raw cell-by-channel matrix \code{Y}. Then the top \code{n.pc}
#' (default 3) principle components are used in Gaussian mixture modelling through 
#' expectation maximisation, which clusters the cellular readouts into \code{nclass} (default 2)
#' different classes. The underlying assumption is that in the dimensionality-reduced PC-space
#' the cells will be phenotypically distinct enough to appear drawn from two separate 
#' classes. 
#' 
#' This is an unsupervised method - other methods may be used (e.g. training a SVM
#' on known cell phenotypes). The cell classes can always be set by alternative methods using
#' \code{cellClass(sp) <- classes}
#' for an \code{SPData} \code{sp} and a vector of cell classes \code{classes}.
#'
#' @param Y A cell by channel matrix of readouts (e.g. from \code{cells(sp)})
#' @param doPCA If true, first performs PCA and uses the top
#' principle components for clustering
#' @param nclass The number of classes into which to cluster the cells (default = 2)
#' @param n.pc The number of principle components to use in Gaussian mixture modelling
#' (default = 3)
#' 
#' @references
#' Chen, W.-C., Maitra, R., Melnykov, V. (2012) EMCluster: EM Algorithm for Model-Based
#' Clustering of Finite Mixture Gaussian Distribution. R Package, URL
#' http://cran.r-project.org/package=EMCluster
#' 
#' @return A vector of length number of cells with numeric values corresponding
#' to distinct classes for each cell.
#'
#' @export
clusterClass <- function(Y, doPCA = TRUE, nclass=2, n.pc =3) {
  d <- NULL
  if(doPCA) {
    ypca <- princomp(Y)
    d <- ypca$scores[,1:n.pc]
  } else {
    d <- Y
  }
  
  emobj <- simple.init(d, nclass=nclass)
  emobj <- shortemcluster(d, emobj)
  ret <- emcluster(d, emobj, assign.class = TRUE)
  ret$class
}

#' Find the cell class with the highest/lowest average channel
#'
#' Given tissue samples with different classes, it may be advantageous to single out a 
#' class that corresponds to a certain phenotype. These are easiest found by considering
#' the cell class with the highest or lowest average (mean/median) expression of a given tissue,
#' e.g. in epithelial cells keratin is over-expressed while in stromal cells vimentin
#' is under-expressed.
#' 
#' @details
#' Most normalization routines normalize each class separately. However, this will make each
#' channel \eqn{N(0,1)} within a cell class, so comparison is just noise. This method normalizes 
#' the selected channel w.r.t. cell size and concentration, then picks out the highest or
#' lowest average (mean/median, set by \code{method}).
#' 
#' \code{channel} greps the channel name to find an index containing the string, but make sure that
#' string exists in a channel and only once, otherwise an error will be thrown.
#' 
#' \code{direction} If \code{"higher"} returns the class with the highest mean/median of the selected
#' channel, while \code{"lower"} will return the class with the lowest mean/median.
#'
#' @param sp An \code{SPData} object
#' @param channel A (partial) channel name 
#' @param direction Whether to look for the highest or lowest average
#' @param method The method for average (mean or median)
#' 
#' @return The class index with the higher average in the channel
#' 
#' @rdname findid-methods
#' @examples
#' \dontrun{
#' ## Keratin is higher in epithelial cells:
#' keratin.class <- findID(sp, "Keratin", direction="higher", method="median")
#' }
#' 
#' @export
findID <- function(sp, channel, direction=c("higher","lower"), method=c("mean","median")) {
  Y <- rawData(sp)
  
  index <- grep(channel, channels(sp))
  if(length(index) == 0) stop(channel + " not found in channels!")
  if(length(index) > 1) stop("More than channel found corresponding to " + channel)
  
  Yk <- Y[,index]
  y.excl <- rowSums(Y[,-index])
  
  norm.fit <- lm(Yk ~ size(sp) + y.excl)
  
  Yk <- residuals(norm.fit)
   
  classes <- lapply(unique(cellClass(sp)), function(i) which(cellClass(sp) == i))
  
  class.means <- sapply(classes, function(cell.list) do.call(method, list((Yk[cell.list]))))
  
  classind <- NULL
  if(direction == "higher") {
    classind <- which.max(class.means)
  } else if(direction == "lower")  {
    classind <- which.min(class.means)
  } 
  return(unique(cellClass(sp))[classind])
}

#' Apply \code{\link{findID}} across an entire \code{\link{SPExp}}
#' 
#' @param SPE An \code{SPExp} object
#' 
#' @return A vector of cell classes corresponding to each sample.
#' 
#' @rdname findid-methods
#' @export
findIDs <- function(SPE, channel, direction=c("higher","lower"), method=c("mean","median")) {
  sapply(SPlist(SPE), findID, channel, direction, method)
}
