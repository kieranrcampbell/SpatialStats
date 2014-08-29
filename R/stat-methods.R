
#' Bind multiple sample data together
#'
#' @description
#' The low power and high noise in some analyses makes it necessary to bind multiple
#' samples together. Given an SPExp object, this method separates out each cell type,
#' calculate neighbour means and binds the samples together, i.e. the individual response
#' matrices for different samples \eqn{Y1}, \eqn{Y2}, \eqn{Y3} are row-bound together to
#' form one large response matrix \eqn{Y = [Y1' Y2' Y3']'}.
#' 
#' @details
#' If \code{normalize = TRUE} it is necessary to introduce sample factors as constructed
#' by \code{\link{ConstructSampleFactors}}, otherwise the regression estimates will be affected
#' by the relative means of different samples.
#'
#' @param SPE The SPExp object to use
#' @param choose.class Choose a particular class of cells. Should be a vector of length 
#' \code{length{SPE}} where the ith entry indicates the class of cells in sample i to pick.
#' @param use.weights Passed to neighbourMeans to optionally weight the regression by relative boundary size
#' @param normalize If true predictor columns are centre-scaled
#' 
#' @return A list with three components:
#' \itemize{
#' \item{\code{X} }{The bound predictor matrix}
#' \item{\code{Y} }{The bound response matrix}
#' \item{\code{sizes} }{A vector with the number of cells selected from each \code{SPE}}
#' }
#'
#' @export
BindSPE <- function(SPE, choose.class=NULL, use.weights=FALSE, normalize=TRUE) {
  if(!is.null(choose.class) && (length(choose.class) != length(SPE))) {
    stop("If choose.class is not null a class must be specified for every sample")
  }  
  
  ## variable selection
  XY <- lapply(1:length(SPE), function(i) {
    sp <- SPE[[i]]
    Y <- cells(sp)
    X <- neighbourMean(sp, use.weights, normalize)
    
    if(!is.null(choose.class)) {
      classID <- choose.class[i]
      chosenCells <- which(cellClass(sp) == classID)
      Y <- Y[chosenCells, ]
      X <- X[chosenCells, ]
    }
    
    list(X=X,Y=Y)
  })
  sizes <- sapply(XY, function(xy) nrow(xy$Y))
  
  all.x <- lapply(XY, function(xy) xy$X)
  all.y <- lapply(XY, function(xy) xy$Y)
  
  X <- do.call(rbind, all.x)
  Y <- do.call(rbind, all.y)
  
  ## normalize Y & X to have mean 0 and unit variance
  if(normalize) {
    m0uv <- function(x) (x - mean(x)) / sd(x)
    Y <- apply(Y, 2, m0uv)
    X <- apply(X, 2, m0uv)
  }
  
  return( list( X=X, Y=Y, sizes=sizes))
}

#' Construct model matrix for integrating multiple samples
#'
#' Given the output of the function BindSPE we would like to control
#' for biases introduced by binning all the data together. ConstructSampleFactors
#' takes the output of BindSPE and returns a model matrix of factors highlighting
#' what sample a given cell is from
#'
#' @param XY An \code{X-Y-sizes} list - the output of \code{BindSPE}
#' @param sample.ids The IDs of the samples that have been bound together - for naming the factors
#' @param intercept If true an intercept is modelled when calling the function \code{model.matrix}
#' but then removed from the resulting matrix so it isn't modelled twice when \code{lm} is
#' later called. If false, the call to \code{model.matrix} is performed without an intercept, and
#' the resulting model.matrix contains all sample factors - this matrix will have one more column
#' that if \code{intercept = FALSE}.
#' 
#' @return A model matrix that correctly accounts for differing
#' means of samples as the result of binding them together
#' using the function \code{BindSPE}. 
#' 
#' @export
ConstructSampleFactors <- function(XY, sample.ids, intercept=TRUE) {
  cell.sizes <- XY$sizes
  Nfactors <- length(cell.sizes) - 1
  factors <- NULL
  #   for(i in 1:Nfactors) {
  #     tcol <- rep(0, sum(cell.sizes))
  #     cs <- cumsum(cell.sizes)
  #     range <- (cs[i] + 1):cs[i+1]
  #     tcol[ range ] <- 1
  #     factors <- cbind(factors, tcol)
  #   }
  sample.factors <- sapply(1:length(cell.sizes), function(i) rep(sample.ids[i], cell.sizes[i]))
  sample.factors <- as.factor(unlist(sample.factors))
  
  m <- NULL
  if(intercept) {
    m <- model.matrix( ~ sample.factors)
    m <- m[,-1] # remove intercept since it's already modelled  
  } else {
    m <- model.matrix( ~ 0 + sample.factors)
  }  
  samples.used <- setdiff(1:length(cell.sizes), intercept*1)
  
  ids.used <- strsplit(colnames(m), "factors")
  ids.used <- sapply(ids.used, function(x) x[2])
  colnames(m) <- ids.used
  return( m )
}

#' Clean up and bonferroni adjust multiple \code{covTest} results
#' 
#' @param cvtests A list of results from \code{covTest}. More than one since
#' we fit a different model for each response variable.
#' @param nvar The original number of predictor variables - this info isn't
#' maintained in the result.
#' 
#' @export
#' 
#' @return An n-by-m matrix of p-values for m predictor and n response variables.
AdjustCovtests <- function(cvtests, nvar) {
  adjusted <- sapply(cvtests, function(cv) {
    cv <- cv$results
    pred.num <- abs(cv[,1])
    
    ## if a predictor appears more than once, remove it from the
    ## table as it may be unreliable
    ta <- table(abs(pred.num))
    multiple <- as.numeric(names(which(ta > 1)))
    P <- cv[,3]
    
    if(length(multiple) > 0) {
      toRemove <- which(pred.num %in% multiple)
     
      pred.num <- pred.num[-toRemove]
      P <- P[-toRemove]
    }
    
    P <- p.adjust(P, method="bonferroni")
    retval <- rep(NA, nvar)
    retval[pred.num] <-  P
    retval[is.na(retval)] <- 1 # if a predictor doesn't appear at all, p-value of 1
    retval
  })
  return( adjusted )
}

#' Given a list of two pathway results, find the overlap between them
#' 
#' The result of the spatial pathway identification is a two-column matrix, where each row represents the numeric
#' identifiers of the causal direction of the pathway, so
#' 1 4
#' 3 2
#' would indicate there is a pathway from component 1 to 4, and from 3 to 2. If we use multiple methods to find these,
#' it is useful to find the overlap - pathways reported from both methods. This method does that, effectively taking the
#' intersection of the two matrices and returning the results in a similar matrix.
#' 
#' @param results A list containg two 'pathway' matrices
#' @param remove Which results to keep: if \code{"same"} then 1 1 would be reported but 1 2 discarded, but if
#' \code{"different"} then 1 2 would be reported but 1 1 discarded (useful for finding pathways that don't 
#' go from same component to same component)
#' 
#' @return A matrix showing common pathways.
#' @export
FindOverlap <- function(results, remove=c("same", "different")) {
  sr <- NULL
  
  if(remove == "different") {
    sr <- lapply(results, function(mat) {
      cross <- which(mat[,1] != mat[,2])
      mat[cross,,drop=FALSE]
    })
  } else {
    sr <- lapply(results, function(mat) {
      cross <- which(mat[,1] == mat[,2])
      mat[cross,,drop=FALSE]
    })
  }
  
  if(0 %in% dim(sr[[1]] || 0 %in% dim(sr[[2]]))) {
    ## one is empty
    return( NULL )
  }
  
  intersect <- apply(sr[[1]], 1, function(row) {
    if(is.matrix(sr[[2]])) {
      r.logic <- apply(sr[[2]], 1, function(r) {
        l <- r == row
        l[1] & l[2]
      })
      any(r.logic)
    } else {
      l <- sr[[2]] == row
      l[1] & l[2]
    }
    
  })
  
  pathways <- sr[[1]][intersect,]
#   if(!is.matrix(pathways)) {
#     pathways <- channels(SPE[[1]])[pathways]
#   } else {
#     pathways <- t(apply(pathways, 1, function(pw) channels(SPE[[1]])[pw]))
#   }
  pathways
}
