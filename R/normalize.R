###################################################################
## Linear model based normalisation for signal against cell size ##
###################################################################


#' Normalization of an \code{SPData} object 
#' 
#' This method performs four steps of of normalization:
#' \enumerate{
#' \item Normalize by cell size
#' \item Normalize by cell concentration
#' \item Scale-centre the rows
#' \item Regenerate the nearest neighbour data from the scale centred data using the
#' nearest neighbour IDs
#' }
#' 
#' @details
#' Each step of normalization can be optionally turned on or off using the appropriate flag.
#' 
#' Setting \code{by.class} to true means each cell class is normalized separately within a sample,
#' while if set to false all cells are normalized together. Note if there are many cell classes
#' and few cells, normalization may become nonsensical - this is not checked for.
#' 
#' \code{norm.concentration} constructs a measure of cell concentration for each channel as the
#' sum across all other channels. It then regresses each channel against the measure of concentration
#' and accepts the residuals as the normalized quantity.
#' 
#' \code{norm.regenerateNN} reconstructs the neighbour data from the newly-normalized cell data.
#' The neighbour information is stored separately in an \code{SPData} object so each can be
#' manipulated / normalized without affecting the other.
#' 
#' @param sp The \code{SPData} instance to use
#' @param by.class Normalize each class of cells separately
#' @param norm.size Regress out a cell size dependency
#' @param norm.concentration Regress out the cell concentration dependency
#' @param norm.centrescale Scale each column of \code{cells(sp)} to have mean 0 and sd 1
#' @param norm.regenerateNN Regenerate nearest neighbour data after normalization
#' 
#' @return An object of class \code{SPData} that has been appropriately normalized.
#'
#' @export
NormalizeSP <- function(sp, by.class = TRUE, 
                        norm.size = TRUE,
                        norm.concentration = TRUE,
                        norm.centrescale = TRUE,
                        norm.regenerateNN = TRUE) {
  raw <- rawData(sp)
  
  if(norm.size) {
    Y <- lmNormalize(raw, size(sp), by.class, cellClass(sp))
  } else {
    Y <- raw
  }
  
  if(norm.concentration) Y <- TotalProteinNormalize(Y, by.class, cell.classes=cellClass(sp))
  if(norm.centrescale)  Y <- preprocess.centre(Y, by.class, cellClass(sp))
  
  if(norm.regenerateNN) {
    X <- lapply(neighbourIDs(sp), function(id) {
      Y[id,,drop=FALSE]
    })
  } else {
    X <- neighbours(sp)
  }
  sp@readouts <- Y
  sp@cellNeighbours <- X
  
  sp
}

# Centre scale a vector to mean 0 and sd 1
# 
# @param x The vector to centre-scale
# @return The centre-scaled vector
centreScale <- function(x) {
  (x - mean(x))/sd(x)
}

## converts the measurements for each channel
## to N(0,1)
## todo: come back and make for more than 2 classes
preprocess.centre <- function(Y, by.class, cell.class=NULL) {
  
  if(!by.class) {
    Y <- apply(Y, 2, centreScale)
  } else { # normalise each cell class separately
    if(is.null(cell.class)) stop("by.class==TRUE and cell.class==NULL")
    
    cell.indices <- lapply(unique(cell.class), function(cl) which(cell.class == cl))
    Y <- apply(Y, 2, function(y) {
      for(i in 1:length(cell.indices)) {
        y[ cell.indices[[i]] ] <- centreScale(y[ cell.indices[[i]] ])
      }
      return ( y )
    })
  }
  Y
}


#' Provides column by column normalisation on Y
#' given cell size s
#'
#' @param Y Measurement by channel matrix
#' @param s Vector of cell sizes
#' @param by.class Normalize each cell class separately
#' @param cell.classes A vector of cell classes as each type is normalized separately. Should be the same length
#' as \code{nrow(Y)}, and is most likely the output of \code{cellClass(sp)} for some \code{SPData} object \code{sp}
lmNormalize <- function(Y, s, by.class, cell.classes=NULL) {
  
  Y.n <- NULL
  if(!by.class) {
    cell.classes <- rep(1, nrow(Y))
  }
  y <- lapply(unique(cell.classes), function(i) Y[cell.classes == i,])
  sizes <- lapply(unique(cell.classes), function(i) s[cell.classes == i])
  
  y.norm <- lapply(unique(cell.classes), function(i) lmY(y[[i]], sizes[[i]]))
  
  Y.n <- matrix(NA, nrow=nrow(Y), ncol=ncol(Y))
  for(i in unique(cell.classes)) {    
    Y.n[cell.classes == i, ] <- y.norm[[i]]
  }
  colnames(Y.n) <- colnames(Y)
  
  if(any(is.na(Y.n))) stop("New matrix contains NAs - cell class fault")
  
  return(Y.n)
}

lmY <- function(Y, s) {
  Y <- apply(Y, 2, function(y) {
    fit <- lm(y ~ s)
    y - fit$fitted
  })
  Y
}


#' Uses lm normalisation to normalize by cell concentration
#'
#' @param Y A cell-by-channel matrix of measurements
#' @param by.class If true normalize each cell class separately
#' @param cell.classes The cell classes as passed by cellClass
#'
#' @return Concentration normalized cell-by-channel matrix
TotalProteinNormalize <- function(Y, by.class, cell.classes=NULL) {
  nchannels <- ncol(Y)
  Y.n <- NULL
  if(!by.class) {
    Y.n <- sapply(1:nchannels, function(i) {
      y <- Y[,i]
      totalP <- rowSums(Y[,-i])
      fit <- lm(y ~ totalP)
      y - fit$fitted
    })
  } else {
    y <- lapply(sort(unique(cell.classes)), function(i) Y[cell.classes == i,])
    y.concnorm <- lapply(y, function(yc) {
      new.yc <- sapply(1:nchannels, function(i) {
        V <- yc[,i]
        Z <- rowSums(yc[,-i])
        fit <- lm(V ~ Z)
        residuals(fit)
      })
      new.yc
    })
    Y.n <- matrix(NA, nrow=nrow(Y), ncol=ncol(Y))
    for(i in sort(unique(cell.classes))) {
      Y.n[cell.classes == i, ]  <- y.concnorm[[i]]
    }
    
  }
  
  colnames(Y.n) <- colnames(Y)
  return(Y.n)
}

