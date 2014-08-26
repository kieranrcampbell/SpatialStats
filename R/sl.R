
############################################################################
## Combines four sample regression then applies a linear model across all ##
## kieran.campbell@sjc.ox.ac.uk                                           ##
############################################################################

library(glmnet)
library(R.utils)
library(devtools)
library(lars)
library(covTest)

##load_all("../SpatialPRo/")

## find tumour ids with something like
# tid <- findHigherID(sp, "Keratin","mean")

dontrun <- function() {
  
  tumour.classes <- findIDs(SPE, "Vimentin", "lower", "median") # get the epitheilial classes in each sample
  ## bind data together and add factors
  XY <- BindSPE(SPE,choose.class = tumour.classes)
  X <- XY$X
  Y <- XY$Y
  
  factors <- ConstructSampleFactors(XY, IDs(SPE))
  X <- cbind(X, factors)
  
  ## randomise?
  ##X <- X[sample(nrow(X)),]
  
  npred <- ncol(X)
  include <- npred:(npred - ncol(factors) + 1)
  
  ## multisample split results
  multisplit.res <- GeneralLassoSig(y,X,B=50, s="usefixed",fixedP = 5, include=include)
  
  ## lars stuff
  lar <- apply(Y, 2, function(y) lars(X, y, "lasso", normalize=FALSE))
  cvtests <- lapply(1:length(lar), function(i) covTest(lar[[i]], X, Y[,i]))
  
  ## clean up the covariance test statistics and apply bonferroni
  cv.results <- AdjustCovtests(cvtests, ncol(X))
  
  ## multiple testing corrections
  nsamp <- 32
  alpha <- 0.05
  
  cv.results <- apply(cv.results, 1, p.adjust, method="BH")
  cv.results <- t(cv.results) # why??
  
  multisplit.res <- apply(multisplit.res, 1, p.adjust, method="BH")
  multisplit.res <- t(multisplit.res)
  
  
  ## results
  covtest.res <- which(cv.results < alpha, arr.ind=TRUE)
  multi.res <- which(multisplit.res < alpha, arr.ind=TRUE)
  rownames(multi.res) <- NULL
  
}



findOverlap <- function(a.results, remove=c("same", "different")) {
  sr <- NULL
  
  if(remove == "different") {
    sr <- lapply(a.results, function(mat) {
      cross <- which(mat[,1] != mat[,2])
      mat[cross,]
    })
  } else {
    sr <- lapply(a.results, function(mat) {
      cross <- which(mat[,1] == mat[,2])
      mat[cross,]
    })
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
  if(!is.matrix(pathways)) {
    pathways <- channels(SPE[[1]])[pathways]
  } else {
    pathways <- t(apply(pathways, 1, function(pw) channels(SPE[[1]])[pw]))
  }
  pathways
}

entire.project <- function() {
  set.seed(123)
  
  SPE@spdata <- lapply(SPlist(SPE), normaliseSP)
  
  immune.ind <- c(5,7,8,10,19,28)
  for(i in 1:length(SPE)) SPE[[i]] <-  SPE[[i]][,-immune.ind]
  
  withWeights <- doJointAnalysis(SPE, useWeights=TRUE)
  noWeights <- doJointAnalysis(SPE, useWeights=FALSE)
  
  pathways.with <- findOverlap(withWeights, remove="different")
  pathways.nowe <- findOverlap(noWeights, remove="different")
  
  pw.same <- findOverlap(noWeights, remove="same")
  
  print(pathways.with)
  print(pathways.nowe)
  
  ## statistical power calculation
  

}