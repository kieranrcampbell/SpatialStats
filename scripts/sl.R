
############################################################################
## Combines four sample regression then applies a linear model across all ##
## kieran.campbell@sjc.ox.ac.uk                                           ##
############################################################################


library(lars)
library(covTest)

library(SpatialPRo)
library(SpatialStats)


# Main analysis begins here -----------------------------------------------


set.seed(123)

SPE@spdata <- lapply(SPlist(SPE), normalizeSP)

#SPE <- SPE[1:4]

spe.list <- list()
immune.ind <- c(5,7,8,10,19,28)
for(i in 1:length(SPE)) spe.list[[i]] <-  SPE[[i]][,-immune.ind]

spe.noimmune <- SPExperiment(getDir(SPE), files(SPE), spe.list)

tumour.classes <- findIDs(spe.noimmune, "Vimentin", "lower", "median") # get the epitheilial classes in each sample
## bind data together and add factors

XY <- BindSPE(spe.noimmune, choose.class = tumour.classes)
X <- XY$X
Y <- XY$Y

factors <- ConstructSampleFactors(XY, IDs(SPE))
X <- cbind(X, factors)

## randomise?
##X <- X[sample(nrow(X)),]

npred <- ncol(X)
include <- npred:(npred - ncol(factors) + 1)

## multisample split results
multisplit.res <- GeneralLassoSig(Y,X,B=100, s="usefixed",fixedP = 5, include=include)

## lars stuff
lar <- apply(Y, 2, function(y) lars(X, y, "lasso", normalize=FALSE))
cvtests <- lapply(1:length(lar), function(i) covTest(lar[[i]], X, Y[,i]))

## clean up the covariance test statistics and apply bonferroni
cv.results <- AdjustCovtests(cvtests, ncol(X))

## multiple testing corrections
nsamp <- 26
alpha <- 0.05

cv.results <- apply(cv.results, 1, p.adjust, method="BH")
cv.results <- t(cv.results) # why??

multisplit.res <- apply(multisplit.res, 1, p.adjust, method="BH")
multisplit.res <- t(multisplit.res)


## results
covtest.res <- which(cv.results < alpha, arr.ind=TRUE)
multi.res <- which(multisplit.res < alpha, arr.ind=TRUE)
rownames(multi.res) <- NULL

all.res <- list(covtest.res, multi.res)

pathways.nowe <- findOverlap(all.res, remove="different")

pw.same <- findOverlap(all.res, remove="same")

print(pathways.nowe)
print(pw.same)


writeClass <- function(sp) {
  cl <- cellClass(sp)
  fname <- paste(ID(sp), ".mat", sep="")
  writeMat(cellclass=cl, con=fname)
}