
############################################################################
## Combines four sample regression then applies a linear model across all ##
## kieranrcampbell@gmail.com                                              ##
############################################################################

doReg <- function() {
library(devtools)
##load_all("..")

load("../data/SPE.Rd")

## phInd for phosphates
## repInd for response variables
load("../../usefulin.Rd")

nsamp <- 4

phNames <- channels(SPE[[1]])[phInd]

## variable selection
XY <- lapply(1:nsamp, function(i) {
    sp <- SPE[[i]]
    tumourID <- findTumourID(sp)
    tumourCells <- which(cellClass(sp) == tumourID)

    Y <- cells(sp)
    X <- neighbourMean(sp, FALSE, TRUE)

    Y <- Y[tumourCells, repInd]
    X <- X[tumourCells, phInd]

    list(X=X,Y=Y)
})

## make x & y into the shape we want
all.x <- lapply(XY, function(xy) xy$X)
all.y <- lapply(XY, function(xy) xy$Y)

x <- do.call(rbind, all.x)
y <- do.call(rbind, all.y)

## remove colinear predictors in x
x <- VIFRemove(x, vif.threshold=2)


cell.sizes <- sapply(XY, function(xy) dim(xy$Y)[1])

ids <- ids(SPE[1:nsamp])
sample.factors <- lapply(1:nsamp, function(i) rep(ids[i], cell.sizes[i]))
sample.factors <- as.factor(unlist(sample.factors))

fit <- lm(y ~ x + sample.factors)

## time for some significance testing
alpha <- 0.05
m <- ncol(x) * ncol(y)
print(sprintf("Adjusted p-value: %f", alpha/m))

s <- summary(fit)

nrep <- ncol(y) ; npred <- ncol(x)

A <- matrix(0,nrow=npred,ncol=nrep)

for(i in 1:nrep) {
    ## looking at which proteins influence i
    coeff <- s[[i]]$coefficients
    p.vals <- as.numeric(coeff[,4])
    p.vals <- p.vals[-c(1,npred+1, npred+2, npred+3)] # remove intercept consideration & sample effects
    A[,i] <- as.numeric(p.vals < alpha/m)
}


rownames(A) <- colnames(x) ; colnames(A) <- colnames(y)

z <- nrep+npred
a <- matrix(0, nrow=z,ncol=z)
colnames(a) <- rownames(a) <- c(colnames(x),colnames(y))
a[1:npred, (npred+1):z] <- A

g <- graph.adjacency(a, mode="directed")
plot(g)

}
