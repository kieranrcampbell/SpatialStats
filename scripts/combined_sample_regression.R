
############################################################################
## Combines four sample regression then applies a linear model across all ##
## kieran.campbell@sjc.ox.ac.uk                                           ##
############################################################################

library(devtools)
##load_all("..")

load("../data/SPE_bad.Rd")
##SPE <- SPE[c(1,3,4)]
SPE <- SPE_bad

## phInd for phosphates
## repInd for response variables
load("../../usefulin.Rd")

nsamp <- 4

phNames <- channels(SPE[[1]])[phInd]

ampkInd <- c(13,14)

## variable selection
XY <- lapply(1:nsamp, function(i) {
    sp <- SPE[[i]]
    tumourID <- findTumourID(sp)
    tumourCells <- which(cellClass(sp) == tumourID)

    Y <- cells(sp)
    X <- neighbourMean(sp, FALSE, TRUE)

    Y <- Y[tumourCells, ]
    X <- X[tumourCells,]

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

aids <- IDs(SPE[1:nsamp])
sample.factors <- lapply(1:nsamp, function(i) rep(aids[i], cell.sizes[i]))
sample.factors <- as.factor(unlist(sample.factors))

## standardize across variables
standardize <- function(x) (x - mean(x)) / sd(x)
y <- apply(y, 2, standardize)
x <- apply(x, 2, standardize)

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
    sample.coeff <- (npred+2):(npred + nsamp)
    p.vals <- p.vals[-c(1,sample.coeff)] # remove intercept consideration & sample effects
    A[,i] <- as.numeric(p.vals < alpha/m)
}


rownames(A) <- colnames(x) ; colnames(A) <- colnames(y)

z <- nrep+npred
a <- matrix(0, nrow=z,ncol=z)
colnames(a) <- rownames(a) <- c(colnames(x),colnames(y))
a[1:npred, (npred+1):z] <- A

g <- graph.adjacency(a, mode="directed")
plot(g)

