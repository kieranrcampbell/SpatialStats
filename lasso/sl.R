
############################################################################
## Combines four sample regression then applies a linear model across all ##
## kieran.campbell@sjc.ox.ac.uk                                           ##
############################################################################

##library(devtools)
library(glmnet)
library(R.utils)
library(devtools)
library(lars)
library(covTest)

load_all("..")

##load("../data/SPE.Rd")
## load("../data/SPE_bad.Rd")

## primary.er1 <- c(2, 4)
## primary.er2 <- c(3)

## sp.list <- c(SPlist(SPE)[primary.er1], SPlist(SPE.bad)[primary.er2])
## ids <- sapply(sp.list, id)
## SPE.primaryl <- SPExperiment(dir="none", files="none", spdata=sp.list, ids=ids)


source("hdimLasso.R")

load("../data/SPE.Rd")
load("../data/SPE_bad.Rd")
##SPE <- SPE.bad

SPE.primary <- SPExperiment("None", "None", spdata=c(SPlist(SPE)[4], SPlist(SPE.bad)),
                        ids = c(IDs(SPE)[4], IDs(SPE.bad)))

SPE <- SPE.primary

bindSPE <- function(SPE, pickTumour=TRUE) {
    ## variable selection
    XY <- lapply(SPlist(SPE), function(sp) {
        if(pickTumour) {
            tumourID <- findTumourID(sp)
            tumourCells <- which(cellClass(sp) == tumourID)
        }

        Y <- cells(sp)
        X <- neighbourMean(sp, TRUE, TRUE)

        Y <- Y[tumourCells, ]
        X <- X[tumourCells, ]

        list(X=X,Y=Y)
    })
    sizes <- sapply(XY, function(xy) nrow(xy$Y))

    all.x <- lapply(XY, function(xy) xy$X)
    all.y <- lapply(XY, function(xy) xy$Y)

    X <- do.call(rbind, all.x)
    Y <- do.call(rbind, all.y)
    return( list( X=X, Y=Y, sizes=sizes))
}

XY <- bindSPE(SPE)
X <- XY$X
Y <- XY$Y

cell.sizes <- XY$sizes
Nfactors <- length(cell.sizes) - 1
factors <- NULL
for(i in 1:Nfactors) {
    tcol <- rep(0, sum(cell.sizes))
    cs <- cumsum(cell.sizes)
    range <- (cs[i] + 1):cs[i+1]
    tcol[ range ] <- 1
    factors <- cbind(factors, tcol)
}

## normalise Y & X to have mean 0 and unit variance
m0uv <- function(x) (x - mean(x)) / sd(x)
Y <- apply(Y, 2, m0uv)
X <- apply(X, 2, m0uv)


colnames(factors) <- paste("sample", IDs(SPE)[1:Nfactors], sep="")

X <- cbind(X, factors)


pdf("img/xy.pdf",width=22,height=15)
par(mfrow=c(4,8))
for(i in 1:32) {
    plot(glmnet(X, Y[,i], standardize=FALSE), label="true",
         main=channels(SPE.primary[[1]])[i], cex.main=1.5, xvar="dev")

}
dev.off()


pdf("img/xy_cv.pdf",width=22,height=15)
par(mfrow=c(4,8))
for(i in 1:32) {
    plot(cv.glmnet(X, Y[,i], standardize=FALSE), main=channels(SPE.primary[[1]])[i], cex.main=1.5)
}
dev.off()

ind <- c(2, 5, 8, 10, 11, 13, 15, 16, 17, 18, 22, 24,30)



all.ps <- apply(Y, 2, function(y) {
    ps <- hdimLasso(y,X,B=50, s="usemin", minP=5)
    as.numeric(ps < 0.05)
})

which(all.ps == 1, arr.ind=TRUE)


## lars stuff
lar <- apply(Y, 2, function(y) lars(X, y, "lasso", normalize=FALSE))
cvtests <- lapply(1:length(lar), function(i) covTest(lar[[i]], X, Y[,i]))

sig <- sapply(cvtests, function(cv) {
    cv <- cv$results
    pred.num <- abs(cv[,1])
    P <- cv[,3]
    P <- p.adjust(P, method="bonferroni")
    retval <- rep(0, length(pred.num))
    retval[pred.num] <-  as.numeric(P < 0.05)
    retval
})

which(sig == 1, arr.ind=TRUE)

cn <- channels(SPE[[1]])
stop("done")
buildGraphfromSigMat <- function(mat, cnames) {
    require(igraph)
    diffRows <- which(mat[,1] != mat[,2])
    mat <- mat[diffRows,]
    from <- cnames[ mat[,1] ]
    to <- cnames[ mat[,2] ]
    from <- paste("NN", from,sep="-")

    edgelist <- rep(NA, 2*length(to))
    type <- rep(0:1, length(to))
    edgelist[c(T,F)] <- from
    edgelist[c(F,T)] <- to
    g <- graph.data.frame(data.frame(from=from,to=to))

    g
}


mat1 <- which(all.ps == 1, arr.ind=TRUE)
mat2 <- which(sig == 1, arr.ind=TRUE)

commonRows <- apply(mat1, 1, function(row) {
    l1 <- mat2 == row
    any(l1[,1] & l1[,2])
})

commonRows <- unlist(commonRows)

interactions <- mat1[commonRows,]

g <- buildGraphfromSigMat(interactions, channels(SPE[[1]]))
