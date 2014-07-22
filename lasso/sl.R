
############################################################################
## Combines four sample regression then applies a linear model across all ##
## kieran.campbell@sjc.ox.ac.uk                                           ##
############################################################################

##library(devtools)
library(glmnet)
library(R.utils)
library(devtools)

load_all("..")

##load("../data/SPE.Rd")
## load("../data/SPE_bad.Rd")

## primary.er1 <- c(2, 4)
## primary.er2 <- c(3)

## sp.list <- c(SPlist(SPE)[primary.er1], SPlist(SPE.bad)[primary.er2])
## ids <- sapply(sp.list, id)
## SPE.primaryl <- SPExperiment(dir="none", files="none", spdata=sp.list, ids=ids)


source("hdimLasso.R")

load("../data/SPE_prim.Rd")
load("../data/SPE.Rd")
load("../data/SPE_bad.Rd")
##SPE <- SPE.bad

SPE.all <- SPExperiment("None", "None", spdata=c(SPlist(SPE), SPlist(SPE.bad)),
                        ids = c(IDs(SPE), IDs(SPE.bad)))

SPE <- SPE.all

## variable selection
XY <- lapply(SPlist(SPE), function(sp) {
    tumourID <- findTumourID(sp)
    tumourCells <- which(cellClass(sp) == tumourID)

    Y <- cells(sp)
    X <- neighbourMean(sp, FALSE, TRUE)

    Y <- Y[tumourCells, ]
    X <- X[tumourCells, ]
    ##print(colMeans(Y))

    list(X=X,Y=Y)
})

allRegress <- function(M) {
    Y.new <- matrix(0, nrow=nrow(M), ncol=ncol(M))
    for(i in 1:32) {
        fit <- glmnet(M[,-i], M[,i])
        cvf <- cv.glmnet(M[,-i], M[,i])
        r <- M[,i] - predict(fit, newx=M[,-i],s=cvf$lambda.1se)
        Y.new[,i] <- r
    }
    colnames(Y.new) <- colnames(M) ; rownames(Y.new) <- rownames(M)
    return( Y.new )
}


all.x <- lapply(XY, function(xy) xy$X)
all.y <- lapply(XY, function(xy) xy$Y)

X <- do.call(rbind, all.x)
Y <- do.call(rbind, all.y)

##Y.new <- allRegress(Y)

cell.sizes <- sapply(XY, function(xy) dim(xy$Y)[1])
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
         main=channels(SPE.primaryl[[1]])[i], cex.main=1.5, xvar="dev")

}
dev.off()


pdf("img/xy_cv.pdf",width=22,height=15)
par(mfrow=c(4,8))
for(i in 1:32) {
    plot(cv.glmnet(X, Y[,i], standardize=FALSE), main=channels(SPE.primaryl[[1]])[i], cex.main=1.5)
}
dev.off()

ind <- c(2, 5, 8, 10, 11, 13, 15, 16, 17, 18, 22, 24,30)


## let's do them separately

## for(i in 1:3) {
##     y <- all.y[[i]]
##     x <- all.x[[i]]

##     pdf(paste("img/cv",i,".pdf",sep=""), width=22, height=15)
##     par(mfrow=c(4,8))
##     for(j in 1:32) {
##         cvf <- cv.glmnet(x,y[,j])
##         plot(cvf, main=channels(sp)[j])
##     }
##     dev.off()
## }


## y <- all.y[[1]]
## library(bnlearn)
## df <- data.frame(y)
## dag <- hc(df, score='bic-g')
## fit <- bn.fit(dag, df)



## lar <- lars(X, Y[,k], "lasso")
## covTest(lar,X, Y[,k])

all.ps <- apply(Y, 2, function(y) {
    ps <- hdimLasso(y,X,B=50, s="lambda.min")
    as.numeric(ps < 0.01)
})
which(all.ps == 1, arr.ind=TRUE)


## lars stuff
lar <- apply(Y, 2, function(y) lars(X, y, "lasso", normalize=FALSE))
cvtests <- lapply(1:length(lar), function(i) covTest(lar[[i]], X, Y[,i]))

sig <- sapply(cvtests, function(cv) {
    cv <- cv$results
    pred.num <- abs(cv[,1])
    P <- cv[,3]
    ##P <- p.adjust(P, method="BH")
    retval <- rep(0, length(pred.num))
    retval[pred.num] <-  as.numeric(P < 0.01)
    retval
})

which(sig == 1, arr.ind=TRUE)
