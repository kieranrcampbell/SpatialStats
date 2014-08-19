
############################################################################
## Combines four sample regression then applies a linear model across all ##
## kieran.campbell@sjc.ox.ac.uk                                           ##
############################################################################

library(glmnet)
library(R.utils)
library(devtools)
library(lars)
library(covTest)

load_all("..")
source("hdimLasso.R")

#' Cleans up results of the covariance test statistic and applies the
#' bonferroni multiple testing correction for each response variable
AdjustCovtests <- function(cvtests, nvar) {
    sapply(cvtests, function(cv) {
        cv <- cv$results
        pred.num <- abs(cv[,1])

        ## if a predictor appears more than once, remove it from the
        ## table as it may be unreliable
        ta <- table(abs(pred.num))
        multiple <- as.numeric(names(which(ta > 1)))
        P <- cv[,3]

        if(length(multiple) > 0) {
            toRemove <- which(pred.num == multiple)
            pred.num <- pred.num[-toRemove]
            P <- P[-toRemove]
        }

        P <- p.adjust(P, method="bonferroni")
        retval <- rep(NA, nvar)
        retval[pred.num] <-  P
        retval[is.na(retval)] <- 1 # if a predictor doesn't appear at all, p-value of 1
        retval
    })
}



doJointAnalysis <- function(SPE, useWeights=TRUE) {
    ## bind data together and add factors
    XY <- BindSPE(SPE, useWeights)
    X <- XY$X
    Y <- XY$Y

    factors <- ConstructSampleFactors(XY, IDs(SPE))
    X <- cbind(X, factors)

    ## randomise?
    ##X <- X[sample(nrow(X)),]

    npred <- ncol(X)
    include <- npred:(npred - ncol(factors) + 1)

    ## multisample split results
    multisplit.res <- apply(Y, 2, function(y) {
        ps <- hdimLasso(y,X,B=50, s="usemin", minP=5, include=include)
        ps
    })


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

    return( list(covtest=covtest.res, multisplit=multi.res ))
}

## ## load data
## load("../data/SPE.rda")
## load("../data/SPE_bad.rda")

## ## construct SPE of good primary tumours
## SPE.primary <- SPExperiment("None", "None", spdata=c(SPlist(SPE)[4], SPlist(SPE.bad)),
##                         ids = c(IDs(SPE)[4], IDs(SPE.bad)))

## SPE <- SPE.primary

load("../data/SPE_report.rda")

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
        ## if(is.matrix(r.logic)) {
        ##     any(r.logic[,1] & r.logic[,2])
        ## } else {
        ##     r.logic[1] && r.logic[2]
        ## }
    })

    pathways <- sr[[1]][intersect,]
    if(!is.matrix(pathways)) {
        pathways <- channels(SPE[[1]])[pathways]
    } else {
        pathways <- t(apply(pathways, 1, function(pw) channels(SPE[[1]])[pw]))
    }
    pathways
}

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

stop("done")
fit.all <- lm(Y ~ X)
s <- summary(fit)
rs <- sapply(s, function(x) x$adj.r.squared)

mrs <-  mean(rs)
nrow(X)

sp <- SPE[[1]]
y <- cells(sp)
x <- neighbourMean(sp, F, T)

tumourID <- findTumourID(sp)
tumourCells <- which(cellClass(sp) == tumourID)
y <- y[tumourCells, ]
x <- x[tumourCells, ]
fit.1 <- lm(y ~ x)
s1 <- summary(fit.1)
nrow(x)

rs1 <- sapply(s1, function(x) x$adj.r.squared)
mrs1 <-  mean(rs1)
