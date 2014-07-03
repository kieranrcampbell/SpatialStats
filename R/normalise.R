
############################################################
## LOESS based normalisation for signal against cell size ##
## kieranrcampbell@gmail.com                              ##
############################################################

#' Provides column by column normalisation on Y
#' given cell size s
#'
#' @param Y Measurement by channel matrix
#' @param s Vector of cell sizes
loessNormalise <- function(Y, s) {

    Y <- apply(Y, 2, function(y) {
        fit <- loess(y ~ s)
        ##j <- order(s)
        ##plot(fit)
        ##lines(s[j], fit$fitted[j],col="red")
        y - fit$fitted
    })
    Y
}


totalProteinNormalise <- function(Y) {
    nchannels <- dim(Y)[2]
    new.Y <- sapply(1:nchannels, function(i) {
        y <- Y[,i]
        totalP <- rowSums(Y[,-i])
        fit <- loess(y ~ totalP)
        y - fit$fitted
    })
    colnames(new.Y) <- colnames(Y)
    new.Y
}


## ## crap script

## pdf("../img/corrplots_after_loess.pdf",width=8, height=8)

## load("../data/sp5.RData")
## sp5 <- sp
## load("../data/sp6.RData")
## sp6 <- sp

## sp5@readouts <- loessNormalise(sp5)
## sp6@readouts <- loessNormalise(sp6)

## channelCorr(sp5, 2)
## channelCorr(sp5, 1)
## channelCorr(sp6, 1)
## channelCorr(sp6, 2)

## dev.off()
