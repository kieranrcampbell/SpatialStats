
############################################################
## LOESS based normalisation for signal against cell size ##
## kieran.campbell@sjc.ox.ac.uk                           ##
############################################################



#' Data normalisation after classification
#'
#' @export
normaliseSP <- function(sp) {
    raw <- rawData(sp)
    Y <- loessNormalise(raw, size(sp), cellClass(sp), TRUE)
    Y <- totalProteinNormalise(Y)
    Y <- preprocess.centre(Y)

    X <- lapply(neighbourIDs(sp), function(id) {
        Y[id,]
    })

    sp@readouts <- Y
    sp@cellNeighbours <- X

    sp
}


#' Provides column by column normalisation on Y
#' given cell size s
#'
#' @param Y Measurement by channel matrix
#' @param s Vector of cell sizes
loessNormalise <- function(Y, s, cell.classes, showPlot=FALSE) {
    if(showPlot) {
        plotLoess(Y, s, cell.classes)
    }

    y <- lapply(1:2, function(i) Y[cell.classes == i,])
    sizes <- lapply(1:2, function(i) s[cell.classes == i])


    y.norm <- lapply(1:2, function(i) loessY(y[[i]], sizes[[i]]))

    Y.n <- matrix(0, nrow=nrow(Y), ncol=ncol(Y))
    Y.n[cell.classes == 1, ] <- y.norm[[1]]
    Y.n[cell.classes == 2, ] <- y.norm[[2]]
    colnames(Y.n) <- colnames(Y)
    Y.n
}

loessY <- function(Y, s) {
    Y <- apply(Y, 2, function(y) {
        ##print(c(length(y), length(s)))
        fit <- loess(y ~ s)
        ##j <- order(s)
        ##plot(fit)
        ##lines(s[j], fit$fitted[j],col="red")
        y - fit$fitted
    })
    Y
}

#' y & s are lists
#'
#'
plotLoess <- function(Y, s, cell.classes) {
    y <- lapply(1:2, function(i) Y[cell.classes == i,])
    sizes <- lapply(1:2, function(i) s[cell.classes == i])

    z1 <- y[[1]] ; z2 <- y[[2]]
    l1 <- loess(z1[,1] ~ sizes[[1]])
    l2 <- loess(z2[,1] ~ sizes[[2]])

    plot(s, Y[,1], col=cell.classes, xlab="Cell size (px)", ylab="Log expr")
    j1 <- order(sizes[[1]])
    j2 <- order(sizes[[2]])
    lines(sizes[[1]][j1], l1$fitted[j1],col="black")
    lines(sizes[[2]][j2], l2$fitted[j2],col="red", cex=2)

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

