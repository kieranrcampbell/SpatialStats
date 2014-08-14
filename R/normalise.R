###################################################################
## Linear model based normalisation for signal against cell size ##
###################################################################


#' Data normalisation after classification
#'
#' @export
normaliseSP <- function(sp) {
    raw <- rawData(sp)

    Y <- lmNormalise(raw, size(sp), cellClass(sp), FALSE)

    Y <- totalProteinNormalise(Y, cell.classes=cellClass(sp))

    Y <- preprocess.centre(Y, cellClass(sp))

    X <- lapply(neighbourIDs(sp), function(id) {
        Y[id,]
    })

    sp@readouts <- Y
    sp@cellNeighbours <- X

    sp
}

## converts the measurements for each channel
## to N(0,1)
## todo: come back and make for more than 2 classes
preprocess.centre <- function(Y, cell.class) {
    Y <- apply(Y, 2, function(y) {
        y1 <- y[cell.class == 1]
        y1 <- (y1 - mean(y1)) / sd(y1)
        y2 <- y[cell.class == 2]
        y2 <- (y2 - mean(y2)) / sd(y2)
        y[cell.class == 1] <- y1
        y[cell.class == 2] <- y2
        y
    })
    Y
}


#' Provides column by column normalisation on Y
#' given cell size s
#'
#' @param Y Measurement by channel matrix
#' @param s Vector of cell sizes
lmNormalise <- function(Y, s, cell.classes, showPlot=FALSE) {
    if(showPlot) {
        plotLm(Y, s, cell.classes)
    }

    y <- lapply(1:2, function(i) Y[cell.classes == i,])
    sizes <- lapply(1:2, function(i) s[cell.classes == i])


    y.norm <- lapply(1:2, function(i) lmY(y[[i]], sizes[[i]]))

    Y.n <- matrix(NA, nrow=nrow(Y), ncol=ncol(Y))
    Y.n[cell.classes == 1, ] <- y.norm[[1]]
    Y.n[cell.classes == 2, ] <- y.norm[[2]]
    colnames(Y.n) <- colnames(Y)
    Y.n
}

lmY <- function(Y, s) {
    Y <- apply(Y, 2, function(y) {
        ##print(c(length(y), length(s)))
        fit <- lm(y ~ s)
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
plotLm <- function(Y, s, cell.classes) {
    y <- lapply(1:2, function(i) Y[cell.classes == i,])
    sizes <- lapply(1:2, function(i) s[cell.classes == i])

    z1 <- y[[1]] ; z2 <- y[[2]]
    l1 <- lm(z1[,1] ~ sizes[[1]])
    l2 <- lm(z2[,1] ~ sizes[[2]])

    plot(s, Y[,1], col=cell.classes, xlab="Cell size (px)", ylab="Log expr")
    j1 <- order(sizes[[1]])
    j2 <- order(sizes[[2]])
    lines(sizes[[1]][j1], l1$fitted[j1],col="black")
    lines(sizes[[2]][j2], l2$fitted[j2],col="red", cex=2)

}

#' Uses lm normalisation to normalise by cell concentration
#'
#' @param separate Normalise each cell class separately
#' @param cell.classes The cell classes as passed by cellClass
#'
totalProteinNormalise <- function(Y, separate=TRUE, cell.classes=NULL) {
    nchannels <- ncol(Y)
    if(!separate) {
        new.Y <- sapply(1:nchannels, function(i) {
            y <- Y[,i]
            totalP <- rowSums(Y[,-i])
            fit <- lm(y ~ totalP)
            y - fit$fitted
        })
        colnames(new.Y) <- colnames(Y)
        return(new.Y)
    } else {
        y <- lapply(1:2, function(i) Y[cell.classes == i,])
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
        Y.n[cell.classes == 1, ] <- y.concnorm[[1]]
        Y.n[cell.classes == 2, ] <- y.concnorm[[2]]
        colnames(Y.n) <- colnames(Y)
        return(Y.n)
    }
}

