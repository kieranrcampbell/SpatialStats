
## generates toy data for SPData and subsequent regression analysis
## four cells with 2,3,2,1 nearest neighbours

getInferredCoefficient <- function(coeff) {
    weights <- list(c(1,1),c(1,1,1),c(1,1),1)

    nchannel <- 1

    meanLogExprs <- 6.74
    sd <- 0.51

    X <- lapply(weights, function(weight) {
        matrix(rnorm(length(weight), sd), ncol=1)
    })

    X <- c(X, rep(0,4))
    weights <- c(weights, rep(1,4))

    noiseSD <- sd

    #cell.base <- c(3,7,5,6,0,0,0,0)

    Y <- sapply(X, function(x) coeff*mean(x) + rnorm(1, 0, noiseSD))
    Y <- Y
    Y <- as.matrix(Y)

    X <- lapply(X, function(x) x + rnorm(length(x),0,noiseSD))

    nn.ids <- list(5:6, 5:7, 7:8,  8)
    nn.ids <- c(nn.ids, 1:4)

    cellClasses <- c(rep(1,4), rep(2, 4))
    boundary <- which(cellClasses > 0)

    responseSubset <- dependentSubset <- 1
    responseClass <- 1

    sp <- new("SPData", protein.names="channel",
              Y = Y, X = X, nn.ids = nn.ids,
              size = rep(1,8),id=1, pos=0, nn.counts=list())

    breg <- weightedSubsetBoundaryRegression(sp, cellClasses, responseClass=1, boundary,
                                             weights, responseSubset, dependentSubset)
    return(breg$fit$coefficients[2])
}

coeff.range <- seq(0,0.5, length.out=10)

N <- 100

coeff.out <- sapply(coeff.range, function(coeff) {
    replicate(N, getInferredCoefficient(coeff))
})

d <- data.frame(coeff.out)
names(d) <- as.character(coeff.range)

dfm <- melt(d)

plt <- ggplot(aes(y=value,x=variable), data=dfm) + geom_boxplot()
ggsave("../test/boxplt.pdf", plot=plt)
