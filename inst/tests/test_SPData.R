## creates small mock dataset

ncells <- 5
nprot <- 10
sID <- 1

Y <- t(replicate(5, rpois(10,5), simplify=T))

protein.names <- paste("prot", 1:nprot, sep="")

X <- lapply(1:ncells, function(i) {
    Y[sample(1:ncells, sample(2:4,1)),]

})

pos <- matrix(sample(1:10, ncells * 2, TRUE), ncol=2)

sp <- SPData(channelNames=protein.names,
             readouts=Y,raw=Y, cellNeighbours=X,
             size=rep(1, ncells), id=sID, weights=list(0),
             pos=pos, cellClass=-1)


test_that("SPData is initialised properly", {

    expect_that(nCells(sp), equals(ncells))
    expect_that(nChannel(sp), equals(nprot))
    expect_that(length(channels(sp)), equals(nprot))
    expect_that(ID(sp), equals(sID))
    expect_that(channels(sp), equals(protein.names))
    expect_equal(pos, xy(sp))
})

test_that("Nearest neighbour averages work", {

    ## simple average
    X.avg <- sapply(X, colMeans)
    X.avg <- t(X.avg)

    X.fromSP <- neighbourMean(sp, FALSE, FALSE)
    colnames(X.fromSP) <- NULL

    expect_that(as.vector(X.fromSP), equals(as.vector(X.avg)))

    w <- sapply(X, function(x) {
        nn <- nrow(x)
        rep(1 / nn, nn)
    })

    ## equal weights currently
    weight(sp) <- w
    Xw.fromSP <- neighbourMean(sp, TRUE, FALSE)
    expect_that(as.vector(Xw.fromSP), equals(as.vector(X.avg)))

    ## now test normalisation
    X.avg <- apply(X.avg, 2, function(x) (x - mean(x)) / sd(x))
    Xw.fromSP.N <- neighbourMean(sp, TRUE, TRUE)
    expect_that(as.vector(X.avg), equals(as.vector(Xw.fromSP.N)))

})
