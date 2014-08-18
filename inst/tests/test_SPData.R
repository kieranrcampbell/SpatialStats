## creates small mock dataset

ncells <- 5
nprot <- 10
sID <- 1

Y <- t(replicate(5, rpois(10,5), simplify=T))

protein.names <- paste("prot", 1:nprot, sep="")

nn.ids <- lapply(1:ncells, function(i) {
  sample(1:ncells, sample(2:4,1))
})

X <- lapply(nn.ids, function(ids) {
  Y[ids,]
})

cell.classes <- c(1,1:4)

pos <- matrix(sample(1:10, ncells * 2, TRUE), ncol=2)

sp <- SPData(channelNames=protein.names,
             readouts=Y,raw=Y, cellNeighbours=X,
             size=rep(1, ncells), id=sID, weights=list(0),
             pos=pos, cellClass=cell.classes, nn.ids = nn.ids)


test_that("SPData is initialised properly", {

    expect_that(nCells(sp), equals(ncells))
    expect_that(nChannel(sp), equals(nprot))
    expect_that(length(channels(sp)), equals(nprot))
    expect_that(ID(sp), equals(sID))
    expect_that(channels(sp), equals(protein.names))
    expect_equal(pos, xy(sp))
    expect_equal(neighbours(object = sp), neighbours(generateNeighbourfromID(sp)))
})

test_that("SPData subsets correctly", {
  ## use neighbourChannel
  cell.subset <- sample(1:nCells(sp), nCells(sp)/2)
  channel.subset <- sample(1:nChannel(sp), nChannel(sp)/2)
  expect_equal(Y[cell.subset,], cells(sp[cell.subset,]))
  expect_equal(Y[,channel.subset], cells(sp[,channel.subset]))  
  
  ## nearest neighbours - does it subset by channel correctly?
  X <- lapply(neighbours(sp), function(x) x[,channel.subset])
  expect_equal(X, neighbourChannel(neighbours(sp), channel.subset))
  
  ## does it still work with single subsetting?
  single.cell <- sample(1:nCells(sp), 1)
  single.channel <- sample(1:nChannel(sp),1)
  sp.cell <- sp[single.cell,]
  sp.channel <- sp[,single.channel]
})

test_that("SPData cell class handling", {
    expect_equal(cell.classes, cellClass(sp))   
    ## neighbourClass 
    
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
