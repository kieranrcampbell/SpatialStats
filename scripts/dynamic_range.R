## plot the dynamic ranges

library("devtools")

load_all("..")

load("../data/spe.RData")

n.exp <- length(ids(spe))

protein.names <- pNames(spe[[1]])

if(n.exp != 21) stop("Need 21 experiments for plotting")

pdf("../img/drange.pdf", height=8, width=12)
par(mfrow=c(3,7), mar=c(1.3,1.3,1.3,1.3))

for(sp in data(spe)) {
    Y <- as.vector(cells(sp))

    boxplot(Y,main=as.character(id(sp)))

}




par(mfrow=c(3,7))

for(sp in data(spe)) {
    Y <- as.vector(cells(sp))

    hist(Y,main=as.character(id(sp)))

}

for(sp in data(spe)) {
    par(mfrow=c(4,8))
    Y <- cells(sp)
    for(i in 1:32) {
        .y <- Y[,i]
        boxplot(.y, main=protein.names[i], xlab="", ylab="", ylim=c(5,10))
    }
    title(as.character(id(sp)), outer=TRUE)
}


dev.off()


