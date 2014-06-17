
## pca & cluster analysis

library(rgl)
library(devtools)
library(gplots)
library(RColorBrewer)
library(EMCluster)


load_all("..")

load("../data/spe.RData")

cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)

exp <- 5
sp <- spe[[exp]]

protein.names <- pNames(sp)

fit <- nnReg(sp)
m <- getPMat(sp, fit, setOne=TRUE, alpha=0.01)

Y <- cells(sp)

y.pc <- princomp(Y)

y.scores <- y.pc$scores

## k <- kmeans(y.scores[,1:3],2)
## k <- pam(y.scores[,1:3], 2, diss=FALSE)

## use k$cluster

d <- y.scores[,1:3]
#d <- Y
emobj <- simple.init(d, nclass=2)
emobj <- shortemcluster(d, emobj)
ret <- emcluster(d, emobj, assign.class=TRUE)

## use ret$class

sp1 <- sp[which(ret$class==1),]
sp2 <- sp[which(ret$class==2),]

fit1 <- nnReg(sp1) ; fit2 <- nnReg(sp2)

s1 <- summary(fit1) ; s2 <- summary(fit2)

m1 <- getPMat(sp1, fit1, setOne=TRUE, alpha=0.01)
m2 <- getPMat(sp2, fit2, setOne=TRUE, alpha=0.01)


## pdf("../img/pca.pdf")

## heatmap.2(m, trace="none", margins=c(10,10), Rowv=FALSE, Colv=FALSE)
## heatmap.2(m1, trace="none", margins=c(10,10))
## heatmap.2(m2, trace="none", margins=c(10,10))

## dev.off()
## ## em stuff

mean1 <- colMeans(cells(sp1))
mean2 <- colMeans(cells(sp2))

names(mean1) <- names(mean2) <- protein.names

k.index <- grep("Keratin", protein.names)

library(ggplot2)

k1 <- data.frame(count=cells(sp1)[,k.index])
k2 <- data.frame(count=cells(sp2)[,k.index])

k1$nm <- 'k1'
k2$nm <- 'k2'

kdist <- rbind(k1, k2)

ggplot(kdist, aes(count, fill = nm)) + geom_density(alpha = 0.2)

md <- mean1 - mean2
mrk <- md[abs(md) > log(2)]

par(mar=c(10,5,10,5))
barplot(mrk, las=2)
#text(cex=1, y=-1.25, names(mrk), xpd=TRUE, srt=45)

