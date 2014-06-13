
## pca & cluster analysis

library(rgl)
library(devtools)
library(gplots)
library(RColorBrewer)

load_all("..")

load("../data/spe.RData")

cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)

exp <- 13
sp <- spe[[exp]]

fit <- nnReg(sp)
m <- getPMat(sp, fit, setOne=TRUE, alpha=0.01)

Y <- cells(sp)

y.pc <- princomp(Y)

y.scores <- y.pc$scores

k <- kmeans(y.scores[,1:3],2)
k <- pam(y.scores[,1:3], 2, diss=FALSE)

## use k$cluster

library(EMCluster)

## d <- y.scores[,1:3]
## emobj <- simple.init(d, nclass=2)
## emobj <- shortemcluster(d, emobj)
## ret <- emcluster(d, emobj, assign.class=TRUE)

## use ret$class

sp1 <- sp[which(k$cluster==1),]
sp2 <- sp[which(k$cluster==2),]

fit1 <- nnReg(sp1) ; fit2 <- nnReg(sp2)

s1 <- summary(fit1) ; s2 <- summary(fit2)

m1 <- getPMat(sp1, fit1, setOne=TRUE, alpha=0.01)
m2 <- getPMat(sp2, fit2, setOne=TRUE, alpha=0.01)


pdf("../img/pca.pdf")

heatmap.2(m, trace="none", margins=c(10,10), Rowv=FALSE, Colv=FALSE)
heatmap.2(m1, trace="none", margins=c(10,10))
heatmap.2(m2, trace="none", margins=c(10,10))

dev.off()

## em stuff


