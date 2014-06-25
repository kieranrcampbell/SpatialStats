library(EMCluster)
library(devtools)
library(rgl)

load_all("..")

load("~/ebi/data/spe.RData")

sp <- spe[[5]]

Y <- cells(sp)

Y <- Y - min(Y) + 1
#Y <- Y / size(sp)


ycor <- cor(t(Y))

y.pc <- princomp(Y)

y.scores <- y.pc$scores


d <- y.scores[,1:3]
#d <- Y
emobj <- simple.init(d, nclass=2)
emobj <- shortemcluster(d, emobj)
ret <- emcluster(d, emobj, assign.class=TRUE)

nnids <- nnID(sp)

clust <- ret$class

cl1 <- which(clust == 1)
cl2 <- which(clust == 2)

boundary <- NULL

for(i in cl1) {
    nn <- nnids[[i]]
    if(length(intersect(nn,cl2)) > 0) boundary <- c(boundary, i)
}

for(i in cl2) {
    nn <- nnids[[i]]
    if(length(intersect(nn, cl1)) > 0) boundary <- c(boundary, i)
}

colours <- ret$class
colours[boundary] <- 3

sp.bound <- sp[boundary,]
Y.b <- cells(sp.bound)
colnames(Y.b) <- pNames(sp.bound)

Y <- Y/size(sp)

protein.names <- pNames(sp)

ecad <- grep("Cadherin", protein.names)
bcat <- grep("bcat", protein.names)
vim <- grep("Vimentin", protein.names)
cd44 <- grep("CD44", protein.names)


cell.class <- ret$class
cell.class[intersect(boundary, which(cell.class == 1))] <- 3
cell.class[intersect(boundary, which(cell.class == 2))] <- 4

nc <- nCells(sp)


emt.counts <- as.vector(Y[,c(vim, bcat, ecad, cd44)])

sizes <- size(sp)
sz <- sizes
sizes <- (sizes - min(sizes)) / (max(sizes) - min(sizes))
#sizes <- 3*sizes



emt.counts <- c(emt.counts, sizes)
labels <- c(rep("Vimentin",nc),rep("BCatenin",nc),
            rep("ECadherin",nc), rep("CD44", nc),
            rep("Cellsize",nc))
labels <- as.factor(labels)
p.data <- data.frame(counts=emt.counts, label=labels, class=as.factor(cell.class))


slug <- grep("Slug", protein.names)
ker <- grep("Keratin", protein.names)
twist <- grep("Twist", protein.names)
nfkb <- grep("NFkB", protein.names)
egfr <- grep("EGFR", protein.names)

emt.counts2 <- as.vector(Y[,c(slug, ker, twist, nfkb, egfr)])

labels2 <- c(rep("Slug",nc),rep("Keratin",nc),rep("Twist",nc),
             rep("NFkB", nc), rep("EGFR", nc))
labels2 <- as.factor(labels2)

p.data2 <- data.frame(counts=emt.counts2, label=labels2, class=as.factor(cell.class))

library(ggplot2)
pdf("../img/emt.pdf", width=10, height=4)
ggplot(aes(y=counts, x=label, fill=class), data=p.data) + geom_boxplot() +
    xlab("") + theme_bw()
ggplot(aes(y=counts, x=label, fill=class), data=p.data2) + geom_boxplot() +
    xlab("") + theme_bw()

dev.off()




