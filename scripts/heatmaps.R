library("devtools")
library(RColorBrewer)

load_all("..")

load("../data/spe.RData")


library(gplots)
cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)

control <- which(ids(spe) %in% c(343, 359))

## mean data heatmap

all.data <- sapply(data(spe), function(sp) {
    Y <- cells(sp)
    colMeans(Y)
})

colnames(all.data) <- as.character(ids(spe))
rownames(all.data) <- pNames(spe[[1]])

png("../img/heatmaps/heatmap_mean_all.png",height=700,width=800)
heatmap.2(all.data, trace="none",scale="col",col=rev(cols), margins=c(5,10),
          main="Mean protein level vs sample", xlab="sample",ylab="protein")
dev.off()


## sample correlation heatmap

samp.cor <- cor(all.data)

s.names <- as.character(ids(spe))

s.names <- gsub("343", "343 (non tumour)", s.names)
s.names <- gsub("359", "359 (non tumour)", s.names)

rownames(samp.cor) <- colnames(samp.cor) <- s.names


png("../img/heatmaps/heatmap_mean_samples.png",height=700,width=800)
heatmap.2(samp.cor, trace="none", col=rev(cols), symm=TRUE, margins=c(10,10),
          main="Correlation between samples \n across mean protein expression")
dev.off()

## protein correlation heatmap

prot.cor <- cor(t(all.data))

p.names <- as.character(pNames(spe[[1]]))

#p.names <- gsub("343", "343 (non tumour)", p.names)
#p.names <- gsub("359", "359 (non tumour)", p.names)

rownames(prot.cor) <- colnames(prot.cor) <- p.names

png("../img/heatmaps/heatmap_mean_proteins.png",height=700,width=800)
heatmap.2(prot.cor, trace="none", col=rev(cols), symm=TRUE, margins=c(10,10),
          main="Correlation between mean protein \n expression across samples")
dev.off()



protein.names <- pNames(spe[[1]])

sig.data <- sapply(data(spe), function(sp) {
    fit <- nnReg(sp)
    coef <- fit$coefficients ## (p+1) by p matrix
    rowMeans(abs(coef[2:33,]))
})

rownames(sig.data) <- pNames(spe[[1]])
colnames(sig.data) <- as.character(ids(spe))


pp.data <- sapply(data(spe), function(sp) {
    fit <- nnReg(sp)
    coef <- fit$coefficients ## (p+1) by p matrix
    #coef[1,] <- 0
    coef[2:33,]
})


dim(pp.data) <- c(32,32,21)

all.means <- apply(pp.data, c(1,2), mean) ## average pathway activation
colnames(all.means) <- protein.names
rownames(all.means) <- c(protein.names)

control.means <- apply(pp.data[,,control], c(1,2), mean) ## average pathway activation
colnames(control.means)  <- protein.names
rownames(control.means)  <- protein.names

tumor.means <- apply(pp.data[,,-control], c(1,2), mean) ## average pathway activation
colnames(tumor.means) <- rownames(tumor.means) <- protein.names

cols <- colorRampPalette(brewer.pal(10, "RdBu"))(256)


png("../img/heatmaps/all_hmap.png", width=800, height=800)
heatmap.2(all.means, trace="none", margins=c(10,10), col=rev(cols), scale=c("col"),
          main="Pathway activation, all samples",ylab="Nearest neighbour protein",
          xlab="Cell protein")
dev.off()

png("../img/heatmaps/control_hmap.png", width=800, height=800)
heatmap.2(control.means, trace="none", margins=c(10,10), scale="col", col=rev(cols),
          main="Pathway activation \n healthy tissue", ylab="Nearest neighbour protein",
          xlab="Cell protein")
dev.off()

png("../img/heatmaps/tumor_hmap.png", width=800, height=800)
heatmap.2(tumor.means, trace="none", margins=c(10,10), scale="col", col=rev(cols),
          main="Pathway activation \n tumour only", ylab="Nearest neighbour protein",
          xlab="Cell protein")
dev.off()

## single cell correlations

## control samples
Y.control <- rbind(cells(spe[[control [1] ]]) / size(spe[[ control[1] ]]),
                   cells(spe[[control [2] ]]) / size(spe[[ control[2] ]]))
ycontrol.cor <- cor(Y.control)
rownames(ycontrol.cor) <- colnames(ycontrol.cor) <- protein.names

Y.t <- NULL
for(i in setdiff(1:(length(files(spe))), control)) {
    Y.t <-rbind(Y.t,cells(spe[[i]]) / size(spe[[i]]))

}

tumor.corr <- cor(Y.t)
rownames(tumor.corr) <- colnames(tumor.corr) <- protein.names

pdf("../img/heatmaps/heatmaps.pdf")
heatmap.2(all.data, trace="none",scale="col",col=rev(cols), margins=c(5,10),
          main="Mean protein \n level vs sample", xlab="sample",ylab="protein")
heatmap.2(samp.cor, trace="none", col=rev(cols), symm=TRUE, margins=c(10,10),
          main="Correlation between samples \n across mean protein expression")
heatmap.2(prot.cor, trace="none", col=rev(cols), symm=TRUE, margins=c(10,10),
          main="Correlation between mean protein \n expression across samples")

heatmap.2(ycontrol.cor, trace="none", col=rev(cols), symm=TRUE, margins=c(10,10),
          main="  Correlation between proteins \n  across single cells (normal)")
heatmap.2(tumor.corr, trace="none", col=rev(cols), symm=TRUE, margins=c(10,10),
          main="  Correlation between proteins \n  across single cells (tumor)")



heatmap.2(all.means, trace="none", margins=c(10,10), col=rev(cols), scale="none",
          main="Pathway activation \n all samples",ylab="Nearest neighbour protein",
          xlab="Cell protein")
heatmap.2(control.means, trace="none", margins=c(10,10), scale="none", col=rev(cols),
          main="Pathway activation \n healthy tissue", ylab="Nearest neighbour protein",
          xlab="Cell protein")
heatmap.2(tumor.means, trace="none", margins=c(10,10), scale="none", col=rev(cols),
          main="Pathway activation \n tumour only", ylab="Nearest neighbour protein",
          xlab="Cell protein")
dev.off()

