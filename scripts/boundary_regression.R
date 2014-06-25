library(EMCluster)
library(devtools)
library(rgl)
library(R.matlab)
library(gplots)

load_all("..")

load("~/ebi/data/spe.RData")

spid <- 5

sp <- spe[[spid]]

Y <- cells(sp)

Y <- Y - min(Y) + 1

ycor <- cor(t(Y))

y.pc <- princomp(Y)

y.scores <- y.pc$scores


d <- y.scores[,1:3]

emobj <- simple.init(d, nclass=2)
emobj <- shortemcluster(d, emobj)
ret <- emcluster(d, emobj, assign.class=TRUE)

nnids <- nnID(sp)

clust <- ret$class

cl1 <- which(clust == 1)
cl2 <- which(clust == 2)


## find out which is tumour and which is stromal
kerindex <- grep("Keratin",pNames(sp))
kerReads <- Y[,kerindex]
tumourID <- which.max( c(mean(kerReads[cl1]), mean(kerReads[cl2])))

boundary <- NULL

for(i in cl1) {
    nn <- nnids[[i]]
    if(length(intersect(nn,cl2)) > 0) boundary <- c(boundary, i)
}

for(i in cl2) {
    nn <- nnids[[i]]
    if(length(intersect(nn, cl1)) > 0) boundary <- c(boundary, i)
}

protein.names <- pNames(sp)
getProteinIds <- function(proteinNames,proteinList) {
    p.indices <- sapply(proteinList, grep, proteinNames)
    p.indices
}


responseSubset <- getProteinIds(protein.names,
                                c("Cadherin",
                                  "bcat",
                                  "Vimentin",
                                  "CD44","Twist",
                                  "Slug", "NFkB",
                                  "EGFR"))

#responseSubset <- getProteinIds(protein.names,c("Slug"))

dependentRemoveSubset <- getProteinIds(protein.names,
                                       c("Keratin",
                                         "ER", "PR",
                                         "Ki67","pSHP",
                                         "Caspase", "Panactin",
                                         "Her2", "CAH9", "CK7","H3"))

dependentRemoveSubset <- unlist(dependentRemoveSubset)

pset <- 1:length(protein.names)

pset <- setdiff(pset, dependentRemoveSubset)

cell.class <- ret$class

weight.name <- paste("../data/boundary_sizes",spid,".mat",sep="")

weights <- readMat(weight.name)
weights <- weights[[1]]
weights <- lapply(weights, as.vector)

regList <- weightedSubsetBoundaryRegression(sp, cell.class, tumourID, boundary,
                                        weights, responseSubset, pset)

fit <- regList$fit

s <- summary(fit)

TTmat <- matrix(0, ncol=length(responseSubset),nrow=length(pset))
STmat <- matrix(0, ncol=length(responseSubset),nrow=length(pset))

colnames(TTmat) <- colnames(STmat) <- protein.names[responseSubset]
rownames(TTmat) <- rownames(STmat) <- protein.names[pset]

for(i in 1:length(s)) {
    coef <- s[[i]]$coefficients
    sig <- as.numeric(coef[-1,4] < 0.05)
    mid <- length(sig)/2
    STmat[,i] <- sig#[1:mid]
    #TTmat[,i] <- sig[-(1:mid)]
}

pdfname <- paste("../img/boundary_hmap_sample",spid,".pdf",sep="")

pdf(pdfname,width=6, height=6)
#heatmap.2(TTmat,dendrogram="none",Rowv=FALSE, Colv=FALSE, trace="none", margins=c(10,10))
heatmap.2(STmat,dendrogram="none",Rowv=FALSE, Colv=FALSE, trace="none", margins=c(10,10))
dev.off()


