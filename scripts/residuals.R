
## calculate residuals
library(devtools)

install("..")



load("../data/spe.RData")


sp <- spe[[ind]]

fit <- nnReg(sp)

res <- residuals(fit)

pdf("../img/residuals.pdf",height=9,width=12)

for(i in 1:length(files(spe))) {
    sp <- spe[[i]]

    fit <- nnReg(sp)

    res <- residuals(fit)

    par(mfrow=c(4,8))
    for(i in 1:32) {
        hist(res[,i],xlab="",ylab="",main=as.character(i))
    }
    title(paste("\n", as.character(id(sp))), outer=TRUE)
}

dev.off()
