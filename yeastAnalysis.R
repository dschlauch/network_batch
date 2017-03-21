## Yeast Data
library(ggplot2)
library(reshape2)
library(varhandle)
yeastExp <- read.table("../data/CompileData/YeastData_NormalizedExpression.txt",sep="\t", header = T, stringsAsFactors = F)
rownames(yeastExp) <- make.unique(yeastExp$Name, sep = ".")
yeastExp <- yeastExp[,-1]
dim(yeastExp)

# Subset data
# yeastExp <- yeastExp[c(T,F),]
# 
# 
# batch1Coex <- cor(t(yeastExp[,c(T,F)]))
# batch2Coex <- cor(t(yeastExp[,c(F,T)]))
# 
# diffCoex <- batch1Coex-batch2Coex
# 
# plot(as.numeric(yeastExp[6,]))
# plot(as.numeric(yeastExp[6,rep(1:25,each=2) + rep(c(0,25),25)]))
# 
# designMatrix <- cbind(1,rep(1:25,each=2),rep(0:1,25))
# 
# naiveCor1 <- cor(t(yeastExp[c(T,F),1:20]))
# naiveCor2 <- cor(t(yeastExp[c(T,F),31:50]))
# source('~/gd/Harvard/Research/network_batch/algorithm.R')
# yeastResult <- themethod(designMatrix,yeastExp[c(T,F),])
# 
# summary(rowSums(yeastExp))
# naiveCorDif <- naiveCor1-naiveCor2
# naiveCorDifVec <- naiveCorDif[row(naiveCorDif)>col(naiveCorDif)]
# 
# differentialCorrelation <- yeastResult$Q%*%diag(yeastResult$estimates[3,])%*%t(yeastResult$Q) 
# differentialCorrelationVec <- differentialCorrelation[row(differentialCorrelation)>col(differentialCorrelation)]
# length(differentialCorrelationVec)
# plot(differentialCorrelationVec[1:10000],naiveCorDifVec[1:10000])


# Existence of batch effect
# top 2000 genes
yeastExpTop2000 <- as.matrix(yeastExp[order(-apply(yeastExp,1,var)),])
batch1Assign <- rep(c(T,F),25)
coexDiff <- cor(t(yeastExp[,batch1Assign])) - cor(t(yeastExp[,!batch1Assign]))
coexvalues <- sort(abs(coexDiff[row(coexDiff)>col(coexDiff)]))
nperm <- 10
diffCoex <- cbind(melt(replicate(nperm,{
    Batch1Perm <- c(replicate(25,{sample(c(T,F))}))
    coex_perm <- abs(cor(t(yeastExpTop2000[,Batch1Perm])) - cor(t(yeastExpTop2000[,!Batch1Perm])))
    sort(coex_perm[row(coex_perm)>col(coex_perm)])
})), coexvalues)

plotSubset <- c(T,rep(F,1000))
plotSubset <- T
png("./figures/yeast_null_dist.png")
ggplot(as.data.frame(diffCoex)[plotSubset,], aes(x=coexvalues)) + 
    geom_point(aes(y=value), alpha=1) + 
    geom_abline(intercept = 0)
dev.off()


# Plot null values AND Observed values
diffCoex <- melt(cbind(replicate(nperm,{
    Batch1Perm <- c(replicate(25,{sample(c(T,F))}))
    coex_perm <- abs(cor(t(yeastExpTop2000[,Batch1Perm])) - cor(t(yeastExpTop2000[,!Batch1Perm])))
    sort(coex_perm[row(coex_perm)>col(coex_perm)])
}), coexvalues))
plotSubset <- c(T,rep(F,1000))
plotSubset <- T
png("./figures/yeast_null_dist_rank.png", width=800, height=800)
ggplot(as.data.frame(diffCoex)[plotSubset,], aes(x=Var1)) + 
    geom_point(aes(y=value, color=(Var2))) + 
    geom_abline(intercept = 0) + theme_bw()
dev.off()
