# Bladder analysis

library(bladderbatch)
library(pamr)
library(limma)
library(sva)
library(ggplot2)
library(reshape2)
library(varhandle)

data(bladderdata)

pheno = pData(bladderEset)
edata = exprs(bladderEset)
mod = model.matrix(~as.factor(cancer), data=pheno)
mod0 = model.matrix(~1,data=pheno)

n.sv = num.sv(edata,mod,method="leek")
svobj = sva(edata,mod,mod0,n.sv=n.sv)

batch = pheno$batch
modcombat = model.matrix(~1, data=pheno)

combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
# 
# pValuesComBat = f.pvalue(combat_edata,mod,mod0)
# qValuesComBat = p.adjust(pValuesComBat,method="BH")
# 
# plot(-log(seq_len(length(pValuesComBat))/length(pValuesComBat)),-log(sort(pValuesComBat)))
# abline(0,1)     


subsetBatch25Cancer <- batch%in%c(2,5) & pheno$cancer=="Cancer"
edataSubset <- combat_edata[,subsetBatch25Cancer]
pdataSubset <- pheno[subsetBatch25Cancer,]
design <- model.matrix(~pdataSubset$batch)

colnames(design) <- c("Intercept","Batch2-Batch5")
fit <- lmFit(edataSubset, design)
fit2 <- eBayes(fit)
output <- topTable(fit2, coef = 2,n=Inf)
output$P.Value
plot(-log(seq_len(nrow(edataSubset))/(nrow(edataSubset))),-log(output$P.Value) )
abline(0,1)
contrast.matrix <- cbind("C1"=c(-1,1,0,rep(0,svobj$n.sv)),"C2"=c(0,-1,1,rep(0,svobj$n.sv)),"C3"=c(-1,0,1,rep(0,svobj$n.sv)))
fitContrasts = contrasts.fit(fit,contrast.matrix)

batches <- pdataSubset$batch==2
combat_edata_subset <- combat_edata[c(T,rep(F,20)),]
coexDiff <- cor(t(combat_edata_subset[,batches])) - cor(t(combat_edata_subset[,!batches]))
coexvalues <- sort(abs(coexDiff[row(coexDiff)>col(coexDiff)]))

nperm <- 10
diffCoex <- cbind(melt(replicate(nperm,{
    Batch1Perm <- sample(batches)
    coex_perm <- abs(cor(t(combat_edata_subset[,Batch1Perm])) - cor(t(combat_edata_subset[,!Batch1Perm])))
    sort(coex_perm[row(coex_perm)>col(coex_perm)])
})), coexvalues)
plotSubset <- c(T,rep(F,1000))
plotSubset <- T
png("./figures/bladder_null_dist.png")
ggplot(as.data.frame(diffCoex)[plotSubset,], aes(x=coexvalues)) + 
    geom_point(aes(y=value), alpha=.5, color="blue") + 
    geom_abline(intercept = 0)
dev.off()


diffCoex <- melt(cbind(replicate(nperm,{
    Batch1Perm <- sample(batches)
    coex_perm <- abs(cor(t(combat_edata_subset[,Batch1Perm])) - cor(t(combat_edata_subset[,!Batch1Perm])))
    sort(coex_perm[row(coex_perm)>col(coex_perm)])
}), coexvalues))
plotSubset <- c(T,rep(F,1000))
plotSubset <- T
png("./figures/bladder_null_dist_rank.png", width=800, height=800)
ggplot(as.data.frame(diffCoex)[plotSubset,], aes(x=Var1)) + 
    geom_point(aes(y=value, color=(Var2))) + 
    geom_abline(intercept = 0) + theme_bw()
dev.off()

