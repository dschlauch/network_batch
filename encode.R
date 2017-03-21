# ENCODE
library(varhandle)
library(ggplot2)
library(reshape2)
library(sva)
library(limma)
load("~/gd/Harvard/Research/data/ENCODE/encode_lcl_rnaseq_counts.RData")

# source('~/gd/Harvard/Research/network_batch/algorithm.R')

head(anno)
colnames(dat)
dat <- dat[rowSums(dat!=0)>152,]

expr <- unfactor(dat)

expr <- expr[order(-apply(expr,1,var))[1:1000],]

yaleBatch <- grepl(pattern = "yale",colnames(dat))

expr = ComBat(dat=expr, batch=yaleBatch, par.prior=TRUE, prior.plots=FALSE)

design <- model.matrix(~yaleBatch)

colnames(design) <- c("Intercept","Batch2-Batch5")
fit <- lmFit(expr, design)
fit2 <- eBayes(fit)
output <- topTable(fit2, coef = 2,n=Inf)
output$P.Value
plot(-log(seq_len(length(output$P.Value))/(length(output$P.Value))),-log(output$P.Value) )
abline(0,1)

coex <- cor(t(expr[,yaleBatch])) - cor(t(expr[,!yaleBatch]))
coexvalues <- sort(abs(coex[row(coex)>col(coex)]))
nperm <- 100

# Plot null values against observed values

diffCoex <- cbind(melt(replicate(nperm,{
    yaleBatchPerm <- sample(yaleBatch)
    coex_perm <- abs(cor(t(expr[,yaleBatchPerm])) - cor(t(expr[,!yaleBatchPerm])))
    sort(coex_perm[row(coex_perm)>col(coex_perm)])
})), coexvalues)
plotSubset <- c(T,rep(F,1000))
plotSubset <- T
png("./figures/encode_null_dist.png")
ggplot(as.data.frame(diffCoex)[plotSubset,], aes(x=coexvalues)) + 
    geom_point(aes(y=value), alpha=.5, color="blue") + 
    geom_abline(intercept = 0)
dev.off()

# Plot null values AND Observed values
diffCoex <- melt(cbind(replicate(nperm,{
    yaleBatchPerm <- sample(yaleBatch)
    coex_perm <- abs(cor(t(expr[,yaleBatchPerm])) - cor(t(expr[,!yaleBatchPerm])))
    sort(coex_perm[row(coex_perm)>col(coex_perm)])
}), coexvalues))
plotSubset <- c(T,rep(F,1000))
plotSubset <- T
png("./figures/encode_null_dist_rank.png", width=800, height=800)
ggplot(as.data.frame(diffCoex)[plotSubset,], aes(x=Var1)) + 
    geom_point(aes(y=value, color=(Var2))) + 
    geom_abline(intercept = 0) + theme_bw()
dev.off()

## Density Plot
nullMatrix <- replicate(nperm,{
    yaleBatchPerm <- sample(yaleBatch)
    coex_perm <- abs(cor(t(expr[,yaleBatchPerm])) - cor(t(expr[,!yaleBatchPerm])))
    sort(coex_perm[row(coex_perm)>col(coex_perm)])
})
colnames(nullMatrix) <- paste0("null",seq_len(ncol(nullMatrix)))
diffCoex <- melt(cbind(nullMatrix, expected=rowMeans(nullMatrix), coexvalues=coexvalues))
diffCoexSubset <- as.data.frame(diffCoex)[plotSubset,]
png("./figures/encode_null_dist_density.png", width=800, height=800)
ggplot(diffCoexSubset) + 
    geom_density(aes(x=value, group=Var2),color="grey",alpha=.05,fill="lightgrey") +
    geom_density(data=diffCoexSubset[diffCoexSubset$Var2=="coexvalues",], aes(x=value,color="blue"),fill="blue", alpha=.2) + 
    geom_density(data=diffCoexSubset[diffCoexSubset$Var2=="expected",], aes(x=value,color="black")) + 
    theme_bw() + theme(plot.title = element_text(hjust = 0.5)) + ggtitle("Distribution of Differential Coexpression after ComBat") +
    scale_colour_manual(name = 'Partition',
                        values =c('blue'='blue','black'='black'), labels = c('Mean Random','Across Batches'))
dev.off()

nullMatrix <- replicate(nperm,{
    yaleBatchPerm <- sample(yaleBatch)
    coex_perm <- abs(cor(t(expr[,yaleBatchPerm])) - cor(t(expr[,!yaleBatchPerm])))
    sort(coex_perm[row(coex_perm)>col(coex_perm)])
})
diffCoex <- cbind(nullMatrix, expected=rowMeans(nullMatrix), coexvalues, rank=seq_len(nrow(nullMatrix)))
colnames(diffCoex)[seq_len(ncol(nullMatrix))] <- paste0("null",seq_len(ncol(nullMatrix)))

png("./figures/encode_null_dist_meanperm.png", width=800, height=800)
ggplot(as.data.frame(diffCoex)) + 
    geom_point(aes(x=expected, y=coexvalues),color="blue") +
    geom_point(aes(x=expected, y=null1),color="black",alpha=.1) + 
    geom_point(aes(x=expected, y=null2),color="black",alpha=.1) + 
    geom_point(aes(x=expected, y=null3),color="black",alpha=.1) + 
    geom_point(aes(x=expected, y=null4),color="black",alpha=.1) + 
    geom_point(aes(x=expected, y=null5),color="black",alpha=.1) + 
    geom_point(aes(x=expected, y=null6),color="black",alpha=.1) + 
    geom_abline(intercept = 0) + theme_bw()
dev.off()

png("./figures/encode_null_dist_meanperm_rank.png", width=800, height=800)
ggplot(as.data.frame(diffCoex)) + 
    geom_point(aes(x=rank,y=coexvalues),color="blue") +
    geom_point(aes(x=rank,y=expected),color="black",alpha=.1) +
    # geom_point(aes(x=rank, y=null1),color="red",alpha=.1) + 
    # geom_point(aes(x=rank, y=null2),color="red",alpha=.1) + 
    geom_abline(intercept = 0) + theme_bw()
dev.off()


df <- replicate(100,{rnorm(20,0)}) %>% cbind(rnorm(20,1)) %>% melt
ggplot(df) + 
    geom_density(aes(x=value, color=factor(Var2)),alpha=.1) +
    geom_density(data=df[df$Var2==101,], aes(x=value),color="black")
