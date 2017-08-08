# ENCODE
library(varhandle)
library(ggplot2)
library(reshape2)
library(sva)
library(limma)
library(dplyr)

set.seed(123)
nperm <- 50
plotSubset <- c(T,rep(F,100))
topNgenes <- 30000

load("~/gd/Harvard/Research/data/ENCODE/encode_lcl_rnaseq_counts.RData")
dat <- dat[rowSums(dat!=0)>152,]
expr <- unfactor(dat)
# select top genes by variance, remove replicates
expr <- expr[order(-apply(expr,1,var))[seq_len(topNgenes)],!grepl("_2_",colnames(dat))]

yaleBatch <- grepl(pattern = "yale",colnames(expr))
expr_combat = ComBat(dat=expr, batch=yaleBatch, par.prior=TRUE, prior.plots=FALSE)


makeDiffCoexDF <- function(exprData){
    
    coex <- cor(t(exprData[,yaleBatch])) - cor(t(exprData[,!yaleBatch]))
    coexvalues <- sort(abs(coex[row(coex)>col(coex)]))
    
    nullMatrix <- replicate(nperm,{
        group1 <- sample(c(rep(T,32),rep(F,31)))
        yaleBatchPerm <- c(rbind(group1, !group1)) 
        coex_perm <- abs(cor(t(exprData[,yaleBatchPerm])) - cor(t(exprData[,!yaleBatchPerm])))
        sort(coex_perm[row(coex_perm)>col(coex_perm)])
    })
    colnames(nullMatrix) <- paste0("null",seq_len(ncol(nullMatrix)))
    diffCoex <- melt(cbind(nullMatrix, expected=rowMeans(nullMatrix), coexvalues=coexvalues))
    as.data.frame(diffCoex)[plotSubset,]
}

diffCoexSubsetRaw    <- cbind(makeDiffCoexDF(expr),type="Uncorrected")
diffCoexSubsetComBat <- cbind(makeDiffCoexDF(expr_combat),type="After Batch Correction")

diffCoexSubset <- rbind(diffCoexSubsetRaw,diffCoexSubsetComBat)

pdf(paste0("./figures/encode_diff_coex_density.pdf"), width=12, height=6)
densplot <- ggplot(diffCoexSubset) + 
    geom_density(aes(x=value, group=Var2),color="grey",alpha=.05,fill="lightgrey") +
    geom_density(data=diffCoexSubset[diffCoexSubset$Var2=="coexvalues",], aes(x=value,color="blue"),fill="blue", alpha=.2) + 
    geom_density(data=diffCoexSubset[diffCoexSubset$Var2=="expected",], aes(x=value,color="black")) + 
    facet_wrap(~type) +
    ggtitle("Differential Coexpression") + xlab("Absolute Differential Coexpression") + ylab("Density") +
    theme_bw(base_size = 20) + 
    theme(plot.title = element_text(hjust = 0.5),legend.justification=c(1,0),
          legend.position=c(.95,.75),
          legend.background = element_rect(fill="gray95", size=.5, linetype="dotted"))+
    scale_colour_manual(name = 'Partition',
                        values =c('blue'='blue','black'='black'), labels = c('Mean Random','Across Batches')) 
print(densplot)
dev.off()

### LIMMA


# Do all genes for LIMMA
expr <- unfactor(dat)
# Remove strange outlier gene
expr <- expr[-which(rownames(expr)=="ENSG00000075624.9"),]
expr <- expr[,!grepl("_2_",colnames(dat))]
yaleBatch <- grepl(pattern = "yale",colnames(expr))
expr_combat = ComBat(dat=expr, batch=yaleBatch, par.prior=TRUE, prior.plots=FALSE)
design <- model.matrix(~yaleBatch)
colnames(design) <- c("Intercept","Batch2-Batch5")
fit <- lmFit(expr, design)
fit_combat <- lmFit(expr_combat, design)
fit2 <- eBayes(fit)
fit2_combat <- eBayes(fit_combat)
output <- topTable(fit2, coef = 2,n=Inf)
output_combat <- topTable(fit2_combat, coef = 2,n=Inf)

sum(output$adj.P.Val<.001)
sum(output_combat$adj.P.Val<.001)

rawDF <- data.frame(pvalues=output$P.Value, 
                logFC=output$logFC, type="Uncorrected")
combatDF <- data.frame(pvalues=output_combat$P.Value,
                logFC=output_combat$logFC, type="After Batch Correction") 
df <- rbind(rawDF,combatDF)

pdf(paste0("./figures/encode_diffexpress.pdf"), width=12, height=6)
diffexpressPlot <- ggplot(df) + 
    geom_point(aes(x=logFC,y=-log10(pvalues)), shape=3) + 
    facet_wrap(~type) + 
    theme_bw(base_size = 20) + theme(plot.title = element_text(hjust = 0.5)) + 
    ggtitle("Differential Expression") + xlab("Log Fold Change") + ylab("-Log p-values")
plot(diffexpressPlot)
dev.off()

