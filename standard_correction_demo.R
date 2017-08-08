library(MASS)
library(sva)
library(limma)


numGenes <- 1000
numBatch1 <- 100
numBatch2 <- 100
module1 <- c(rep(1,200),rep(0,800))
module2 <- c(rep(0,200),rep(1,200),rep(0,600))
meanExpression1 <- rnorm(numGenes, 7)
meanExpression2 <- rnorm(numGenes, 7)
nperm <- 100
plotSubset <- c(T,rep(F,10))

sigma1 <- tcrossprod(module1)
diag(sigma1) <- sample(4:15,numGenes,replace = T)
sigma2 <- tcrossprod(module2)
diag(sigma2) <- sample(4:15,numGenes,replace = T)
data <- t(rbind(mvrnorm(numBatch1, mu = meanExpression1, Sigma = sigma1),mvrnorm(numBatch2, mu=meanExpression2, Sigma=sigma2)))

copdData      <- read.table("~/gd/Harvard/Research/data/COPDGene/COPDGene_GSExpressionData.txt",row.names=1,header=T)
data <- as.matrix(cbind(copdData[sample(nrow(copdData),numGenes),1:numBatch1],copdData[sample(nrow(copdData),numGenes),1:numBatch2]))

# mod0 = model.matrix(~rep(1,numBatch1+numBatch2))[,1,drop=F]
# n.sv = num.sv(data, mod0, method="leek")
# svobj = sva(data, mod0, n.sv=n.sv)

batches <- c(rep("Batch1",numBatch1),rep("Batch2",numBatch2))
expr_combat = ComBat(dat=data, batch=batches, par.prior=TRUE, prior.plots=FALSE)



makeDiffCoexDF <- function(exprData){
    coex <- cor(t(exprData[,batches=="Batch1"])) - cor(t(exprData[,batches=="Batch2"]))
    coexvalues <- sort(abs(coex[row(coex)>col(coex)]))
    
    nullMatrix <- replicate(nperm,{
        cat(".")
        batchPerm <- sample(c(rep(T,numBatch2),rep(F,numBatch2)))
        coex_perm <- abs(cor(t(exprData[,batchPerm])) - cor(t(exprData[,!batchPerm])))
        sort(coex_perm[row(coex_perm)>col(coex_perm)])
    })
    colnames(nullMatrix) <- paste0("null",seq_len(ncol(nullMatrix)))
    diffCoex <- melt(cbind(nullMatrix, expected=rowMeans(nullMatrix), coexvalues=coexvalues))
    as.data.frame(diffCoex)[plotSubset,]
}

diffCoexSubsetRaw    <- cbind(makeDiffCoexDF(data),type="Uncorrected")
diffCoexSubsetComBat <- cbind(makeDiffCoexDF(expr_combat),type="After Batch Correction")

diffCoexSubset <- rbind(diffCoexSubsetRaw,diffCoexSubsetComBat)

pdf(paste0("./figures/demo_diff_coex_density.pdf"), width=12, height=6)
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
design <- model.matrix(~batches)
colnames(design) <- c("Intercept","Batch2")
fit <- lmFit(data, design)
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

pdf(paste0("./figures/demo_diffexpress.pdf"), width=12, height=6)
diffexpressPlot <- ggplot(df) + 
    geom_point(aes(x=logFC,y=-log10(pvalues)), shape=3) + 
    facet_wrap(~type) + 
    theme_bw(base_size = 20) + theme(plot.title = element_text(hjust = 0.5)) + 
    ggtitle("Differential Expression") + xlab("Log Fold Change") + ylab("-Log p-values")
plot(diffexpressPlot)
dev.off()
