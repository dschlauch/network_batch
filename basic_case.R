# test case

library(MASS)
source('~/gd/Harvard/Research/network_batch/algorithm.R')
source('~/gd/Harvard/Research/network_batch/diagnosis_plots.R')

g=1000
n=2000
outputDir <- "figures/simple7_nonstandardized"
dir.create(outputDir)

Sigma1aVec <- c(rep(1,g/5),rep(0,g/5),rep(0,g/5),rep(0,g/5),rep(0,g/5))
Sigma1a <- tcrossprod(Sigma1aVec)
diag(Sigma1a) <- 1

Sigma2aVec <- cbind(c(rep(1,g/5),rep(0,g/5),rep(0,g/5),rep(0,g/5),rep(0,g/5)),
                    c(rep(0,g/5),rep(1,g/5),rep(0,g/5),rep(0,g/5),rep(0,g/5)))
Sigma2a <- tcrossprod(Sigma2aVec)
diag(Sigma2a) <- 1

Sigma1bVec <- c(rep(0,g/5),rep(1,g/5),rep(0,g/5),rep(0,g/5),rep(0,g/5))
Sigma1b <- tcrossprod(Sigma1bVec)
diag(Sigma1b) <- 1

Sigma2bVec <- c(rep(0,g/5),rep(0,g/5),rep(0,g/5),rep(0,g/5),rep(0,g/5))
Sigma2b <- tcrossprod(Sigma2bVec)
diag(Sigma2b) <- 1


x1a <- mvrnorm(n=n, mu=rep(0,g),Sigma = Sigma1a)
x2a <- mvrnorm(n=n, mu=rep(0,g),Sigma = Sigma2a)
x1b <- mvrnorm(n=n, mu=rep(0,g),Sigma = Sigma1b)
x2b <- mvrnorm(n=n, mu=rep(0,g),Sigma = Sigma2b)
data <- t(rbind(x1a,x2a,x1b,x2b))

trueGeneLabels <- c(rep("Real",g*.4),rep("Batch",g*.4),rep("BG",g*.2))
X <- cbind(1,c(rep(0,2*n),rep(1,2*n)),c(rep(0,n),rep(1,n),rep(0,n),rep(1,n)))

simpleRes <- themethod(X,data,N = 6, standardize = F)

differentialCorrelationNaive <- round(cor(t(data[,X[,2]==1]))-cor(t(data[,X[,2]==0])),4)

coex <- cor(t(data))
diag(coex) <- NA
png(paste0(outputDir,'/coex_heatmap.png'), width = 1600, height = 1200)
heatmap.2(coex, Rowv = F, Colv = F, trace = "none", 
          labRow=trueGeneLabels, col="bluered", dendrogram = "none", RowSideColors = cbPalette[as.factor(trueGeneLabels)])
dev.off()

plotEigenvectors(simpleRes, 
                 trueGeneLabels,
                 dir=outputDir, numEigenvectors=6)

differentialCorrelation <- simpleRes$Q%*%diag(simpleRes$estimates[2,])%*%t(simpleRes$Q)

png(paste0(outputDir,'/est_diff_coex.png'), width = 1600, height = 1200)
heatmap.2(differentialCorrelation, Rowv = F, Colv = F, trace = "none", 
          labRow=trueGeneLabels, col="bluered", dendrogram = "none", RowSideColors = cbPalette[as.factor(trueGeneLabels)])
dev.off()


png(paste0(outputDir,'/naive_est_diff_coex.png'), width = 1600, height = 1200)
heatmap.2(differentialCorrelationNaive, Rowv = F, Colv = F, trace = "none", 
          labRow=trueGeneLabels, col="bluered", dendrogram = "none", RowSideColors = cbPalette[as.factor(trueGeneLabels)])
dev.off()
differentialCorrelation[1:5,1:5]
# mean(insilico_MasterDF$newMeth[insilico_MasterDF$labels=="Real effect"])
# mean(insilico_MasterDF$newMeth[insilico_MasterDF$labels=="Batch effect"])
# 
# diagnosticPlots(insilico_MasterDF, outputDir)
