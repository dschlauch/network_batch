library(MASS)
library(gplots)
library(ROCR)
library(ggplot2)
library(rARPACK)

source('~/gd/Harvard/Research/network_batch/diagnosis_plots.R')
source('~/gd/Harvard/Research/network_batch/algorithm.R')
source('~/gd/Harvard/Research/network_batch/simulateData.R')
source('~/gd/Harvard/Research/network_batch/generateMasterDF.R')

seed <- sample(10000,1)
set.seed(seed)
numGenes <- 4000 
numSamples <-200
addedError <- 8
batchEffectMultiplier <- 2

note <- "_testing"
outputDir <- paste0("figures/InSilico_",addedError,"Error_",batchEffectMultiplier,"Batch_UniformEffects_",numSamples,"samples_",numGenes,"genes","_",seed,note)
dir.create(outputDir,showWarnings = F)

mu <- rnorm(numGenes,mean = 9)
batches <- c(rep(0,numSamples/2),rep(1,numSamples/2))
caseControl <- c(rep(0,numSamples*4/10),rep(1,numSamples/2),rep(0,numSamples*1/10))
X <- cbind(rep(1,numSamples),batches, caseControl)
# blockSeq <- c("A","B","B","B","B","B","B","C","C","C","D","D","D","D","E","E","F","F","F","F","G","G","G","G","G","G","G")
blockSeq <- sample(LETTERS[1:10],50,replace=T)

# heatmap.2(data, trace = "none", col = "bluered")
study <- simulateStudy(numGenes=numGenes, numSamples=numSamples, addedError=addedError, 
                       blockSeq=blockSeq, mu=mu, caseControl=caseControl, batches=batches, batchEffectMultiplier)

# Recreate the truth
batchMat <- tcrossprod(study$trueEffects$batch1Effect) - tcrossprod(study$trueEffects$batch2Effect) 
realMat <- tcrossprod(study$trueEffects$casesEffect) - tcrossprod(study$trueEffects$controlsEffect)

truePairwiseLabels <- rep("Background",choose(numGenes,2))
truePairwiseLabels[batchMat[row(batchMat) > col(batchMat)]!=0] <- "Batch effect"
truePairwiseLabels[realMat[row(realMat) > col(realMat)]!=0] <- "Real effect"

trueGeneLabels <- rep("Background", numGenes)
trueGeneLabels[study$batchEffectedGenes] <- "Batch"
trueGeneLabels[study$realEffectedGenes] <- "Real"

coex <- cor(t(study$data))
diag(coex) <- NA
png(paste0(outputDir,'/coex_heatmap.png'), width = 1600, height = 1200)
heatmap.2(coex[c(T,F,F,F),c(T,F,F,F)], Rowv = F, Colv = F, trace = "none", 
          labRow=trueGeneLabels[c(T,F,F,F)], col="bluered", dendrogram = "none", RowSideColors = cbPalette[as.factor(trueGeneLabels[c(T,F,F,F)])])
dev.off()

insilico_result <- themethod(X, study$data, absolute = F, eigen_function = eigs_sym)
differentialCorrelationNaive <- cor(t(insilico_result$G_standard[,caseControl==1]))-cor(t(insilico_result$G_standard[,caseControl==0]))
differentialCorrelationNaivewBatch <- 
    (cor(t(insilico_result$G_standard[,caseControl==1&batches==1]))-cor(t(insilico_result$G_standard[,caseControl==0&batches==1])) +
    cor(t(insilico_result$G_standard[,caseControl==1&batches==0]))-cor(t(insilico_result$G_standard[,caseControl==0&batches==0])))/2

insilico_MasterDF <- generateMasterDF(insilico_result, differentialCorrelationNaive, differentialCorrelationNaivewBatch, truePairwiseLabels, maxPoints=800000)

plotEigenvectors(insilico_result, 
                 trueGeneLabels,
                 dir=outputDir, numEigenvectors=6)

mean(insilico_MasterDF$newMeth[insilico_MasterDF$labels=="Real effect"])
mean(insilico_MasterDF$newMeth[insilico_MasterDF$labels=="Batch effect"])

diagnosticPlots(insilico_MasterDF, outputDir)
# hist(differentialCorrelationNaive)