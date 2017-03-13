library(MASS)
library(gplots)
library(ROCR)
library(ggplot2)
library(rARPACK)

source('~/gd/Harvard/Research/network_batch/diagnosis_plots.R')
source('~/gd/Harvard/Research/network_batch/algorithm.R')
source('~/gd/Harvard/Research/network_batch/simulateData.R')
source('~/gd/Harvard/Research/network_batch/generateMasterDF.R')

outputDir <- "figures/InSilico_extraBatch"
seed <- 0
set.seed(seed)
numGenes <- 4000 
numSamples <-400
mu <- rnorm(numGenes,mean = 9)
batches <- c(rep(0,numSamples/2),rep(1,numSamples/2))
addedError <- 20
caseControl <- c(rep(0,numSamples*3/8),rep(1,numSamples/2),rep(0,numSamples*1/8))
X <- cbind(rep(1,numSamples),batches, caseControl)
blockSeq <- c("A","B","B","B","B","B","B","C","C","C","D","D","D","D","E","E","F","F","F","F","G","G","G","G","G","G","G")
blockSeq <- sample(LETTERS[1:10],50,replace=T)

# heatmap.2(data, trace = "none", col = "bluered")
study <- simulateStudy(numGenes=numGenes, numSamples=numSamples, addedError=addedError, 
                       blockSeq=blockSeq, mu=mu, caseControl=caseControl, batches=batches)

# Recreate the truth
batchMat <- tcrossprod(study$trueEffects$batch1Effect) - tcrossprod(study$trueEffects$batch2Effect) 
realMat <- tcrossprod(study$trueEffects$casesEffect) - tcrossprod(study$trueEffects$controlsEffect)

trueLabels <- rep("Background",choose(numGenes,2))
trueLabels[batchMat[row(batchMat) > col(batchMat)]==1] <- "Batch effect"
trueLabels[realMat[row(realMat) > col(realMat)]>0] <- "Real effect"
trueLabels[realMat[row(realMat) > col(realMat)]<0] <- "Negative effect"

insilico_result <- themethod(X, study$data, absolute = F, eigen_function = eigs_sym)
differentialCorrelationNaive <- cor(t(insilico_result$G_standard[,caseControl==1]))-cor(t(insilico_result$G_standard[,caseControl==0]))
differentialCorrelationNaivewBatch <- 
    cor(t(insilico_result$G_standard[,caseControl==1&batches==1]))-cor(t(insilico_result$G_standard[,caseControl==0&batches==1])) +
    cor(t(insilico_result$G_standard[,caseControl==1&batches==0]))-cor(t(insilico_result$G_standard[,caseControl==0&batches==0]))

insilico_MasterDF <- generateMasterDF(insilico_result, differentialCorrelationNaive, differentialCorrelationNaivewBatch, trueLabels, maxPoints=800000)

plotEigenvectors(insilico_result, study$trueEffects, outputDir)


diagnosticPlots(insilico_MasterDF, outputDir)
# hist(differentialCorrelationNaive)