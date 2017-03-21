# Load the Immuno-Navigator Macrophage data

# How data was obtained...
# macrophageExp <- read.table("~/gd/Harvard/Research/data/ImmunoNavigator/expression_1.txt",sep="\t", header = T, stringsAsFactors = F)
# saveRDS(macrophageExp, "~/gd/Harvard/Research/data/ImmunoNavigator/macrophageExp.rds")

library(gplots)

source('~/gd/Harvard/Research/network_batch/diagnosis_plots.R')
source('~/gd/Harvard/Research/network_batch/algorithm.R')
source('~/gd/Harvard/Research/network_batch/simulateData.R')
source('~/gd/Harvard/Research/network_batch/generateMasterDF.R')
source('~/gd/Harvard/Research/network_batch/run_wgcna.R')

seed <- sample(10000,1)
set.seed(seed)

numSamples <- 600
numGenes <- 4000
batchProp <- .75

samplesA1 <- 1:(batchProp*numSamples/2)
samplesB1 <- (batchProp*numSamples/2+1):(numSamples/2)
samplesB2 <- (numSamples/2+1):(batchProp*numSamples/2 + numSamples/2)
samplesA2 <- (batchProp*numSamples/2 + numSamples/2+1):(numSamples)

batches <- rep(0,numSamples)
batches[c(samplesA1,samplesB1)]<-1
caseControls <- rep(0,numSamples)
caseControls[c(samplesA1,samplesA2)] <- 1

note <- "_firsttry"
outputDir <- paste0("figures/macrophage_",numSamples,"samples_",numGenes,"genes_",batchProp,"batchprop_",seed,note)
dir.create(outputDir,showWarnings = F)

macrophageExp <- readRDS("~/gd/Harvard/Research/data/ImmunoNavigator/macrophageExp.rds")

# Subset number of samples and scramble it to avoid any possibly confounding by index for samples or genes
macrophageExp <- macrophageExp[!is.na(rowSums(macrophageExp)),sample(ncol(macrophageExp))[1:numSamples]]
macrophageExp <- macrophageExp[sample(nrow(macrophageExp)),]

expressionData1 <- as.matrix(macrophageExp[1:numGenes,])
expressionData1 <- expressionData1-rowMeans(expressionData1)

geneModules <- getClustersFromWGCNA(expressionData1)
orderedModules <- sort(geneModules)
modulesToScramble <- c(1,2)
expressionData1Scrambled <- expressionData1[order(geneModules),]
expressionData1Original <- expressionData1[order(geneModules),]

# Recreate the truth
batchGenes <- as.numeric(orderedModules==modulesToScramble[1])
realGenes <- as.numeric(orderedModules==modulesToScramble[2])
batchMat <- tcrossprod(batchGenes) 
realMat <- tcrossprod(realGenes)

trueLabels <- rep("Background",choose(numGenes,2))
trueLabels[batchMat[row(batchMat) > col(batchMat)]!=0] <- "Batch effect"
trueLabels[realMat[row(realMat) > col(realMat)]!=0] <- "Real effect"

trueGeneLabels <- rep("Background", numGenes)
trueGeneLabels[batchGenes==1] <- "Batch"
trueGeneLabels[realGenes==1] <- "Real"

coex <- cor(t(expressionData1Original))
diag(coex) <- NA
png(paste0(outputDir,'/coex_heatmap_pam_modules.png'), width = 1600, height = 1200)
heatmap.2(coex[c(T,F,F,F),c(T,F,F,F)], Rowv = F, Colv = F, trace = "none", 
          labRow=sort(geneModules)[c(T,F,F,F)], col="bluered", dendrogram = "none", 
          RowSideColors = cbPalette[as.factor(trueGeneLabels[c(T,F,F,F)])])
dev.off()


# Uncorrelate module 1
expressionData1Scrambled[orderedModules==modulesToScramble[1], batches==0] <- 
    t(apply(expressionData1Scrambled[orderedModules==modulesToScramble[1],batches==0],1, sample))
expressionData1Scrambled[orderedModules==modulesToScramble[2], caseControls==0] <- 
    t(apply(expressionData1Scrambled[orderedModules==modulesToScramble[2],caseControls==0],1, sample))

coex <- cor(t(expressionData1Scrambled))
diag(coex) <- NA
png(paste0(outputDir,'/coex_heatmap_pam_modules_scrambledModules.png'), width = 1600, height = 1200)
heatmap.2(coex[c(T,F,F,F),c(T,F,F,F)], Rowv = F, Colv = F, trace = "none", 
          labRow=sort(geneModules)[c(T,F,F,F)], col="bluered", dendrogram = "none", 
          RowSideColors = cbPalette[as.factor(trueGeneLabels[c(T,F,F,F)])])
dev.off()

G_star <- expressionData1Original-rowMeans(expressionData1Original)
eigenG <- eigs_sym(tcrossprod((G_star/sqrt(rowSums(G_star^2)))),10)

expressionData <- expressionData1Scrambled

X <- cbind(rep(1,numSamples),batches, caseControls)


# macrophage_result <- themethod(X, as.matrix(expressionData), absolute=F, eigen_function = eigs_sym, N=100)
macrophage_result <- themethod(X, as.matrix(expressionData), absolute=F, eigen_function = eigs_sym, N=10, eigenG=eigenG)



# Post-estimation analysis ------------------------------------------------




differentialCorrelationNaive <- cor(t(macrophage_result$G_standard[,caseControls==1]))-cor(t(macrophage_result$G_standard[,caseControls==0]))
differentialCorrelationNaivewBatch <- 
    (cor(t(macrophage_result$G_standard[,caseControls==1&batches==1]))-cor(t(macrophage_result$G_standard[,caseControls==0&batches==1])) +
         cor(t(macrophage_result$G_standard[,caseControls==1&batches==0]))-cor(t(macrophage_result$G_standard[,caseControls==0&batches==0])))/2


macrophage_MasterDF <- generateMasterDF(macrophage_result, differentialCorrelationNaive, differentialCorrelationNaivewBatch, trueLabels, maxPoints=800000)

plotEigenvectors(macrophage_result, trueGeneLabels, outputDir)

mean(macrophage_MasterDF$newMeth[macrophage_MasterDF$labels=="Real effect"])
mean(macrophage_MasterDF$newMeth[macrophage_MasterDF$labels=="Batch effect"])

diagnosticPlots(macrophage_MasterDF, outputDir)



