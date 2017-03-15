library(gplots)

source('~/gd/Harvard/Research/network_batch/diagnosis_plots.R')
source('~/gd/Harvard/Research/network_batch/algorithm.R')
source('~/gd/Harvard/Research/network_batch/simulateData.R')
source('~/gd/Harvard/Research/network_batch/generateMasterDF.R')
source('~/gd/Harvard/Research/network_batch/run_wgcna.R')

seed <- sample(10000,1)
set.seed(seed)

numSamples <- 200
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

# numBatchGenes <- 600
# numRealGenes <- 600
# numBackgroundGenes <- 2800
# nBackgroundHubs<-4
# numGenes <- numBackgroundGenes+numBatchGenes+numRealGenes

note <- "_modulescrambling"
outputDir <- paste0("figures/ECLIPSE_",numSamples,"samples_",numGenes,"genes_",batchProp,"batchprop_",seed,note)
dir.create(outputDir,showWarnings = F)

# Load the ECLIPSE data
eclipseExp <- read.table("~/gd/Harvard/Research/data/Eclipse/ECLIPSE_Blood_Exp.txt",row.names=1,header=T)

# Subset number of samples and scramble it to avoid any possibly confounding by index for samples or genes
eclipseExp <- eclipseExp[!is.na(rowSums(eclipseExp)),sample(ncol(eclipseExp))[1:numSamples]]
eclipseExp <- eclipseExp[sample(nrow(eclipseExp)),]

expressionData1 <- as.matrix(eclipseExp[1:numGenes,])
expressionData1 <- expressionData1-rowMeans(expressionData1)
# expressionData2 <- as.matrix(eclipseExp[(numGenes+1):(2*numGenes),])
# expressionData2 <- expressionData2-rowMeans(expressionData2)

# backgroundIndices <- 1:numBackgroundGenes
# batchIndices <- (numBackgroundGenes+1):(numBackgroundGenes+numBatchGenes)
# realIndices <- (numBackgroundGenes+numBatchGenes+1):(numBackgroundGenes+numBatchGenes+numRealGenes)

# expressionData <- matrix(NA, nrow=nrow(expressionData1), ncol=ncol(expressionData1))
# expressionData[backgroundIndices, c(T,F)] <- expressionData1[backgroundIndices, c(T,F)]
# expressionData[backgroundIndices, c(F,T)] <- expressionData2[backgroundIndices, c(F,T)]
# expressionData[batchIndices, batches==0] <- expressionData1[batchIndices, batches==0]
# expressionData[batchIndices, batches==1] <- expressionData2[batchIndices, batches==1]
# expressionData[realIndices, caseControls==0] <- expressionData1[realIndices, caseControls==0]
# expressionData[realIndices, caseControls==1] <- expressionData2[realIndices, caseControls==1]

geneModules <- getClustersFromWGCNA(expressionData1)
orderedModules <- sort(geneModules)
modulesToScramble <- c(1,2)
expressionData1Scrambled <- expressionData1[order(geneModules),]
expressionData1Original <- expressionData1[order(geneModules),]

coex <- cor(t(expressionData1Original))
diag(coex) <- NA
png(paste0(outputDir,'/coex_heatmap_pam_modules.png'), width = 1600, height = 1200)
heatmap.2(coex[c(T,F,F,F),c(T,F,F,F)], Rowv = F, Colv = F, trace = "none", labRow=sort(geneModules)[c(T,F,F,F)], col="bluered", dendrogram = "none")
dev.off()

# Uncorrelate module 1
expressionData1Scrambled[orderedModules==modulesToScramble[1], batches==0] <- 
    t(apply(expressionData1Scrambled[orderedModules==modulesToScramble[1],batches==0],1, sample))
expressionData1Scrambled[orderedModules==modulesToScramble[2], caseControls==0] <- 
    t(apply(expressionData1Scrambled[orderedModules==modulesToScramble[2],caseControls==0],1, sample))

coex <- cor(t(expressionData1Scrambled))
diag(coex) <- NA
png(paste0(outputDir,'/coex_heatmap_pam_modules_scrambledModules.png'), width = 1600, height = 1200)
heatmap.2(coex[c(T,F,F,F),c(T,F,F,F)], Rowv = F, Colv = F, trace = "none", labRow=sort(geneModules)[c(T,F,F,F)], col="bluered", dendrogram = "none")
dev.off()

# if (FALSE){
    G_star <- expressionData1Original-rowMeans(expressionData1Original)
    eigenG <- eigs_sym(tcrossprod((G_star/sqrt(rowSums(G_star^2)))),200)
# }
expressionData <- expressionData1Scrambled
# expressionData <- cbind(
#     expressionData[c(backgroundIndices, batch1Indices, realAIndices),samplesA1],                 #A1
#     expressionData[c(backgroundIndices, batch1Indices, realBIndices),samplesB1],               #B1
#     expressionData[c(backgroundIndices, batch2Indices, realBIndices),samplesB2],                  #B2
#     expressionData[c(backgroundIndices, batch2Indices, realAIndices),samplesA2])    #A2

# backgroundGEXP <- expressionData[backgroundIndices,1:numSamples]
# hubs <- matrix(backgroundIndices,ncol=nBackgroundHubs)
# shuffle <- identity
# backgroundGEXP <- do.call(rbind,apply(hubs,2,function(x){
#     nd <- backgroundGEXP[x,sample(ncol(backgroundGEXP))]
#     colnames(nd) <-colnames(backgroundGEXP) 
#     nd
# }))
# batchGEXP <- cbind(
#     expressionData[batch1Indices, shuffle(samplesA1)],                 #A1
#     expressionData[batch1Indices, shuffle(samplesB1)],               #B1
#     expressionData[batch2Indices, shuffle(samplesB2)],                  #B2
#     expressionData[batch2Indices, shuffle(samplesA2)])
# realGEXP <- cbind(
#     expressionData[realAIndices, shuffle(samplesA1)],                 #A1
#     expressionData[realBIndices, shuffle(samplesB1)],               #B1
#     expressionData[realBIndices, shuffle(samplesB2)],                  #B2
#     expressionData[realAIndices, shuffle(samplesA2)])
# colnames(batchGEXP) <- colnames(backgroundGEXP)
# colnames(realGEXP) <- colnames(backgroundGEXP)
# expressionData <- rbind(backgroundGEXP, batchGEXP, realGEXP)
# 
# expressionData[,samplesA1] <- expressionData[,samplesA1]-rowMeans(expressionData[,samplesA1])
# expressionData[,samplesB1] <- expressionData[,samplesB1]-rowMeans(expressionData[,samplesB1])
# expressionData[,samplesB2] <- expressionData[,samplesB2]-rowMeans(expressionData[,samplesB2])
# expressionData[,samplesA2] <- expressionData[,samplesA2]-rowMeans(expressionData[,samplesA2])

X <- cbind(rep(1,numSamples),batches, caseControls)
# 
# batchArr <- c(rep(0,numBackgroundGenes),rep(1,numBatchGenes),rep(0,numRealGenes))
# realArr <- c(rep(0,numBackgroundGenes+numBatchGenes),rep(1,numBatchGenes))
# batchMat <- batchArr%*%t(batchArr)
# realMat <- realArr%*%t(realArr)

# eclipse_result <- themethod(X, as.matrix(expressionData), absolute=F, eigen_function = eigs_sym, N=100)
eclipse_result <- themethod(X, as.matrix(expressionData), absolute=F, eigen_function = eigs_sym, N=100, eigenG=eigenG)



# Post-estimation analysis ------------------------------------------------




differentialCorrelationNaive <- cor(t(eclipse_result$G_standard[,caseControls==1]))-cor(t(eclipse_result$G_standard[,caseControls==0]))
differentialCorrelationNaivewBatch <- 
    (cor(t(eclipse_result$G_standard[,caseControls==1&batches==1]))-cor(t(eclipse_result$G_standard[,caseControls==0&batches==1])) +
    cor(t(eclipse_result$G_standard[,caseControls==1&batches==0]))-cor(t(eclipse_result$G_standard[,caseControls==0&batches==0])))/2

# Recreate the truth
# batchGenes <- c(rep(0,numBackgroundGenes), rep(1,numBatchGenes), rep(0,numRealGenes))
# realGenes <- c(rep(0,numBackgroundGenes), rep(0,numBatchGenes), rep(1,numRealGenes))
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

eclipse_MasterDF <- generateMasterDF(eclipse_result, differentialCorrelationNaive, differentialCorrelationNaivewBatch, trueLabels, maxPoints=800000)

plotEigenvectors(eclipse_result, trueGeneLabels, outputDir)

mean(eclipse_MasterDF$newMeth[eclipse_MasterDF$labels=="Real effect"])
mean(eclipse_MasterDF$newMeth[eclipse_MasterDF$labels=="Batch effect"])

diagnosticPlots(eclipse_MasterDF, outputDir)



# 
# ######################
# 
# 
# X_bar <- matrix(colMeans(X),nrow=1)
# fitValues <- c(X_bar%*%eclipseRes$estimates)
# correctedCorrelation <- eclipseRes$Q%*%diag(fitValues)%*%t(eclipseRes$Q)
# 
# differentialCorrelation <- eclipseRes$Q%*%diag(eclipseRes$estimates[3,])%*%t(eclipseRes$Q) 
# 
# differentialCorrelationNaive <- cor(t(eclipseRes$G_standard[,caseControls==1]))-cor(t(eclipseRes$G_standard[,caseControls==0]))
# differentialCorrelationNaivewBatch <- 
#     (cor(t(eclipseRes$G_standard[,caseControls==1&batches==1]))-cor(t(eclipseRes$G_standard[,caseControls==0&batches==1])) +
#     cor(t(eclipseRes$G_standard[,caseControls==1&batches==0]))-cor(t(eclipseRes$G_standard[,caseControls==0&batches==0])))/2
# 
# 
# # batchMat <- as.numeric(study$blocks=="B") %*% t(as.numeric(study$blocks=="B"))
# # realMat <- as.numeric(study$blocks=="C") %*% t(as.numeric(study$blocks=="C"))
# 
# labels <- rep("Background",choose(nrow(eclipseRes$G_standard),2))
# labels[batchMat[row(batchMat) > col(batchMat)]==1] <- "Batch effect"
# labels[realMat[row(realMat) > col(realMat)]==1] <- "Real effect"
# 
# batch1 <- abs(cor(t(as.matrix(eclipseRes$G_standard)[,c(samplesA1,samplesB1)])))
# batch2 <- abs(cor(t(as.matrix(eclipseRes$G_standard)[,c(samplesA2,samplesB2)])))
# real1 <- abs(cor(t(as.matrix(eclipseRes$G_standard)[,c(samplesA1,samplesA2)])))
# real2 <- abs(cor(t(as.matrix(eclipseRes$G_standard)[,c(samplesB1,samplesB2)])))
# group1A <- abs(cor(t(as.matrix(eclipseRes$G_standard)[,samplesA1])))
# coex <- cor(t(as.matrix(eclipseRes$G_standard)))
# 
# diag(batch1) <- NA
# diag(batch2) <- NA
# diag(real1) <- NA
# diag(real2) <- NA
# diag(group1A) <- NA
# diag(coex) <- NA
# 
# # heatmap.2(coex, Rowv = F, Colv = F, trace = "none", labRow=c(rep("-----",numBackgroundGenes),rep("-",numBatchGenes),rep("-----",numRealGenes)), col="bluered", dendrogram = "none")
# # heatmap.2(batch1, Rowv = F, Colv = F, trace = "none", labRow=c(rep("-----",numBackgroundGenes),rep("-",numBatchGenes),rep("-----",numRealGenes)), col="bluered", dendrogram = "none")
# # heatmap.2(batch2, Rowv = F, Colv = F, trace = "none", labRow=c(rep("-----",numBackgroundGenes),rep("-",numBatchGenes),rep("-----",numRealGenes)), col="bluered", dendrogram = "none")
# # heatmap.2(real1, Rowv = F, Colv = F, trace = "none", labRow=c(rep("-----",numBackgroundGenes),rep("-",numBatchGenes),rep("-----",numRealGenes)), col="bluered", dendrogram = "none")
# # heatmap.2(real2, Rowv = F, Colv = F, trace = "none", labRow=c(rep("-----",numBackgroundGenes),rep("-",numBatchGenes),rep("-----",numRealGenes)), col="bluered", dendrogram = "none")
# # heatmap.2(group1A, Rowv = F, Colv = F, trace = "none", labRow=c(rep("-----",numBackgroundGenes),rep("-",numBatchGenes),rep("-----",numRealGenes)), col="bluered", dendrogram = "none")
# # heatmap.2(differentialCorrelationNaive, Rowv = F, Colv = F, trace = "none", labRow=c(rep("-----",numBackgroundGenes),rep("-",numBatchGenes),rep("-----",numRealGenes)), col="bluered", dendrogram = "none")
# # heatmap.2(differentialCorrelation, Rowv = F, Colv = F, trace = "none", labRow=c(rep("-----",numBackgroundGenes),rep("-",numBatchGenes),rep("-----",numRealGenes)), col="bluered", dendrogram = "none")
# 
# allCorrelations <- data.frame(newMeth=differentialCorrelation[row(differentialCorrelation) > col(differentialCorrelation)],
#                               naiveMeth=differentialCorrelationNaive[row(differentialCorrelationNaive) > col(differentialCorrelationNaive)],
#                               naiveWBatch=differentialCorrelationNaivewBatch[row(differentialCorrelationNaivewBatch) > col(differentialCorrelationNaivewBatch)],
#                               labels=labels)
# allCorrelations <- allCorrelations[sample(nrow(allCorrelations)),]
# 
# # ggplot(allCorrelations[1:50000,]) + geom_point(aes(x=naiveMeth,y=newMeth, color=factor(labels)), alpha=.5, size=1) +
# # 2/7/17 note
# 
# diagnosticPlots(allCorrelations, "figures/ECLIPSE_findBatch")

