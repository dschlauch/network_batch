library(gplots)
# set.seed(1234)
numSamples <- 120

batchProp <- .75


samplesA1 <- 1:(batchProp*numSamples/2)
samplesB1 <- (batchProp*numSamples/2+1):(numSamples/2)
samplesB2 <- (numSamples/2+1):(batchProp*numSamples/2 + numSamples/2)
samplesA2 <- (batchProp*numSamples/2 + numSamples/2+1):(numSamples)


numBatchGenes <- 300
numRealGenes <- 300
numBackgroundGenes <- 1000
nBackgroundHubs<-4


eclipseExp <- read.table("~/gd/Harvard/Research/data/Eclipse/ECLIPSE_Blood_Exp.txt",row.names=1,header=T)
eclipseExp <- eclipseExp[!is.na(rowSums(eclipseExp)),sample(ncol(eclipseExp))[1:numSamples]]
expressionData <- eclipseExp[sample(nrow(eclipseExp), (2*numBatchGenes+2*numRealGenes+numBackgroundGenes)),]


backgroundIndices <- 1:numBackgroundGenes
batch1Indices <- (numBackgroundGenes+1):(numBackgroundGenes+numBatchGenes)
batch2Indices <- (numBackgroundGenes+numBatchGenes+1):(numBackgroundGenes+2*numBatchGenes)
realAIndices <- (numBackgroundGenes+2*numBatchGenes+1):(numBackgroundGenes+2*numBatchGenes+numRealGenes)
realBIndices <- (numBackgroundGenes+2*numBatchGenes+numRealGenes+1):(numBackgroundGenes+2*numBatchGenes+2*numRealGenes)

# expressionData <- cbind(
#     expressionData[c(backgroundIndices, batch1Indices, realAIndices),samplesA1],                 #A1
#     expressionData[c(backgroundIndices, batch1Indices, realBIndices),samplesB1],               #B1
#     expressionData[c(backgroundIndices, batch2Indices, realBIndices),samplesB2],                  #B2
#     expressionData[c(backgroundIndices, batch2Indices, realAIndices),samplesA2])    #A2

backgroundGEXP <- expressionData[backgroundIndices,1:numSamples]
hubs <- matrix(backgroundIndices,ncol=nBackgroundHubs)
backgroundGEXP <- do.call(rbind,apply(hubs,2,function(x){
    nd <- backgroundGEXP[x,sample(ncol(backgroundGEXP))]
    colnames(nd) <-colnames(backgroundGEXP) 
    nd
}))
batchGEXP <- cbind(
    expressionData[batch1Indices, sample(samplesA1)],                 #A1
    expressionData[batch1Indices, sample(samplesB1)],               #B1
    expressionData[batch2Indices, sample(samplesB2)],                  #B2
    expressionData[batch2Indices, sample(samplesA2)])
realGEXP <- cbind(
    expressionData[realAIndices, sample(samplesA1)],                 #A1
    expressionData[realBIndices, sample(samplesB1)],               #B1
    expressionData[realBIndices, sample(samplesB2)],                  #B2
    expressionData[realAIndices, sample(samplesA2)])
colnames(batchGEXP) <- colnames(backgroundGEXP)
colnames(realGEXP) <- colnames(backgroundGEXP)
expressionData <- rbind(backgroundGEXP, batchGEXP, realGEXP)

expressionData[,samplesA1] <- expressionData[,samplesA1]-rowMeans(expressionData[,samplesA1])
expressionData[,samplesB1] <- expressionData[,samplesB1]-rowMeans(expressionData[,samplesB1])
expressionData[,samplesB2] <- expressionData[,samplesB2]-rowMeans(expressionData[,samplesB2])
expressionData[,samplesA2] <- expressionData[,samplesA2]-rowMeans(expressionData[,samplesA2])

batches <- rep(0,numSamples)
batches[c(samplesA1,samplesB1)]<-1
caseControls <- rep(0,numSamples)
caseControls[c(samplesA1,samplesA2)] <- 1
X <- cbind(rep(1,numSamples),batches, caseControls)

batchArr <- c(rep(0,numBackgroundGenes),rep(1,numBatchGenes),rep(0,numRealGenes))
realArr <- c(rep(0,numBackgroundGenes+numBatchGenes),rep(1,numBatchGenes))
batchMat <- batchArr%*%t(batchArr)
realMat <- realArr%*%t(realArr)

source('~/gd/Harvard/Research/network_batch/algorithm.R')
eclipseRes <- themethod(X, as.matrix(expressionData), absolute=T, eigen_function = eigs_sym, N=100)



# Post-estimation analysis ------------------------------------------------





X_bar <- matrix(colMeans(X),nrow=1)
fitValues <- c(X_bar%*%eclipseRes$estimates)
correctedCorrelation <- eclipseRes$Q%*%diag(fitValues)%*%t(eclipseRes$Q)

differentialCorrelation <- eclipseRes$Q%*%diag(eclipseRes$estimates[3,])%*%t(eclipseRes$Q) 

differentialCorrelationNaive <- cor(t(eclipseRes$G_standard[,caseControls==1]))-cor(t(eclipseRes$G_standard[,caseControls==0]))
differentialCorrelationNaivewBatch <- 
    (cor(t(eclipseRes$G_standard[,caseControls==1&batches==1]))-cor(t(eclipseRes$G_standard[,caseControls==0&batches==1])) +
    cor(t(eclipseRes$G_standard[,caseControls==1&batches==0]))-cor(t(eclipseRes$G_standard[,caseControls==0&batches==0])))/2


# batchMat <- as.numeric(study$blocks=="B") %*% t(as.numeric(study$blocks=="B"))
# realMat <- as.numeric(study$blocks=="C") %*% t(as.numeric(study$blocks=="C"))

labels <- rep("Background",choose(nrow(eclipseRes$G_standard),2))
labels[batchMat[row(batchMat) > col(batchMat)]==1] <- "Batch effect"
labels[realMat[row(realMat) > col(realMat)]==1] <- "Real effect"

batch1 <- abs(cor(t(as.matrix(eclipseRes$G_standard)[,c(samplesA1,samplesB1)])))
batch2 <- abs(cor(t(as.matrix(eclipseRes$G_standard)[,c(samplesA2,samplesB2)])))
real1 <- abs(cor(t(as.matrix(eclipseRes$G_standard)[,c(samplesA1,samplesA2)])))
real2 <- abs(cor(t(as.matrix(eclipseRes$G_standard)[,c(samplesB1,samplesB2)])))
group1A <- abs(cor(t(as.matrix(eclipseRes$G_standard)[,samplesA1])))
coex <- cor(t(as.matrix(eclipseRes$G_standard)))

diag(batch1) <- NA
diag(batch2) <- NA
diag(real1) <- NA
diag(real2) <- NA
diag(group1A) <- NA
diag(coex) <- NA

# heatmap.2(coex, Rowv = F, Colv = F, trace = "none", labRow=c(rep("-----",numBackgroundGenes),rep("-",numBatchGenes),rep("-----",numRealGenes)), col="bluered", dendrogram = "none")
# heatmap.2(batch1, Rowv = F, Colv = F, trace = "none", labRow=c(rep("-----",numBackgroundGenes),rep("-",numBatchGenes),rep("-----",numRealGenes)), col="bluered", dendrogram = "none")
# heatmap.2(batch2, Rowv = F, Colv = F, trace = "none", labRow=c(rep("-----",numBackgroundGenes),rep("-",numBatchGenes),rep("-----",numRealGenes)), col="bluered", dendrogram = "none")
# heatmap.2(real1, Rowv = F, Colv = F, trace = "none", labRow=c(rep("-----",numBackgroundGenes),rep("-",numBatchGenes),rep("-----",numRealGenes)), col="bluered", dendrogram = "none")
# heatmap.2(real2, Rowv = F, Colv = F, trace = "none", labRow=c(rep("-----",numBackgroundGenes),rep("-",numBatchGenes),rep("-----",numRealGenes)), col="bluered", dendrogram = "none")
# heatmap.2(group1A, Rowv = F, Colv = F, trace = "none", labRow=c(rep("-----",numBackgroundGenes),rep("-",numBatchGenes),rep("-----",numRealGenes)), col="bluered", dendrogram = "none")
# heatmap.2(differentialCorrelationNaive, Rowv = F, Colv = F, trace = "none", labRow=c(rep("-----",numBackgroundGenes),rep("-",numBatchGenes),rep("-----",numRealGenes)), col="bluered", dendrogram = "none")
# heatmap.2(differentialCorrelation, Rowv = F, Colv = F, trace = "none", labRow=c(rep("-----",numBackgroundGenes),rep("-",numBatchGenes),rep("-----",numRealGenes)), col="bluered", dendrogram = "none")

allCorrelations <- data.frame(newMeth=differentialCorrelation[row(differentialCorrelation) > col(differentialCorrelation)],
                              naiveMeth=differentialCorrelationNaive[row(differentialCorrelationNaive) > col(differentialCorrelationNaive)],
                              naiveWBatch=differentialCorrelationNaivewBatch[row(differentialCorrelationNaivewBatch) > col(differentialCorrelationNaivewBatch)],
                              labels=labels)
allCorrelations <- allCorrelations[sample(nrow(allCorrelations)),]

# ggplot(allCorrelations[1:50000,]) + geom_point(aes(x=naiveMeth,y=newMeth, color=factor(labels)), alpha=.5, size=1) +
# 2/7/17 note

diagnosticPlots(allCorrelations, "figures/ECLIPSE_findBatch")

