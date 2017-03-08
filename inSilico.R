library(MASS)
library(gplots)
library(ROCR)
library(ggplot2)


set.seed(1234)
numGenes <- 2000 
numSamples <-400
mu <- rnorm(numGenes,mean = 9)
batches <- c(rep(0,numSamples/2),rep(1,numSamples/2))
addedError <- 20
caseControl <- c(rep(0,numSamples*3/8),rep(1,numSamples/2),rep(0,numSamples*1/8))
X <- cbind(rep(1,numSamples),batches, caseControl)
# caseControl <- c(rep(0,numSamples*1/4),rep(1,numSamples/2),rep(0,numSamples*1/4))

simulateStudy <- function(numGenes, numSamples, addedError){
    start <- Sys.time()
    blockSeq <- c("A","B","B","B","B","B","B","C","C","C","D","D","D","D","E","E","F","F","F","F","G","G","G","G","G","G","G")
    blocks <- rep(blockSeq, length.out=numGenes)
    blocks <- sort(blocks)
    blockA <- as.numeric(blocks=="A")
    blockB <- as.numeric(blocks=="B") # Batch 1 only
    blockC <- as.numeric(blocks=="C") # Control only
    blockD <- as.numeric(blocks=="D")
    # blockE <- as.numeric(blocks=="E") # Batch 2 only
    blockF <- as.numeric(blocks=="F") # Cases only
    blockG <- as.numeric(blocks=="G") 
    
    batch1 <- cbind(blockA, blockB, blockD)
    batch2 <- cbind(blockA, blockD)#, blockE)      
    Sigma1 <- batch1%*%t(batch1)
    Sigma2 <- batch2%*%t(batch2)
    SigmaCase <- blockF %*% t(blockF) + blockD%*%t(blockD) - blockF%*%t(blockD)
    # SigmaCase <- diag(numGenes)
    SigmaControl <- blockC %*% t(blockC)
    
    
    Sigmas <- list(case1=Sigma1+SigmaCase, 
                   case2=Sigma2+SigmaCase, 
                   control1=Sigma1+SigmaControl, 
                   control2=Sigma2+SigmaControl)
    Sigmas <- lapply(Sigmas, function(x){
        diag(x) <- addedError # adding a bit of random error
        x
    })
    counts <- list(case1=sum(batches==0&caseControl==0),
                   case2=sum(batches==1&caseControl==0),
                   control1=sum(batches==0&caseControl==1),
                   control2=sum(batches==1&caseControl==1))
    data <- cbind(t(mvrnorm(counts[['case1']],mu=mu, Sigma = Sigmas[['case1']])),
                  t(mvrnorm(counts[['control1']],mu=mu, Sigma = Sigmas[['control1']])),
                  t(mvrnorm(counts[['control2']],mu=mu, Sigma = Sigmas[['control2']])),
                  t(mvrnorm(counts[['case2']],mu=mu, Sigma = Sigmas[['case2']])))
    data[data<0] <-0
    
    print(paste("Data generation in",round(as.numeric(difftime(Sys.time(), start,units = "secs")),1), "seconds"))
    list(data=data, blocks=blocks)
    
}
# heatmap.2(data, trace = "none", col = "bluered")
study <- simulateStudy(numGenes=numGenes, numSamples=numSamples, addedError=addedError)


source('~/gd/Harvard/Research/network_batch/algorithm.R')
result <- themethod(X, study$data, absolute = F)

plot(result$D)
evsDF <- data.frame(result$Q)
names(evsDF) <- paste0("EV",1:ncol(evsDF))
evsDF$Gene <- 1:numGenes
evsDF$blocks <- study$blocks
ggplot(evsDF) + geom_point(aes(x=Gene,y=EV1, shape=blocks))
ggplot(evsDF) + geom_point(aes(x=Gene,y=EV2, shape=blocks))
ggplot(evsDF) + geom_point(aes(x=Gene,y=EV3, shape=blocks))
ggplot(evsDF) + geom_point(aes(x=Gene,y=EV4, shape=blocks))
ggplot(evsDF) + geom_point(aes(x=Gene,y=EV5, shape=blocks))

plot(result$estimates[1,1:20])
plot(result$estimates[2,1:20])
plot(result$estimates[3,1:20])

plot(result$D[1:20])

X_bar <- matrix(colMeans(X),nrow=1)
fitValues <- c(X_bar%*%result$estimates)
correctedCorrelation <- result$Q%*%diag(fitValues)%*%t(result$Q)
plot(result$D[1:50])
plot(fitValues[1:50])
# heatmap.2(correctedCorrelation)

differentialCorrelation <- result$Q%*%diag(result$estimates[3,])%*%t(result$Q) 

differentialCorrelationNaive <- cor(t(result$G_standard[,caseControl==1]))-cor(t(result$G_standard[,caseControl==0]))
differentialCorrelationNaivewBatch <- 
    cor(t(result$G_standard[,caseControl==1&batches==1]))-cor(t(result$G_standard[,caseControl==0&batches==1])) +
    cor(t(result$G_standard[,caseControl==1&batches==0]))-cor(t(result$G_standard[,caseControl==0&batches==0]))

# Recreate the truth
batchMat <- as.numeric(study$blocks=="B") %*% t(as.numeric(study$blocks=="B"))
realMat <- as.numeric(study$blocks=="C") %*% t(as.numeric(study$blocks=="C"))
realMat2 <- as.numeric(study$blocks%in%c("D")) %*% t(as.numeric(study$blocks%in%c("D"))) + as.numeric(study$blocks%in%c("F")) %*% t(as.numeric(study$blocks%in%c("F")))
realMatNeg <- as.numeric(study$blocks%in%c("F")) %*% t(as.numeric(study$blocks%in%c("D")))
realMat[realMat>0] <- 1

labels <- rep("Background",choose(numGenes,2))
labels[batchMat[row(batchMat) > col(batchMat)]==1] <- "Batch effect"
labels[realMat[row(realMat) > col(realMat)]==1] <- "Real effect"
labels[realMat2[row(realMat2) > col(realMat2)]==1] <- "Real effect 2"
labels[realMatNeg[row(realMatNeg) > col(realMatNeg)]==1] <- "Negative effect"

allCorrelations <- data.frame(newMeth=differentialCorrelation[row(differentialCorrelation) > col(differentialCorrelation)],
                              naiveMeth=differentialCorrelationNaive[row(differentialCorrelationNaive) > col(differentialCorrelationNaive)],
                              naiveWBatch=differentialCorrelationNaivewBatch[row(differentialCorrelationNaivewBatch) > col(differentialCorrelationNaivewBatch)],
                              labels=labels)

coex <- cor(t(study$data[,caseControl==0]))
diag(coex) <- NA
# heatmap.2(coex, trace = "none", col="bluered", Colv = F, Rowv = F, dendrogram = "none", breaks=c(-.1,-.01,.01,.1))
hist(coex[row(coex)>col(coex)])

allCorrelations <- allCorrelations[sample(nrow(allCorrelations)),]

cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot(allCorrelations[allCorrelations$labels!="Background",]) + geom_point(aes(x=naiveMeth,y=newMeth, color=factor(labels)), alpha=.5, size=.4) +
    ggtitle("Comparison of our method vs naive approach") + ylab("Our Method") + xlab("Naive Approach") + theme_bw()+ scale_colour_manual(values=cbPalette)

ggplot(allCorrelations[1:20000,]) + geom_point(aes(x=naiveMeth,y=naiveWBatch, color=factor(labels)), alpha=1, size=.5) +
    ggtitle("Comparison of our method vs naive approach") + ylab("Naive with Batch") + xlab("Naive Approach")

ggplot(allCorrelations) + geom_density(aes(naiveMeth, color=factor(labels))) + theme_bw() +
    ggtitle("Estimated correlation difference using Naive Method") + scale_colour_manual(values=cbPalette)
ggplot(allCorrelations) + geom_density(aes(naiveWBatch, color=factor(labels))) + theme_bw() +
    ggtitle("Estimated correlation difference using Naive Method with Batch")
ggplot(allCorrelations) + geom_density(aes(newMeth, color=factor(labels))) + theme_bw() +
    ggtitle("Estimated correlation difference using Our Method")

onlyEffects <- allCorrelations[labels%in%c("Batch effect","Real effect","Real effect 2"),]

# Consider absolute correlations
onlyEffects[,1:3] <- abs(onlyEffects[,1:3])

methodPred  <- prediction(onlyEffects$newMeth, onlyEffects$labels%in%c("Real effect","Real effect 2"))
roc.methodPred  <- performance(methodPred, measure = c("tpr","auc"), x.measure = "fpr")
auc.methodPred  <- performance(methodPred, "auc")@y.values[[1]]

methodPredNaive  <- prediction(onlyEffects$naiveMeth, onlyEffects$labels%in%c("Real effect","Real effect 2"))
roc.methodPred.naive  <- performance(methodPredNaive, measure = c("tpr","auc"), x.measure = "fpr")
auc.methodPred.naive  <- performance(methodPredNaive, "auc")@y.values[[1]]

methodPredNaiveWBatch  <- prediction(onlyEffects$naiveWBatch, onlyEffects$labels%in%c("Real effect","Real effect 2"))
roc.methodPred.naive.w.batch  <- performance(methodPredNaiveWBatch, measure = c("tpr","auc"), x.measure = "fpr")
auc.methodPred.naive.w.batch  <- performance(methodPredNaiveWBatch, "auc")@y.values[[1]]

plot(roc.methodPred, main="Predicting True Covariance vs Batch Covariance", col = 2, lwd=3)
lines(roc.methodPred.naive@x.values[[1]], roc.methodPred.naive@y.values[[1]], col = 4, lwd=3)
lines(roc.methodPred.naive.w.batch@x.values[[1]], roc.methodPred.naive.w.batch@y.values[[1]], col = 5, lwd=3)

legend("bottomright", c(paste("Our Method",round(auc.methodPred,4)), 
                        paste("Naive",round(auc.methodPred.naive,4)), 
                        paste("Naive Batch",round(auc.methodPred.naive.w.batch,4))), 
       lty=1,lwd=1,col=c(2,4,5),title="Area under ROC curve")
abline(0,1)
# hist(differentialCorrelationNaive)