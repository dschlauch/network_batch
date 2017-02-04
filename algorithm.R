library(MASS)
library(gplots)
library(ggplot2)
set.seed(1234)
numGenes <- 2000 
numSamples <-400
mu <- rnorm(numGenes,mean = 9)
batches <- c(rep(0,numSamples/2),rep(1,numSamples/2))
caseControl <- c(rep(0,numSamples*3/8),rep(1,numSamples/2),rep(0,numSamples*1/8))
blockSeq <- c("A","B","B","B","C","C","C","D","D","D","D","E","E","F","F","G","G","G","G","G","G","G","G")
blocks <- rep(blockSeq, length.out=numGenes)
blockA <- as.numeric(blocks=="A")
blockB <- as.numeric(blocks=="B") # Batch 1 only
blockC <- as.numeric(blocks=="C") # Control only
blockD <- as.numeric(blocks=="D")
# blockE <- as.numeric(blocks=="E") # Batch 2 only
# blockF <- as.numeric(blocks=="F") # Cases only

batch1 <- cbind(blockA, blockB, blockD)
batch2 <- cbind(blockA, blockD)#, blockE)      
Sigma1 <- batch1%*%t(batch1)
Sigma2 <- batch2%*%t(batch2)
# SigmaCase <- blockF %*% t(blockF)
SigmaCase <- diag(numGenes)
SigmaControl <- blockC %*% t(blockC)


Sigmas <- list(case1=Sigma1+SigmaCase, 
               case2=Sigma2+SigmaCase, 
               control1=Sigma1+SigmaControl, 
               control2=Sigma2+SigmaControl)
Sigmas <- lapply(Sigmas, function(x){
    diag(x) <- 10 # adding a bit of random error
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
# heatmap.2(data, trace = "none", col = "bluered")

G_star <- data-rowMeans(data)
G_standard <- (G_star/sqrt(rowSums(G_star^2)))

eigenG <- eigen(G_standard%*%t(G_standard))
Q <- eigenG$vectors
D <- diag(eigenG$values)
plot(diag(D))
evsDF <- data.frame(Q)
names(evsDF) <- paste0("EV",1:ncol(evsDF))
evsDF$Gene <- 1:numGenes
evsDF$blocks <- blocks
ggplot(evsDF) + geom_point(aes(x=Gene,y=EV1, shape=blocks))
ggplot(evsDF) + geom_point(aes(x=Gene,y=EV2, shape=blocks))
ggplot(evsDF) + geom_point(aes(x=Gene,y=EV3, shape=blocks))
ggplot(evsDF) + geom_point(aes(x=Gene,y=EV4, shape=blocks))
ggplot(evsDF) + geom_point(aes(x=Gene,y=EV5, shape=blocks))

X <- cbind(rep(1,numSamples),batches, caseControl)
# X <- cbind(rep(1,numSamples))
# X <- cbind(rep(1,10),c(rep(1,5),rep(0,5)))
# X <- cbind(rep(1,10),c(rep(1,5),rep(0,5)), rnorm(10))


hatmat <- ginv(t(X)%*%X)%*%t(X)
est <- t(sapply(seq_len(nrow(hatmat)), function(hatmatRow){
    diag(ginv(Q)%*%(G_standard)%*%(numSamples*diag(hatmat[hatmatRow,]))%*%t(G_standard)%*%t(ginv(Q)))
}))

plot(est[1,1:20])
plot(est[2,1:20])
plot(est[3,1:20])

plot(eigenG$values[1:20])

X_bar <- matrix(colMeans(X),nrow=1)
fitValues <- c(X_bar%*%est)
correctedCorrelation <- Q%*%diag(fitValues)%*%t(Q)
plot(diag(D)[1:50])
plot(fitValues[1:50])
# heatmap.2(correctedCorrelation)

differentialCorrelation <- Q%*%diag(est[3,])%*%t(Q) 

differentialCorrelationNaive <- cor(t(G_standard[,caseControl==1]))-cor(t(G_standard[,caseControl==0]))


batchMat <- blockB %*% t(blockB)
realMat <- blockC %*% t(blockC)

labels <- rep("Background",choose(numGenes,2))
labels[batchMat[row(batchMat) > col(batchMat)]==1] <- "Batch effect"
labels[realMat[row(realMat) > col(realMat)]==1] <- "Real effect"

allCorrelations <- data.frame(newMeth=differentialCorrelation[row(differentialCorrelation) > col(differentialCorrelation)],
                              naiveMeth=differentialCorrelationNaive[row(differentialCorrelationNaive) > col(differentialCorrelationNaive)],
                              labels=labels)

ggplot(allCorrelations[1:20000,]) + geom_point(aes(x=naiveMeth,y=newMeth, color=factor(labels)), alpha=.4) +
    ggtitle("Comparison of our method vs naive approach") + ylab("Our Method") + xlab("Naive Approach")

ggplot(allCorrelations) + geom_density(aes(naiveMeth, color=factor(labels))) + theme_bw() +
    ggtitle("Estimated correlation difference using Naive Method")
ggplot(allCorrelations) + geom_density(aes(newMeth, color=factor(labels))) + theme_bw() +
    ggtitle("Estimated correlation difference using Our Method")


# hist(differentialCorrelationNaive)
