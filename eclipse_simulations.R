library(MASS)
library(ggplot2)
library(data.table)

numGenes <- 1000
numSamples <- 200
gexp <- matrix(rnorm(numGenes*numSamples,10), numGenes, numSamples)


eclipseExp <- read.table("~/gd/Harvard/Research/data/Eclipse/ECLIPSE_Blood_Exp.txt",row.names=1,header=T)

SigmaA <- cor(t(eclipseExp[1:numGenes,]))
muA <- rowMeans(eclipseExp[1:numGenes,])
SigmaB <- cor(t(eclipseExp[(numGenes+1):(2*numGenes),]))
muB <- rowMeans(eclipseExp[(numGenes+1):(2*numGenes),])
genesA <- t(mvrnorm(n = numSamples/2, muA, SigmaA, tol = 1e-6, empirical = FALSE, EISPACK = FALSE))
genesB <- t(mvrnorm(n = numSamples/2, muB, SigmaB, tol = 1e-6, empirical = FALSE, EISPACK = FALSE))
gexp <- cbind(genesA, genesB)



batch <- c(rep(0,numSamples/2), rep(1,numSamples/2))
gexp[,batch==0] <- gexp[,batch==0]-rowMeans(gexp[,batch==0])
gexp[,batch==1] <- gexp[,batch==1]-rowMeans(gexp[,batch==1])


G <- cor(t(gexp))
n_star <- 100
m <- 2 # number of resamplings per s_i
s <- rep(10:90,m)

H <- lapply(s, function(s_i){
    sample_indicesA <- sample(100, s_i)
    sample_indicesB <- 100+sample(100, n_star-s_i)
    uncorrected <- cor(t(gexp[,c(sample_indicesA,sample_indicesB)]))
    uncorrected
})

H_corrected <- lapply(s, function(s_i){
    sample_indicesA <- sample(100, s_i)
    sample_indicesB <- 100+sample(100, n_star-s_i)
    corrected <- (cor(t(gexp[,sample_indicesA])) + cor(t(gexp[,sample_indicesB])))/2
    corrected
})


benchmark <- G[row(G)>col(G)]

validate_nets_cor <- function(H){
    unlist(lapply(H, function(H_i){
        cor(benchmark, H_i[row(H_i)>col(H_i)])
    }))
}

require(ROCR)
validate_nets_aucroc <- function(H){
    unlist(lapply(H, function(H_i){
        methodPred  <- prediction(H_i[row(H_i)>col(H_i)],benchmark>.1)
        roc.methodPred  <- performance(methodPred, measure = c("tpr","auc"), x.measure = "fpr")
        auc.methodPred  <- performance(methodPred, "auc")@y.values[[1]]
    }))
}

validation_cor_uncorrected <- validate_nets_cor(H)
validation_cor_corrected <- validate_nets_cor(H_corrected)
validation_aucroc_uncorrected <- validate_nets_aucroc(H)
validation_aucroc_corrected <- validate_nets_aucroc(H_corrected)

res_table <- data.table(validation_cor_uncorrected, validation_cor_corrected, validation_aucroc_uncorrected, validation_aucroc_corrected, batch_prop=s/n_star)

png("./validation_cor.png", width=1200)
ggplot(res_table)  + ggtitle("Correlation of networks with whole dataset by batch proportion") +
    geom_point(aes(x=batch_prop, y=validation_cor_uncorrected), alpha=.3, col="blue") + stat_smooth(aes(x=batch_prop, y=validation_cor_uncorrected), col="blue") +
    geom_point(aes(x=batch_prop, y=validation_cor_corrected), alpha=.3, col="red") + stat_smooth(aes(x=batch_prop, y=validation_cor_corrected), col="red")
dev.off()

png("./validation_aucroc.png", width=1200)
ggplot(res_table)  + ggtitle("AUC-ROC of networks with whole dataset by batch proportion") +
    geom_point(aes(x=batch_prop, y=validation_aucroc_uncorrected), alpha=.3, col="blue") + stat_smooth(aes(x=batch_prop, y=validation_aucroc_uncorrected), col="blue") +
    geom_point(aes(x=batch_prop, y=validation_aucroc_corrected), alpha=.3, col="red") + stat_smooth(aes(x=batch_prop, y=validation_aucroc_corrected), col="red")
dev.off()





xBatch <- cbind(1,gexp[2,],gexp[2,]*batch)
resBatch <- ginv(crossprod(xBatch))%*%t(xBatch)%*%t(gexp)

x <- cbind(1,gexp[2,])
res <- ginv(crossprod(x))%*%t(x)%*%t(gexp)

resultDT <- data.table(cbind(t(resBatch)[,2:3],t(res)[,2]))
names(resultDT) <- c("Overall", "batch", "naive")

newDT <- data.table(c(resultDT$Overall,resultDT$batch,resultDT$naive),rep

png("./singleGeneCorrection.png", width=1200)
ggplot(resultDT) +geom_point(aes(y=Overall,x=1:1000))+geom_point(aes(y=batch,x=1:1000),col="red") + geom_point(aes(y=naive,x=1:1000),col="blue") +
    geom_segment(x=1:1000,xend=1:1000,aes(y=Overall,yend=batch), alpha=.3) + ggtitle("Correlations")
dev.off()

newDT <- data.table(c(resultDT$Overall,resultDT$batch,resultDT$naive),rep(c("Overall","batchDif","naive"),each=length(resultDT$Overall)))
names(newDT) <- c("correlation","type")
png("./singleGeneCorrection.png", width=1200)
ggplot(newDT) +geom_point(aes(y=correlation,x=rep(1:1000,3),color=type))+
    geom_segment(x=1:1000,xend=1:1000,aes(y=Overall,yend=batch), alpha=.3)
dev.off()