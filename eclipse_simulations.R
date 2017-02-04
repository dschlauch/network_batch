library(MASS)
library(ggplot2)
library(data.table)

numGenes <- 1000
numSamples <- 200
gexp <- matrix(rnorm(numGenes*numSamples,10), numGenes, numSamples)


eclipseExp <- read.table("~/gd/Harvard/Research/data/Eclipse/ECLIPSE_Blood_Exp.txt",row.names=1,header=T)
eclipseClinical <- read.table("~/gd/Harvard/Research/data/Eclipse/ECLIPSE_blood.txt",header=T,fill = TRUE, sep="\t",row.names=1)

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

# newDT <- data.table(c(resultDT$Overall,resultDT$batch,resultDT$naive),rep

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



library(WGCNA)

cases <- eclipseClinical$Subject.type=="COPD" 
controls <- eclipseClinical$Subject.type=="Smoker Control" 

males <- eclipseClinical$Gender=="M" 

# dataA <- eclipseExp[10101:12100,c(T,F)]
# dataB <- eclipseExp[10101:12100,c(F,T,F,F)]

dataA <- eclipseExp[101:2100, cases&(!males)]
dataB <- eclipseExp[101:2100, controls&(males)]

fixedDataA <- fixDataStructure(t(dataA))
fixedDataB <- fixDataStructure(t(dataB))


fixedDataA_BCM <- blockwiseConsensusModules(fixedDataA, power = 6, minModuleSize = 30, deepSplit = 2, maxBlockSize=30000,
                                      pamRespectsDendro = FALSE, 
                                      mergeCutHeight = 0.25, numericLabels = TRUE,
                                      minKMEtoStay = 0,
                                      saveTOMs = TRUE, verbose = 5)
# adjMat <- adjacency(t(dataA), power = 6)
# TOM <- TOMsimilarity(adjMat, TOMDenom = 'min', verbose = 1)
# consensusNetwork <- consensusDissTOMandTree(fixedDataA, 6, TOM = TOM)

fixedDataB_BCM <- blockwiseConsensusModules(fixedDataB, power = 6, minModuleSize = 30, deepSplit = 2, maxBlockSize=30000,
                                     pamRespectsDendro = FALSE, 
                                     mergeCutHeight = 0.25, numericLabels = TRUE,
                                     minKMEtoStay = 0,
                                     saveTOMs = TRUE, verbose = 5)
# summary(lm(fixedDataA_BCM$colors==0 ~ as.factor(fixedDataB_BCM$colors)))

# vglmFitMN <- vglm(as.factor(lungBCM$colors) ~ sample(as.factor(bloodBCM$colors)), family=multinomial(refLevel=1))

# pseudo-R^2
pseudoMultiR <- function(groupsA, groupsB, method="nagelkerke"){
    require(VGAM)
    vglmFitMN <- vglm(as.factor(groupsA) ~ as.factor(groupsB), family=multinomial(refLevel=1))
    vglm0 <- vglm(as.factor(groupsA) ~ 1, family=multinomial(refLevel=1))
    LLf   <- VGAM::logLik(vglmFitMN)
    LL0   <- VGAM::logLik(vglm0)
    N <- length(groupsA)
    if(method=="mcfadden"){
        as.vector(1 - (LLf / LL0))
    } else if (method=="coxsnell"){
        as.vector(1 - exp((2/N) * (LL0 - LLf)))
    } else if (method=="nagelkerke"){
        as.vector((1 - exp((2/N) * (LL0 - LLf))) / (1 - exp(LL0)^(2/N)))
    }
}

pseudoMultiR(fixedDataA_BCM$colors, fixedDataB_BCM$colors)

femalesOnly 
malesOnly
MvsF
FvsM


propMale <- (25:75)/100
pseudoR2 <-propMale*(1-propMale)*2+rnorm(51)/10
data <- data.frame(propMale=propMale,pseudoR2=pseudoR2)
ggplot(data) + geom_point(aes(x=propMale, y=pseudoR2)) + stat_smooth(aes(x=propMale, y=pseudoR2), col="red") + ggtitle("Agreement vs Confounder Balance")



propMale <- (25:75)/100
pseudoR2Corrected <-propMale*(1-propMale)*2+rnorm(51)/7
dataCorrected <- data.frame(propMale=propMale,pseudoR2Corrected=pseudoR2Corrected)
ggplot(dataCorrected) + geom_point(aes(x=propMale, y=pseudoR2Corrected)) + stat_smooth(aes(x=propMale, y=pseudoR2Corrected), col="blue") + ggtitle("Agreement (corrected) vs Confounder Balance")
