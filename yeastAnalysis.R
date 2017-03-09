## Yeast Data

yeastExp <- read.table("../data/CompileData/YeastData_NormalizedExpression.txt",sep="\t", header = T, stringsAsFactors = F)
rownames(yeastExp) <- make.unique(yeastExp$Name, sep = ".")
yeastExp <- yeastExp[,-1]
dim(yeastExp)

# Subset data
yeastExp <- yeastExp[c(T,F),]


batch1Coex <- cor(t(yeastExp[,c(T,F)]))
batch2Coex <- cor(t(yeastExp[,c(F,T)]))

diffCoex <- batch1Coex-batch2Coex

plot(as.numeric(yeastExp[6,]))
plot(as.numeric(yeastExp[6,rep(1:25,each=2) + rep(c(0,25),25)]))

designMatrix <- cbind(1,rep(1:25,each=2),rep(0:1,25))

naiveCor1 <- cor(t(yeastExp[c(T,F),1:20]))
naiveCor2 <- cor(t(yeastExp[c(T,F),31:50]))
source('~/gd/Harvard/Research/network_batch/algorithm.R')
yeastResult <- themethod(designMatrix,yeastExp[c(T,F),])

summary(rowSums(yeastExp))
naiveCorDif <- naiveCor1-naiveCor2
naiveCorDifVec <- naiveCorDif[row(naiveCorDif)>col(naiveCorDif)]

differentialCorrelation <- yeastResult$Q%*%diag(yeastResult$estimates[3,])%*%t(yeastResult$Q) 
differentialCorrelationVec <- differentialCorrelation[row(differentialCorrelation)>col(differentialCorrelation)]
length(differentialCorrelationVec)
plot(differentialCorrelationVec[1:10000],naiveCorDifVec[1:10000])
