# LTCDNM Data
library(WGCNA)
source('~/gd/Harvard/Research/network_batch/algorithm.R')
ltcdnmExp <- read.table("~/gd/Harvard/Research/data/LTCOPD/LTCOPD_exp.txt",row.names=1,header=T)
rowSD <- apply(ltcdnmExp,1,sd)
ltcdnmExp <- ltcdnmExp[order(-rowSD)[1:5000],]

ltcdnmClin <- read.table("~/gd/Harvard/Research/data/LTCOPD/LTCOPD_clinical.txt",header=T,fill = TRUE, sep="\t",row.names=1)
# filterCovariates <- ltcdnmClin$COPD%in%c("case","cont")&ltcdnmClin$Gender%in%c("M","F")
# ltcdnmExp <- ltcdnmExp[,filterCovariates]
# ltcdnmClin <- ltcdnmClin[filterCovariates,]
designMatrix <- cbind(1,as.numeric(ltcdnmClin$COPD=="case"), as.numeric(ltcdnmClin$Gender=="M"))

ltcdnmResult <- themethod(designMatrix,ltcdnmExp, absolute = F, eigen_function = eigs_sym, N=100)


X_bar <- matrix(colMeans(designMatrix),nrow=1)
fitValues <- c(X_bar%*%ltcdnmResult$estimates)

X_cases <- c(1, 1, X_bar[3])
X_controls <- c(1, 0, X_bar[3])
X_Male <- c(1, X_bar[2], 1)
X_Non_male <- c(1, X_bar[2], 0)

correctedCasesCorrelation <- abs(ltcdnmResult$Q%*%diag(c(X_cases%*%ltcdnmResult$estimates))%*%t(ltcdnmResult$Q))
correctedControlsCorrelation <- abs(ltcdnmResult$Q%*%diag(c(X_controls%*%ltcdnmResult$estimates))%*%t(ltcdnmResult$Q))

# differentialCorrelationCOPD <- ltcdnmResult$Q%*%diag(ltcdnmResult$estimates[2,])%*%t(ltcdnmResult$Q)
# differentialCorrelationGender <- ltcdnmResult$Q%*%diag(ltcdnmResult$estimates[3,])%*%t(ltcdnmResult$Q)

correctedCasesCorrelation[correctedCasesCorrelation>1] <-1
correctedControlsCorrelation[correctedControlsCorrelation>1] <-1
TOMCases <- TOMsimilarity(correctedCasesCorrelation, TOMDenom = 'min', verbose = 1)
TOMControls <- TOMsimilarity(correctedControlsCorrelation, TOMDenom = 'min', verbose = 1)
multiExpr <- list(S1=list(data=t(ltcdnmExp[,ltcdnmClin$COPD=="case"])),S2=list(data=t(ltcdnmExp[,ltcdnmClin$COPD!="case"])))
consensusNetwork <- consensusDissTOMandTree(multiExpr, 
                                            6, TOM = list(cases=TOMCases,controls=TOMControls))

plotDendroAndColors(consensusNetwork$consTree, cbind(labels2colors(dat1$allLabels), 
                                                     labels2colors(dat2$allLabels)),c("cases","controls"), dendroLabels=FALSE)


# GOenrichmentAnalysis
modules = cutreeDynamic(dendro = consensusNetwork$consTree, distM = w, deepSplit = 4, pamRespectsDendro = FALSE,
                        minClusterSize = 30)
GOenrichmentAnalysis(moduleColors, allLLIDs, organism="human", nBestP=5)
