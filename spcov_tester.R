# 11/21/16
# This script will generate some batched data and run sparse covariance estimates on using spcov

library(spcov)
library(MASS)
library(gplots)
library(ggplot2)

sampleSizes <- rep(50,4)

numGenes <- 1000

# Sigma_Treat_A <- diag(numGenes)
# Sigma_Treat_B <- diag(numGenes)
# Sigma_Control_A <- diag(numGenes)
# Sigma_Control_B <- diag(numGenes)

# Use real data (ECLIPSE)

eclipseExp <- read.table("~/gd/Harvard/Research/data/Eclipse/ECLIPSE_Blood_Exp.txt",row.names=1,header=T)
set.seed(1000)
all <- sample(nrow(eclipseExp),numGenes+300)

Sigma_Treat_A   <- cor(t(eclipseExp[all[1:numGenes],]))
Sigma_Treat_B   <- cor(t(eclipseExp[all[c(1:(numGenes-100),numGenes+1:100)],]))
Sigma_Control_A <- cor(t(eclipseExp[all[c(1:(numGenes-100),numGenes+101:200)],]))
Sigma_Control_B <- cor(t(eclipseExp[all[c(1:(numGenes-100),numGenes+201:300)],]))
mu=rep(0,numGenes)

# generate the gene expression data
gexp_Treat_A <- mvrnorm(n=sampleSizes[1], mu=mu, Sigma = Sigma_Treat_A)
gexp_Treat_B <- mvrnorm(n=sampleSizes[2], mu=mu, Sigma = Sigma_Treat_B)
gexp_Control_A <- mvrnorm(n=sampleSizes[3], mu=mu, Sigma = Sigma_Control_A)
gexp_Control_B <- mvrnorm(n=sampleSizes[4], mu=mu, Sigma = Sigma_Control_B)


estimated_Treat_A <- cor(gexp_Treat_A)
estimated_Treat_B <- cor(gexp_Treat_B)
estimated_Control_A <- cor(gexp_Control_A)
estimated_Control_B <- cor(gexp_Control_B)
# Plot the true correlation matrices

diag(Sigma_Treat_A) <- NA
diag(Sigma_Treat_B) <- NA
diag(Sigma_Control_A) <- NA
diag(Sigma_Control_B) <- NA
rownames(Sigma_Treat_A) <- paste("Gene",1:numGenes)
rownames(Sigma_Treat_B) <- paste("Gene",1:numGenes)
rownames(Sigma_Control_A) <- paste("Gene",1:numGenes)
rownames(Sigma_Control_B) <- paste("Gene",1:numGenes)
colnames(Sigma_Treat_A) <- paste("Gene",1:numGenes)
colnames(Sigma_Treat_B) <- paste("Gene",1:numGenes)
colnames(Sigma_Control_A) <- paste("Gene",1:numGenes)
colnames(Sigma_Control_B) <- paste("Gene",1:numGenes)

combinedSigma <- abs(Sigma_Treat_A) + abs(Sigma_Treat_B) + abs(Sigma_Control_A) + abs(Sigma_Control_B)
hm2 <- heatmap.2(combinedSigma, col = "bluered", trace = "none")
geneOrder <- hm2$rowInd

heatmap.2(abs(Sigma_Treat_A[geneOrder,geneOrder]), col = "bluered", trace = "none", Rowv = FALSE, Colv = FALSE, dendrogram = "none")
heatmap.2(abs(Sigma_Treat_B[geneOrder,geneOrder]), col = "bluered", trace = "none", Rowv = FALSE, Colv = FALSE, dendrogram = "none")
heatmap.2(abs(Sigma_Control_A[geneOrder,geneOrder]), col = "bluered", trace = "none", Rowv = FALSE, Colv = FALSE, dendrogram = "none")
heatmap.2(abs(Sigma_Control_B[geneOrder,geneOrder]), col = "bluered", trace = "none", Rowv = FALSE, Colv = FALSE, dendrogram = "none")

heatmap.2(abs(Sigma_Treat_A-Sigma_Treat_B), col = "bluered", trace = "none", Rowv = FALSE, Colv = FALSE, dendrogram = "none")

gexp_Treat_A <- mvrnorm(sampleSizes[1], mu = rep(0,numGenes), Sigma = Sigma_Treat_A)
gexp_Treat_B <- mvrnorm(sampleSizes[1], mu = rep(0,numGenes), Sigma = Sigma_Treat_B)
gexp_Control_A <- mvrnorm(sampleSizes[1], mu = rep(0,numGenes), Sigma = Sigma_Control_A)
gexp_Control_B <- mvrnorm(sampleSizes[1], mu = rep(0,numGenes), Sigma = Sigma_Control_B)



spcov()