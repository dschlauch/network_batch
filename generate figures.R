#  Plots for paper

# Demo of simple case, four plots in one
library(sva)
library(MASS)
library(ggplot2)
library(limma)
library(gridExtra)

numSamples <- 1000
gene1 <- c(rnorm( numSamples,7),rnorm(numSamples,9))
gene2 <- c(rnorm(numSamples,7),rnorm(numSamples,9))
batch <- c(rep("A",numSamples),rep("B",numSamples))

df <- data.frame(gene1,gene2,batch)

plotA <- ggplot(df, aes(x=gene1,y=gene2)) +geom_point(aes(col=batch),size=3, alpha=.5) +ggtitle("Raw data") +theme_bw() 
correct <- c(rep(0,numSamples),rep(2,numSamples))
bc <- data.frame(gene1=gene1-correct,gene2=gene2-correct,batch)
plotB <- ggplot(bc, aes(x=gene1,y=gene2)) +geom_point(aes(col=batch),size=3, alpha=.5) +ggtitle("Batch corrected data") +theme_bw() 

genes <- mvrnorm(n = numSamples, c(7,7), matrix(c(1,.95,.95,1),nrow=2), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
gene1 <- c(genes[,1],rnorm(numSamples,9))
gene2 <- c(genes[,2],rnorm(numSamples,9))
df2 <- data.frame(gene1=gene1,gene2=gene2,batch)
plotC <- ggplot(df2, aes(x=gene1,y=gene2)) +geom_point(aes(col=batch),size=3, alpha=.5) +ggtitle("Raw data") +theme_bw() 
df2 <- data.frame(gene1=gene1-correct,gene2=gene2-correct,batch)
plotD <- ggplot(df2, aes(x=gene1,y=gene2)) +geom_point(aes(col=batch),size=3, alpha=.5) +ggtitle("Batch corrected data") +theme_bw() 


png("~/gd/Harvard/Research/R_workspace/batch_manuscript/figures/simulated_example.png", width=1200, height=1200)
grid.arrange(plotA, plotB, plotC, plotD, ncol=2)
dev.off()
