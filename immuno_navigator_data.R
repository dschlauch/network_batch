## Yeast Data
library(data.table)
killerTCells <- read.table("~/Downloads/expression_22.txt",sep="\t", header = T, stringsAsFactors = F)
macrophages <- read.table("~/gd/Harvard/Research/data/ImmunoNavigator/expression_1.txt",sep="\t", header = T, stringsAsFactors = F)
allMouse <- read.table("~/Downloads/expression_0.txt",sep="\t", header = T, stringsAsFactors = F)
cellTypesMouse <- read.table("~/gd/Harvard/Research/data/ImmunoNavigator/sample_to_cell_types.txt",header = F, stringsAsFactors = F, sep="\t")

library(ggfortify)
autoplot(prcomp(t(killerTCells)))
covmat <- cov(t(killerTCells))

autoplot(princomp(covmat=cov(killerTCells)), type="lines")

as <- princomp()
biplot(as)


batchAssignmentsMouse <- read.table("./sample_to_batches.txt",row.names=1)

macrophageBatch <- batchAssignmentsMouse[colnames(macrophages),]
table(macrophageBatch)



intersect(batchAssignmentsMouse$V1,colnames(killerTCells))
colnames(killerTCells)%in%batchAssignmentsMouse$V1



batchAssignmentsMouse <- read.table("./sample_to_batches.txt")

## Macrophages
macrophages <- read.table("~/gd/Harvard/Research/data/ImmunoNavigator/expression_1.txt",sep="\t", header = T, stringsAsFactors = F)
macrophageBatch <- batchAssignmentsMouse[colnames(macrophages),]
includedBatches <- rownames(as.matrix(table(macrophageBatch)))[table(macrophageBatch)>10]

# Filter the NAs and below threshold
macrophages <- macrophages[,macrophageBatch%in%includedBatches]
macrophageBatch <- batchAssignmentsMouse[colnames(macrophages),]

# Make the design matrix
designMatrix <- cbind(1,do.call(cbind,lapply(unique(macrophageBatch),function(x){
    as.numeric(x==macrophageBatch)
})))

## Mature B Cells
matureBCells <- read.table("~/gd/Harvard/Research/data/ImmunoNavigator/expression_3.txt",sep="\t", header = T, stringsAsFactors = F)
matureBCellsBatch <- batchAssignmentsMouse[colnames(matureBCells),]
table(matureBCellsBatch)
includedBatches <- rownames(as.matrix(table(matureBCellsBatch)))[table(matureBCellsBatch)>10]

# Filter the NAs and below threshold
matureBCells <- matureBCells[,matureBCellsBatch%in%includedBatches]
matureBCellsBatch <- batchAssignmentsMouse[colnames(matureBCells),]

# Make the design matrix
designMatrix <- cbind(1,do.call(cbind,lapply(unique(matureBCellsBatch),function(x){
    as.numeric(x==matureBCellsBatch)
})))

## CD4Cells
CD4Cells <- read.table("~/gd/Harvard/Research/data/ImmunoNavigator/expression_2.txt",sep="\t", header = T, stringsAsFactors = F)
CD4CellsBatch <- batchAssignmentsMouse[colnames(CD4Cells),]
table(CD4CellsBatch)
includedBatches <- rownames(as.matrix(table(CD4CellsBatch)))[table(CD4CellsBatch)>10]

# Filter the NAs and below threshold
CD4Cells <- CD4Cells[,CD4CellsBatch%in%includedBatches]
CD4CellsBatch <- batchAssignmentsMouse[colnames(CD4Cells),]

# Make the design matrix
designMatrix <- cbind(1,do.call(cbind,lapply(unique(CD4CellsBatch),function(x){
    as.numeric(x==CD4CellsBatch)
})))


source('~/gd/Harvard/Research/network_batch/algorithm.R')

macrophagesResult <- themethod(designMatrix, macrophages[1:10000,])
