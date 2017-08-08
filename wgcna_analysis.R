library(GOstats)
library(org.Hs.eg.db)
library(WGCNA)
library(AnnotationDbi)


source('~/gd/Harvard/Research/network_batch/algorithm.R')
# Load data
numGenes <- 12000
exprFile <- "~/gd/Harvard/Research/data/COPDGene/COPDGene_GSExpressionData.txt"
clinicalFile <- "~/gd/Harvard/Research/data/COPDGene/COPDGene_clinical.txt"
clinicalFile <- "~/gd/Harvard/Research/data/COPDGene/GSE42057_series_matrix_clinical.txt"
clinicalData <- t(read.table(clinicalFile,header=F,fill = TRUE, sep="\t",stringsAsFactors = F))
colnames(clinicalData) <- clinicalData[1,]
clinicalData <- clinicalData[-1,]
colnames(clinicalData) <- gsub("!","",colnames(clinicalData))

dataset <- list()
dataset$exp      <- read.table(exprFile,row.names=1,header=T)
colnames(dataset$exp) <- substring(colnames(dataset$exp),1,10)
dataset$exp <- dataset$exp[,clinicalData[,"Sample_geo_accession"]]

dataset$clinical <- data.frame(clinicalData)
copdStatus <- as.numeric(grepl("COPD subject",dataset$clinical$Sample_source_name_ch1))
gender <- as.numeric(gsub("gender: ","",as.character(dataset$clinical$gender)))
packyears <- as.numeric(gsub("ats_packyears: ","",as.character(dataset$clinical$ats_packyears)))
age <- as.numeric(gsub("age_enroll: ","",as.character(dataset$clinical$age_enroll)))

# Subset data
# k.data=softConnectivity(t(dataset$exp),power=6)-1 
# kRank = rank(-k.data)
# topConnectedGenes <- dataset$exp[kRank<=numGenes,]
topConnectedGenes <- dataset$exp[sample(nrow(dataset$exp),numGenes),]

# Compute Differential coexpression
# diffco <- cor(t(topConnectedGenes[,variable1]))-cor(t(topConnectedGenes[,!variable1]))
# ADJdataDiffCo <- diffco^2
# ADJdataDiffCo[ADJdataDiffCo>1] <- 1
# Compute DiffCoex

# Compute Our differential coexpression
design <- cbind(1,copdStatus, gender, packyears, age)
result <- themethod(design, as.matrix(topConnectedGenes), absolute=F, eigen_function = eigs_sym)

ourMethoddifferentialCorrelation <- result$Q%*%diag(result$estimates[2,])%*%t(result$Q)
ADJdataOurs <- abs(ourMethoddifferentialCorrelation)

# Run WGCNA, get N modules for each dc method
# dissTOMdiffco=TOMdist(ADJdataDiffCo) 
dissTOMOurs=TOMdist(ADJdataOurs) 
# hierTOMdiffco = hclust(as.dist(dissTOMdiffco),method="average")
hierTOMoursCC = hclust(as.dist(dissTOMOurs),method="average")

par(mfrow=c(2,1),mar=c(2,4,1,1))
plot(hierTOMoursCC, main="COPDGene Case-Control", labels=F, xlab="", sub="");
colorhdataOneCC= cutreeStaticColor(hierTOMoursCC,cutHeight = 0.85, minSize = 10) 
plotColorUnderTree(hierTOMoursCC,colors=data.frame(module=colorhdataOneCC))
title("Module membership data set I") 
# Run TopGO on modules
colors <- unique(colorhdataOneCC)

module1Genes <- rownames(topConnectedGenes)[colorhdataOneCC==("turquoise")]
module2Genes <- rownames(topConnectedGenes)[colorhdataOneCC==("blue")]
module3Genes <- rownames(topConnectedGenes)[colorhdataOneCC==("brown")]

write.table(rownames(topConnectedGenes), file="./allGeneList.txt", quote = F, row.names = F, col.names = F)
write.table(module1Genes, file="./module1OursGeneList_casecontrol.txt", quote = F, row.names = F, col.names = F)
write.table(module2Genes, file="./module2OursGeneList_casecontrol.txt", quote = F, row.names = F, col.names = F)
write.table(module3Genes, file="./module3OursGeneList_casecontrol.txt", quote = F, row.names = F, col.names = F)


# Gender results
ourMethoddifferentialCorrelation <- result$Q%*%diag(result$estimates[3,])%*%t(result$Q)
ADJdataOurs <- abs(ourMethoddifferentialCorrelation)

dissTOMOurs=TOMdist(ADJdataOurs) 
hierTOMours = hclust(as.dist(dissTOMOurs),method="average")

# plot(hierTOMours,labels=F,main="Dendrogram, ours, 5000 most connected") 

par(mfrow=c(2,1),mar=c(2,4,1,1))
plot(hierTOMours, main="COPDGene Gender", labels=F, xlab="", sub="");
colorhdataOneGender= cutreeStaticColor(hierTOMours,cutHeight = 0.92, minSize = 10) 
plotColorUnderTree(hierTOMours,colors=data.frame(module=colorhdataOneGender))
title("Module membership data set I") 

module1GenesGender <- rownames(topConnectedGenes)[colorhdataOneGender==("turquoise")]
module2GenesGender <- rownames(topConnectedGenes)[colorhdataOneGender==("blue")]

write.table(module1GenesGender, file="./module1OursGeneList_Gender.txt", quote = F, row.names = F, col.names = F)
write.table(module2GenesGender, file="./module2OursGeneList_Gender.txt", quote = F, row.names = F, col.names = F)

# Packyears results
ourMethoddifferentialCorrelation <- result$Q%*%diag(result$estimates[4,])%*%t(result$Q)
ADJdataOurs <- abs(ourMethoddifferentialCorrelation)

dissTOMOurs=TOMdist(ADJdataOurs) 
hierTOMours = hclust(as.dist(dissTOMOurs),method="average")

plot(hierTOMours,labels=F,main="Dendrogram, ours, 5000 most connected") 

par(mfrow=c(2,1),mar=c(2,4,1,1))
plot(hierTOMours, main="COPDGene Pack-years", labels=F, xlab="", sub="");
colorhdataOnePY= cutreeStaticColor(hierTOMours,cutHeight = 0.992, minSize = 10) 
plotColorUnderTree(hierTOMours,colors=data.frame(module=colorhdataOnePY))
title("Module membership data set I") 

module1GenesPY <- rownames(topConnectedGenes)[colorhdataOnePY==("turquoise")]
module2GenesPY <- rownames(topConnectedGenes)[colorhdataOnePY==("blue")]

write.table(module1GenesPY, file="./module1OursGeneList_Packyears.txt", quote = F, row.names = F, col.names = F)
write.table(module2GenesPY, file="./module2OursGeneList_Packyears.txt", quote = F, row.names = F, col.names = F)


# Age results
ourMethoddifferentialCorrelation <- result$Q%*%diag(result$estimates[5,])%*%t(result$Q)
ADJdataOurs <- abs(ourMethoddifferentialCorrelation)

dissTOMOurs=TOMdist(ADJdataOurs) 
hierTOMours = hclust(as.dist(dissTOMOurs),method="average")

plot(hierTOMours,labels=F,main="Dendrogram, ours, 5000 most connected") 

par(mfrow=c(2,1),mar=c(2,4,1,1))
plot(hierTOMours, main="COPDGene Age", labels=F, xlab="", sub="");
colorhdataOneAge= cutreeStaticColor(hierTOMours,cutHeight = 0.99, minSize = 10) 
plotColorUnderTree(hierTOMours,colors=data.frame(module=colorhdataOneAge))
title("Module membership data set I") 

module1GenesAge <- rownames(topConnectedGenes)[colorhdataOneAge==("turquoise")]
module2GenesAge <- rownames(topConnectedGenes)[colorhdataOneAge==("blue")]

write.table(module1GenesAge, file="./module1OursGeneList_Age.txt", quote = F, row.names = F, col.names = F)
write.table(module2GenesAge, file="./module2OursGeneList_Age.txt", quote = F, row.names = F, col.names = F)


# Standard diffco results
diffco <- cor(t(topConnectedGenes[,copdStatus==1]))-cor(t(topConnectedGenes[,copdStatus==0]))
ADJdataDiffCo <- diffco^2
ADJdataDiffCo[ADJdataDiffCo>1]<-1
dissTOMOurs=TOMdist(ADJdataDiffCo) 
hierTOMours = hclust(as.dist(dissTOMOurs),method="average")

plot(hierTOMours,labels=F,main="Dendrogram, ours, 5000 most connected") 

par(mfrow=c(2,1),mar=c(2,4,1,1))
plot(hierTOMours, main="COPDGene Case-control (diffcoex)", labels=F, xlab="", sub="");
colorhdataOneDiffco= cutreeStaticColor(hierTOMours,cutHeight = 0.9, minSize = 10) 
plotColorUnderTree(hierTOMours,colors=data.frame(module=colorhdataOneDiffco))
title("Module membership data set I") 

module1GenesDiffCo <- rownames(topConnectedGenes)[colorhdataOneDiffco==("turquoise")]
module2GenesDiffCo <- rownames(topConnectedGenes)[colorhdataOneDiffco==("blue")]

write.table(module1Genes, file="./module1OursGeneList_Standard.txt", quote = F, row.names = F, col.names = F)
write.table(module2Genes, file="./module2OursGeneList_Standard.txt", quote = F, row.names = F, col.names = F)


# Plot all top modules under Case-control tree

pdf("./figures/module_membership.pdf")
par(mfrow=c(6,1),mar=c(1,4,1,1))
layout(c(1,1,1,1,2,3,4,5,6))
plot(hierTOMoursCC, main="COPDGene Case-Control", labels=F, xlab="", sub="")
plotBar <- function(x, label="module"){
    xModule1 <- ifelse(x=="turquoise","red","black")
    plotColorUnderTree(hierTOMoursCC,colors=data.frame(module=xModule1), label)
}
plotBar(colorhdataOneCC,"Case-control")
plotBar(colorhdataOneGender, "Gender")
plotBar(colorhdataOneAge, "Age")
plotBar(colorhdataOnePY, "Pack-years")
plotBar(colorhdataOneDiffco, "Diffco")
# plotBar(colorhdataOneRandom, "Random")
dev.off()

# GO analysis
goAnalysis <- function(module,filename){
    genemap <- AnnotationDbi::select(org.Hs.eg.db, module, "ENTREZID", "SYMBOL")
    genemap <- genemap[!duplicated(genemap$ENTREZID),]
    univmap <- AnnotationDbi::select(org.Hs.eg.db, rownames(topConnectedGenes), "ENTREZID", "SYMBOL")
    univmap <- univmap[!duplicated(univmap$ENTREZID),]
    
    params <- new("GOHyperGParams",
        geneIds=genemap$ENTREZID,
        universeGeneIds=univmap$ENTREZID,
        annotation="org.Hs.eg.db",
        ontology="BP",
        pvalueCutoff=1,
        conditional=TRUE,
        testDirection="over")
    
    hyperGRes <- hyperGTest(params)
    summary(hyperGRes)[1:20,]
    summary(hyperGRes)[grepl("anatomical",summary(hyperGRes)$Term),]
    summary(hyperGRes)[grepl("localization",summary(hyperGRes)$Term),]
    summary(hyperGRes)[grepl("morph",summary(hyperGRes)$Term),]
    saveRDS(hyperGRes,file = filename)
    hyperGRes
}
hyperGRes <- goAnalysis(module1Genes, "GO_casecontrol.rds")
goAnalysis(module1GenesGender, "GO_gender.rds")
goAnalysis(module1GenesAge, "GO_age.rds")
goAnalysis(module1GenesPY, "GO_packyears.rds")
goAnalysis(module1GenesDiffCo, "GO_diffco.rds")

# write.table(genemap$ENTREZID, file="./module1ENTREZ.txt", quote = F, row.names = F, col.names = F)
# write.table(univmap$ENTREZID, file="./allENTREZ.txt", quote = F, row.names = F, col.names = F)

randomGenes <- sample(rownames(topConnectedGenes),1000)
write.table(randomGenes, file="./randomModule.txt", quote = F, row.names = F, col.names = F)


plot(result$estimates[2,1:20])
points(result$estimates[1,1:20], col=2)
points(result$estimates[3,1:20], col=3)
points(result$estimates[4,1:20], col=4)
points(result$estimates[5,1:20], col=5)

