
exprFile <- "~/gd/Harvard/Research/data/GTEx/GTEx_expr.txt"
clinicalFile <- "~/gd/Harvard/Research/data/GTEx/GTEx_clinical.txt"
#GTEX labels
casesString <- "cells_ebv-transformed_lymphocytes"
controlsString <- "skin"
phenotypeName <- "our_subtypes"

geneSubset <- sample(24369,2000)
gtexExp      <- read.table(exprFile,row.names=1,header=T)[geneSubset,]
gtexClinical <- read.table(clinicalFile,header=T,fill = TRUE, sep="\t",row.names=1)

bloodFilter <- gtexClinical$SMTSD=="Whole Blood"
celllinebloodFilter <- gtexClinical$SMTSD=="Cells - EBV-transformed lymphocytes"
lungFilter <- gtexClinical$SMTSD=="Lung"


centers <- as.character(gtexClinical$SMCENTER)
sum(centers=="C1")

gexp <- gtexExp[, bloodFilter|celllinebloodFilter]

batch <- c(rep(0,54), rep(1,191))
gexp[,batch==0] <- gexp[,batch==0]-rowMeans(gexp[,batch==0])
gexp[,batch==1] <- gexp[,batch==1]-rowMeans(gexp[,batch==1])



G <- cor(t(gexp))
G[is.na(G)] <- 0
n_star <- 120
m <- 2 # number of resamplings per s_i
s <- rep(10:44,m)

H <- lapply(s, function(s_i){
    sample_indicesA <- sample(54, s_i)
    sample_indicesB <- 54+sample(191, n_star-s_i)
    uncorrected <- cor(t(gexp[,c(sample_indicesA,sample_indicesB)]))
    uncorrected[is.na(uncorrected)] <- 0
    uncorrected
})

H_corrected <- lapply(s, function(s_i){
    sample_indicesA <- sample(54, s_i)
    sample_indicesB <- 54+sample(191, n_star-s_i)
    corrected <- (cor(t(gexp[,sample_indicesA])) + cor(t(gexp[,sample_indicesB])))/2
    corrected[is.na(corrected)] <- 0
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

png("./gtex_validation_cor.png", width=1200)
ggplot(res_table)  + ggtitle("GTEx - Correlation of networks with whole dataset by batch proportion") +
    geom_point(aes(x=batch_prop, y=validation_cor_uncorrected), alpha=.3, col="blue") + stat_smooth(aes(x=batch_prop, y=validation_cor_uncorrected), col="blue") +
    geom_point(aes(x=batch_prop, y=validation_cor_corrected), alpha=.3, col="red") + stat_smooth(aes(x=batch_prop, y=validation_cor_corrected), col="red")
dev.off()

png("./gtex_validation_aucroc.png", width=1200)
ggplot(res_table)  + ggtitle("GTEx - AUC-ROC of networks with whole dataset by batch proportion") +
    geom_point(aes(x=batch_prop, y=validation_aucroc_uncorrected), alpha=.3, col="blue") + stat_smooth(aes(x=batch_prop, y=validation_aucroc_uncorrected), col="blue") +
    geom_point(aes(x=batch_prop, y=validation_aucroc_corrected), alpha=.3, col="red") + stat_smooth(aes(x=batch_prop, y=validation_aucroc_corrected), col="red")
dev.off()

library(WGCNA)
gtexExp[] <- lapply(gtexExp, as.numeric)
fixedBloodData <- fixDataStructure(t(gtexExp[,bloodFilter]))
fixedLungData <- fixDataStructure(t(gtexExp[,lungFilter]))

bloodBCM <- blockwiseConsensusModules(fixedBloodData, power = 6, minModuleSize = 30, deepSplit = 2, maxBlockSize=30000,
                                      pamRespectsDendro = FALSE, 
                                      mergeCutHeight = 0.25, numericLabels = TRUE,
                                      minKMEtoStay = 0,
                                      saveTOMs = TRUE, verbose = 5)

lungBCM <- blockwiseConsensusModules(fixedLungData, power = 6, minModuleSize = 30, deepSplit = 2, maxBlockSize=30000,
                                     pamRespectsDendro = FALSE, 
                                     mergeCutHeight = 0.25, numericLabels = TRUE,
                                     minKMEtoStay = 0,
                                     saveTOMs = TRUE, verbose = 5)

save(bloodBCM, lungBCM, file = "blockwiseConsensusModules_blood_lung.RData")

load("blockwiseConsensusModules_blood_lung.RData")

# consMEs = bloodBCM$multiMEs;
moduleLabels = bloodBCM$colors[bloodBCM$goodGenes]
# Convert the numeric labels to color labels
moduleColorsBlood = labels2colors(moduleLabels)
consTree = bloodBCM$dendrograms[[1]];
png("./dendro_blood.png")
plotDendroAndColors(consTree, moduleColorsBlood,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors (GTEx Blood)")
dev.off()
# consMEs = lungBCM$multiMEs;
moduleLabels = lungBCM$colors[lungBCM$goodGenes]
# Convert the numeric labels to color labels
moduleColorsLung = labels2colors(moduleLabels)
consTree = lungBCM$dendrograms[[1]]; 
png("./dendro_lung.png")
plotDendroAndColors(consTree, moduleColorsLung,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors (GTEx Lung)")
dev.off()

library(igraph)
V(g)$x <- c(1, 1, 1, 2, 2, 2, 2)
V(g)$y <- c(3, 2, 1, 3.5, 2.5, 1.5, 0.5)
allColorsBlood <- labels2colors(bloodBCM$colors)
allColorsLung <- labels2colors(lungBCM$colors)
adjtable <- table(allColorsBlood, allColorsLung)
# adjtable[] <- as.numeric(scale(adjtable)>.2)
colnames(adjtable) <- paste0("LungModule",1:5)
rownames(adjtable) <- paste0("BloodModule",1:8)
# graph.bipartite(c(0,0,0,0,0,1,1,1,1,1,1,1,1), melt(adjtable))
df <- melt(adjtable)
png("./bipartiteBloodLung.png")
g <- graph_from_data_frame(df[df$value>8,], directed = FALSE)
plot(g)
dev.off()
png("./bipartiteBloodLung_CenterA.png")
plot(graph_from_data_frame(df[df$value>0,], directed = FALSE))
dev.off()
png("./bipartiteBloodLung_CenterB.png")
plot(graph_from_data_frame(df[df$value>3,], directed = FALSE))
dev.off()
png("./bipartiteBloodLung_CenterC.png")
plot(graph_from_data_frame(df[df$value>20,], directed = FALSE))
dev.off()
