library(MASS)
library(gplots)
library(ROCR)
library(ggplot2)
library(rARPACK)
library(grid)
library(gridExtra)
library(reshape2)

# Creates all the diagnosis plots for a result object
diagnosticPlots <- function(differentialCorrelationsDF, dir="."){
    dir.create(dir)
    
    
    cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    
    plotOursVsNaive <- ggplot(differentialCorrelationsDF[differentialCorrelationsDF$labels!="Background",]) + 
        geom_point(aes(x=naiveMeth,y=newMeth, color=factor(labels)), alpha=.5, size=2) +
        ggtitle("Pairwise Differential Coexpression estimates") + ylab("Our Method") + xlab("Naive Approach") + theme_bw(base_size = 30)+
        guides(color=guide_legend(title="True Interaction")) + scale_colour_manual(values=cbPalette)
    
    png(paste0(dir,'/OursVsNaive.png'), width = 1600, height = 1200)
    print(plotOursVsNaive)
    dev.off()
    
    naiveDensity <- ggplot(differentialCorrelationsDF) + 
        geom_density(aes(naiveMeth, color=factor(labels), fill=factor(labels)), alpha=0.3) + 
        theme_bw(base_size = 30) +
        ggtitle("Estimated Pairwise Differential Coexpression using Naive Method") + xlab("Pairwise Differential Coexpression") +
        guides(color=guide_legend(title="True Interaction")) + scale_colour_manual(values=cbPalette)
    naiveWBatchDensity <- ggplot(differentialCorrelationsDF) + 
        geom_density(aes(naiveWBatch, color=factor(labels), fill=factor(labels)), alpha=0.3) + 
        theme_bw(base_size = 30) +
        ggtitle("Estimated Pairwise Differential Coexpression using Naive Method with Batch") + xlab("Pairwise Differential Coexpression") +
        guides(color=guide_legend(title="True Interaction")) + scale_colour_manual(values=cbPalette)
    ourMethodDensity <- ggplot(differentialCorrelationsDF) + 
        geom_density(aes(newMeth, color=factor(labels), fill=factor(labels)), alpha=0.3) + 
        theme_bw(base_size = 30) +
        ggtitle("Estimated Pairwise Differential Coexpression using Our Method") + xlab("Pairwise Differential Coexpression") +
        guides(color=guide_legend(title="True Interaction")) + scale_colour_manual(values=cbPalette)
    
    png(paste0(dir,'/NaiveDensity.png'), width = 1600, height = 1200)
    print(naiveDensity)
    dev.off()
    
    png(paste0(dir,'/NaiveWBatchDensity.png'), width = 1600, height = 1200)
    print(naiveWBatchDensity)
    dev.off()
    
    png(paste0(dir,'/ourMethodDensity.png'), width = 1600, height = 1200)
    print(ourMethodDensity)
    dev.off()
    
    onlyEffects <- differentialCorrelationsDF[differentialCorrelationsDF$labels!="Background",]
    
    # Consider absolute correlations
    # onlyEffects[,1:3] <- abs(onlyEffects[,1:3])
    ########
    
    plotROC <- function(corrDF, positive, plottitle="Title"){
        library(ROCR)
        corrDF$newMeth <- abs(corrDF$newMeth)
        corrDF$naiveMeth <- abs(corrDF$naiveMeth)
        corrDF$naiveWBatch <- abs(corrDF$naiveWBatch)
        
        methodPred  <- prediction(corrDF$newMeth, corrDF$labels%in%positive)
        roc.methodPred  <- performance(methodPred, measure = c("tpr","auc"), x.measure = "fpr")
        auc.methodPred  <- performance(methodPred, "auc")@y.values[[1]]
        
        methodPredNaive  <- prediction(corrDF$naiveMeth, corrDF$labels%in%positive)
        roc.methodPred.naive  <- performance(methodPredNaive, measure = c("tpr","auc"), x.measure = "fpr")
        auc.methodPred.naive  <- performance(methodPredNaive, "auc")@y.values[[1]]
        
        methodPredNaiveWBatch  <- prediction(corrDF$naiveWBatch, corrDF$labels%in%positive)
        roc.methodPred.naive.w.batch  <- performance(methodPredNaiveWBatch, measure = c("tpr","auc"), x.measure = "fpr")
        auc.methodPred.naive.w.batch  <- performance(methodPredNaiveWBatch, "auc")@y.values[[1]]
        
        plot(roc.methodPred, main=plottitle, col = 2, lwd=3)
        lines(roc.methodPred.naive@x.values[[1]], roc.methodPred.naive@y.values[[1]], col = 4, lwd=3)
        lines(roc.methodPred.naive.w.batch@x.values[[1]], roc.methodPred.naive.w.batch@y.values[[1]], col = 5, lwd=3)
        
        legend("bottomright", c(paste("Our Method",round(auc.methodPred,4)), paste("Naive",round(auc.methodPred.naive,4)), paste("Naive Batch",round(auc.methodPred.naive.w.batch,4))), 
               lty=1,lwd=1,col=c(2,4,5),title="Area under ROC curve")
        abline(0,1)
    }
    unique(differentialCorrelationsDF$labels)
    png(paste0(dir,'/OursVsOthersROC.png'), width = 800, height = 400)
    par(mfrow=c(1,2))
    plotROC(differentialCorrelationsDF, positive=c("Real effect","Negative effect"), "Real Effects vs All Pairs")
    plotROC(onlyEffects, positive=c("Real effect","Negative effect"), "Real Effects vs Batch Effects")
    dev.off()
    par(mfrow=c(1,1))
    
}

plotEigenvectors <- function(ourMethodResult, trueEffects, dir=".", numEigenvectors=6){
    cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    
    dir.create(dir)
    numGenes <- nrow(ourMethodResult$G_standard)
    geneLabels <- rep("Background", numGenes)
    
    geneLabels[(rowSums(trueEffects$batch1Effect)-rowSums(trueEffects$batch2Effect))!=0] <- "Batch"
    geneLabels[(rowSums(trueEffects$casesEffect)-rowSums(trueEffects$controlsEffect))!=0] <- "Real"
    nrows<-nrow(ourMethodResult$estimates)
    
    plottingDF <- data.frame(ourMethodResult$Q[,seq_len(numEigenvectors)], geneLabels)
    names(plottingDF)[seq_len(numEigenvectors)] <- paste("Eigenvector",seq_len(numEigenvectors))
    plottingDFMelt <- melt(plottingDF, id.vars="geneLabels")
    plottingDFMelt$geneLabels <- geneLabels
    plottingDFMelt$Gene <- seq_len(numGenes)
    
    eigenvectorPlots <- ggplot(plottingDFMelt) + geom_point(aes(x=Gene,y=value,color=geneLabels)) +
        facet_wrap(~variable) + theme_bw(base_size = 30) +
        guides(color=guide_legend(title="Gene Group")) + scale_colour_manual(values=cbPalette)
    
    numVectors <- min(ncol(ourMethodResult$estimates),20)
    est <- t(ourMethodResult$estimates[,1:20])
    colnames(est)[1:3] <- c("Intercept","Batch","Case-Control") 
    eigenvalueDF <- melt(est,value.name = "Eigenvalue")
    
    eigenvaluePlots <- ggplot(eigenvalueDF) + geom_point(aes(x=Var1,y=Eigenvalue)) + 
        facet_wrap(~Var2) + theme_bw(base_size = 30) +
        ylab("Eigenvalue Coefficent") + xlab("Top 20 Eigenvectors") + scale_colour_manual(values=cbPalette)
    
    png(paste0(dir,'/EigenvectorPlots.png'), width = 1600, height = 1200)
    multiplot(eigenvectorPlots,eigenvaluePlots)
    dev.off()
}




multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
        print(plots[[1]])
        
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}