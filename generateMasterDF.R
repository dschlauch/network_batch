generateMasterDF <- function(ourMethodResult, correlationResult, correlationWBatchResult, trueLabels, maxPoints=1000000){
    
    X_bar <- matrix(colMeans(X),nrow=1)
    fitValues <- c(X_bar%*%ourMethodResult$estimates)
    correctedCorrelation <- ourMethodResult$Q%*%diag(fitValues)%*%t(ourMethodResult$Q)
    differentialCorrelation <- ourMethodResult$Q%*%diag(ourMethodResult$estimates[3,])%*%t(ourMethodResult$Q)
    
    

    
    allCorrelations <- data.frame(newMeth=differentialCorrelation[row(differentialCorrelation) > col(differentialCorrelation)],
                                  naiveMeth=correlationResult[row(correlationResult) > col(correlationResult)],
                                  naiveWBatch=correlationWBatchResult[row(correlationWBatchResult) > col(correlationWBatchResult)],
                                  labels=trueLabels)
    if(maxPoints){
        numberOfPointsToPlot <- min(maxPoints, nrow(allCorrelations))
        allCorrelations <- allCorrelations[sample(nrow(allCorrelations),numberOfPointsToPlot),]
    }
    allCorrelations
    # coex <- cor(t(study$data[,caseControl==0]))
    # diag(coex) <- NA
    # # heatmap.2(coex, trace = "none", col="bluered", Colv = F, Rowv = F, dendrogram = "none", breaks=c(-.1,-.01,.01,.1))
    # hist(coex[row(coex)>col(coex)])
}