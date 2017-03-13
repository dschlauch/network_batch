
simulateStudy <- function(numGenes, numSamples, addedError, blockSeq, mu, caseControl, batches){
    start <- Sys.time()
    blocks <- sort(rep(blockSeq, length.out=numGenes))
    blockA <- as.numeric(blocks=="A") # All Samples
    blockB <- as.numeric(blocks=="B") # Batch 1 only
    blockC <- as.numeric(blocks=="C") # Batch 2 only
    blockD <- as.numeric(blocks=="D") # Cases only
    blockE <- as.numeric(blocks=="E") # Controls only
    blockF <- as.numeric(blocks=="F") # Cases and negative with D 
    blockG <- as.numeric(blocks=="G") # No Samples
    
    batch1Effect   <- cbind(blockA, blockB)
    batch2Effect   <- cbind(blockA, blockC)
    casesEffect    <- cbind(blockD - blockF)
    controlsEffect <- cbind(blockE)
    
    SigmaBatch1 <- 4*tcrossprod(batch1Effect)
    SigmaBatch2 <- 4*tcrossprod(batch2Effect)
    SigmaCase <- tcrossprod(casesEffect)
    SigmaControl <- tcrossprod(controlsEffect)
    
    
    Sigmas <- list(Batch1Control=SigmaBatch1+SigmaControl, 
                   Batch1Case=SigmaBatch1+SigmaCase, 
                   Batch2Case=SigmaBatch2+SigmaCase, 
                   Batch2Control=SigmaBatch2+SigmaControl)
    Sigmas <- lapply(Sigmas, function(x){
        x<-x/addedError
        diag(x) <- 1 # adding a bit of random error
        x
    })
    counts <- list(Batch1Control=sum(batches==0&caseControl==0),
                   Batch1Case=sum(batches==0&caseControl==1),
                   Batch2Case=sum(batches==1&caseControl==1),
                   Batch2Control=sum(batches==1&caseControl==0))
    data <- cbind(t(mvrnorm(counts[['Batch1Control']],mu=mu, Sigma = Sigmas[['Batch1Control']])),
                  t(mvrnorm(counts[['Batch1Case']],mu=mu, Sigma = Sigmas[['Batch1Case']])),
                  t(mvrnorm(counts[['Batch2Case']],mu=mu, Sigma = Sigmas[['Batch2Case']])),
                  t(mvrnorm(counts[['Batch2Control']],mu=mu, Sigma = Sigmas[['Batch2Control']])))
    data[data<0] <-0
    
    print(paste("Data generation in",round(as.numeric(difftime(Sys.time(), start,units = "secs")),1), "seconds"))
    list(data=data, blocks=blocks, trueEffects=list(
        batch1Effect=batch1Effect,
        batch2Effect=batch2Effect,
        casesEffect=casesEffect,
        controlsEffect=controlsEffect))
    
}