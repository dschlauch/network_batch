library(MASS)
library(gplots)
library(ggplot2)
library(rARPACK)

themethod <- function(X, expressionData, absolute=F, eigen_function=eigen, N=ncol(expressionData)){
    start <- Sys.time()
    numSamples <- ncol(expressionData)
    if(absolute){
        corAdj <- abs
    } else {
        corAdj <- identity
    }
    N <- min(N,nrow(expressionData))
    G_star <- expressionData-rowMeans(expressionData)
    G_standard <- (G_star/sqrt(rowSums(G_star^2)))
    G_standard <- as.matrix(G_standard)
    eigenG <- eigen_function(corAdj(tcrossprod(G_standard)),N)
    Q <- eigenG$vectors
    D <- diag(eigenG$values)
    
    hatmat <- ginv(crossprod(X))%*%t(X)
    Qinv <- ginv(Q)
    QinvG <- Qinv%*%(G_standard)
    
    est <- t(sapply(seq_len(nrow(hatmat)), function(hatmatRow){
        if(absolute){
            print(hatmatRow)
            middle <- Reduce("+",
                             lapply(seq_len(ncol(G_standard)), function(i){abs(tcrossprod(G_standard[,i]))*hatmat[hatmatRow,i]}))
            numSamples*diag(Qinv%*%middle%*%t(Qinv))
        } else {
            diag(QinvG%*%(numSamples*diag(hatmat[hatmatRow,]))%*%t(QinvG))
        }
    }))
    print(paste("Computation performed in",round(as.numeric(difftime(Sys.time(), start,units = "secs")),1), "seconds"))
    list(estimates=est, Q=Q, D=eigenG$values, G_standard=G_standard)
}

