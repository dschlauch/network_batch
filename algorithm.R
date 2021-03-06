library(MASS)
library(gplots)
library(ggplot2)
library(rARPACK)

themethod <- function(X, expressionData, absolute=F, eigen_function=eigs_sym, N=ncol(expressionData), standardize=T, eigenG=NULL){
    start <- Sys.time()
    numSamples <- ncol(expressionData)
    if(absolute){
        corAdj <- abs
    } else {
        corAdj <- identity
    }
    N <- min(N,nrow(expressionData))
    
    if (standardize){
        G_star <- expressionData-rowMeans(expressionData)
        G <- (G_star/sqrt(rowSums(G_star^2)))
        G <- as.matrix(G)
    } else {
        G <- expressionData
        G <- (G/sqrt(rowSums(G^2)))
        G <- as.matrix(G)
    }
    if(is.null(eigenG)){
        eigenG <- eigen_function(corAdj(tcrossprod(G)),N)
    }
    Q <- eigenG$vectors
    D <- diag(eigenG$values)
    
    hatmat <- ginv(crossprod(X))%*%t(X)
    Qinv <- ginv(Q)
    
    #####
    QinvG <- Qinv%*%(G)
    
    est <- t(sapply(seq_len(nrow(hatmat)), function(hatmatRow){
        if(absolute){
            print(hatmatRow)
            middle <- Reduce("+",
                             lapply(seq_len(ncol(G)), function(i){abs(tcrossprod(G[,i]))*hatmat[hatmatRow,i]}))
            numSamples*diag(Qinv%*%middle%*%t(Qinv))
        } else {
            diag(QinvG%*%(numSamples*diag(hatmat[hatmatRow,]))%*%t(QinvG))
        }
    }))
    print(paste("Computation performed in",round(as.numeric(difftime(Sys.time(), start,units = "secs")),1), "seconds"))
    list(estimates=est, Q=Q, D=eigenG$values, G=G)
}

