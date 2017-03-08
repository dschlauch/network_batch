library(MASS)
library(gplots)
library(ggplot2)


themethod <- function(X, expression, absolute=F){
    start <- Sys.time()
    numSamples <- ncol(expression)
    if(absolute){
        corAdj <- abs
    } else {
        corAdj <- identity
    }
    G_star <- expression-rowMeans(expression)
    G_standard <- (G_star/sqrt(rowSums(G_star^2)))
    G_standard <- as.matrix(G_standard)
    eigenG <- eigen(corAdj(tcrossprod(G_standard)))
    Q <- eigenG$vectors
    D <- diag(eigenG$values)
    
    hatmat <- ginv(crossprod(X))%*%t(X)
    Qinv <- ginv(Q)
    QinvG <- Qinv%*%(G_standard)
    
    est <- t(sapply(seq_len(nrow(hatmat)), function(hatmatRow){
        if(absolute){
            diag(Qinv%*%abs((G_standard)%*%(numSamples*diag(hatmat[hatmatRow,]))%*%t(G_standard))%*%t(Qinv))
        } else {
            diag(QinvG%*%(numSamples*diag(hatmat[hatmatRow,]))%*%t(QinvG))
        }
    }))
    print(paste("Computation performed in",round(as.numeric(difftime(Sys.time(), start,units = "secs")),1), "seconds"))
    list(estimates=est, Q=Q, D=eigenG$values, G_standard=G_standard)
}

