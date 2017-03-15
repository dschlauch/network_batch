# WGCNA
library(WGCNA)
library(cluster)

getClustersFromWGCNA <- function(expressionData){
    k=softConnectivity(datE=t(expressionData),power=6) 
    # sizeGrWindow(10,5)
    # par(mfrow=c(1,2))
    # hist(k)
    # scaleFreePlot(k, main="Check scale free topology\n")
    
    datExpr=t(expressionData)[, rank(-k,ties.method="first" )<=3600]
    
    ADJ1 <- abs(cor(t(expressionData))^6)
    dissADJ=1-ADJ1
    
    dissTOM=TOMdist(ADJ1)
    collectGarbage()
    # pam6=pam(as.dist(dissADJ), 6)
    # # Cross-tabulte the detected and the true (simulated) module membership:
    # table(pam6$clustering)

    pamTOM8=pam(as.dist(dissTOM), 8)
    pamTOM8$clustering
    # hierADJ=hclust(as.dist(dissADJ), method="average" )
    # colorStaticADJ=as.character(cutreeStaticColor(hierADJ, cutHeight=.99, minSize=100))
    # colorStaticADJ
}
