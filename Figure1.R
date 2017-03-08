library(gplots)
library(RColorBrewer)
columns <- 10
rows <- 50
mu1 <- rnorm(rows,7,2)
mu2 <- rnorm(rows,7,2)
data <- cbind(t(sapply(mu1,function(mu) rnorm(columns,mu))),t(sapply(mu2,function(mu) rnorm(columns,mu))))

png("./figures/expression_heatmap.png", width=300, height=400)
heatmap.2(data, trace = "none", col = "bluered",key = F, dendro="none",labRow = FALSE, labCol = FALSE, Rowv = F)
dev.off()
system("convert ./figures/expression_heatmap.png -trim ./figures/expression_heatmap.png")

data.corrected <- cbind(data[,1:columns]-mu1, data[,(columns+1):(2*columns)]-mu2)

png("./figures/expression_heatmap_corrected.png", width=300, height=400)
heatmap.2(data.corrected, trace = "none", col = "bluered",key = F, dendro="none",labRow = FALSE, labCol = FALSE, Rowv = F)
dev.off()
system("convert ./figures/expression_heatmap_corrected.png -trim ./figures/expression_heatmap_corrected.png")


# vec <- rnorm(rows,mu1)-7
coexpression.mat <- cor(t(data[,1:10]))
# coexpression.mat <- coexpression.mat[rows:1,]

png("./figures/coexpression.png", width=400, height=400)
heatmap.2(t(coexpression.mat), trace = "none", col=brewer.pal(9,"Greens"),key = F, dendro="none",labRow = FALSE, labCol = FALSE)
dev.off()
system("convert ./figures/coexpression.png -trim ./figures/coexpression.png")
system("convert ./figures/coexpression.png -rotate 90 ./figures/coexpression.png")



# vec <- rnorm(rows,mu1)-7
coexpression.mat.corrected <- cor(t(data[,11:20]))
# coexpression.mat <- coexpression.mat[rows:1,]

png("./figures/coexpression.corrected.png", width=400, height=400)
heatmap.2(t(coexpression.mat.corrected), trace = "none", col=brewer.pal(9,"Greens"),key = F, dendro="none",labRow = FALSE, labCol = FALSE)
dev.off()
system("convert ./figures/coexpression.corrected.png -trim ./figures/coexpression.corrected.png")
system("convert ./figures/coexpression.corrected.png -rotate 90 ./figures/coexpression.corrected.png")
