#Fig4 Pearson correlation
pear <- hclust(as.dist(1-cor(rv_matrix1, method="pearson")), method="average")

par(oma = c(5,1,1,5))
heatmap.2(rv_matrix1, col=color(), scale="none", Rowv=as.dendrogram(pear), Colv="Rowv", trace="none", density.info="none")
