#Additional computation following the random forest imputation

#Full Z transformed dataset stored in variable Z

# Cluster -----------------------------------------------------------------

library(pheatmap)

row.names(Z) <- Metabs[,1]

#Pearson correlation

pheatmap(Z, cluster_rows=TRUE, cluster_cols=TRUE, clustering_distance_rows="correlation", clustering_distance_cols="correlation", clustering_method="complete", legend=TRUE, show_rownames=T, show_colnames=T)

#Euclidean distance
pheatmap(Z, cluster_rows=TRUE, cluster_cols=TRUE, clustering_distance_rows="euclidean", clustering_distance_cols="euclidean", clustering_method="complete", legend=TRUE, show_rownames=T, show_colnames=T)