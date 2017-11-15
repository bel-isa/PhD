#10 Clusters
RVclu1 <- pam(rv_matrix1, k=10)

RVcolors <- sample(rainbow(10)) 
RVcolors <- RVcolors[as.vector(RVclu1$clustering)]

pcl <- as.vector(RVclu1$clustering)

pam1 <- cbind(rv_matrix1, pcl)

pam2 <- pam1[order(pam1[,34]),]

pam3 <- t(pam2[,-34])

pam4 <- cbind(pam3, pcl)

pam5 <- pam4[order(pam4[,34]),]

pam6 <- pam5[,-34]

RVcolors <- sample(rainbow(10)) 
RVcolors <- RVcolors[as.vector(pam2[,34])]

heatmap.2(pam6, col=color(), scale="none", Rowv=FALSE, Colv=FALSE, dendrogram="none", ColSideColors=RVcolors, RowSideColors=RVcolors, trace="none", density.info="none")

#9 Clusters
RVclu1 <- pam(rv_matrix1, k=9)

RVcolors <- sample(rainbow(9)) 
RVcolors <- RVcolors[as.vector(RVclu1$clustering)]

pcl <- as.vector(RVclu1$clustering)

pam1 <- cbind(rv_matrix1, pcl)

pam2 <- pam1[order(pam1[,34]),]

pam3 <- t(pam2[,-34])

pam4 <- cbind(pam3, pcl)

pam5 <- pam4[order(pam4[,34]),]

pam6 <- pam5[,-34]

RVcolors <- sample(rainbow(9)) 
RVcolors <- RVcolors[as.vector(pam2[,34])]

heatmap.2(pam6, col=color(), scale="none", Rowv=FALSE, Colv=FALSE, dendrogram="none", ColSideColors=RVcolors, RowSideColors=RVcolors, trace="none", density.info="none")

#8 Clusters
RVclu1 <- pam(rv_matrix1, k=8)

RVcolors <- sample(rainbow(8)) 
RVcolors <- RVcolors[as.vector(RVclu1$clustering)]

pcl <- as.vector(RVclu1$clustering)

pam1 <- cbind(rv_matrix1, pcl)

pam2 <- pam1[order(pam1[,34]),]

pam3 <- t(pam2[,-34])

pam4 <- cbind(pam3, pcl)

pam5 <- pam4[order(pam4[,34]),]

pam6 <- pam5[,-34]

RVcolors <- sample(rainbow(8)) 
RVcolors <- RVcolors[as.vector(pam2[,34])]

heatmap.2(pam6, col=color(), scale="none", Rowv=FALSE, Colv=FALSE, dendrogram="none", ColSideColors=RVcolors, RowSideColors=RVcolors, trace="none", density.info="none")

#7 Clusters
RVclu1 <- pam(rv_matrix1, k=7)

RVcolors <- sample(rainbow(7)) 
RVcolors <- RVcolors[as.vector(RVclu1$clustering)]

pcl <- as.vector(RVclu1$clustering)

pam1 <- cbind(rv_matrix1, pcl)

pam2 <- pam1[order(pam1[,34]),]

pam3 <- t(pam2[,-34])

pam4 <- cbind(pam3, pcl)

pam5 <- pam4[order(pam4[,34]),]

pam6 <- pam5[,-34]

RVcolors <- sample(rainbow(7)) 
RVcolors <- RVcolors[as.vector(pam2[,34])]

heatmap.2(pam6, col=color(), scale="none", Rowv=FALSE, Colv=FALSE, dendrogram="none", ColSideColors=RVcolors, RowSideColors=RVcolors, trace="none", density.info="none")

#6 Clusters
RVclu1 <- pam(rv_matrix1, k=6)

RVcolors <- sample(rainbow(6)) 
RVcolors <- RVcolors[as.vector(RVclu1$clustering)]

pcl <- as.vector(RVclu1$clustering)

pam1 <- cbind(rv_matrix1, pcl)

pam2 <- pam1[order(pam1[,34]),]

pam3 <- t(pam2[,-34])

pam4 <- cbind(pam3, pcl)

pam5 <- pam4[order(pam4[,34]),]

pam6 <- pam5[,-34]

RVcolors <- sample(rainbow(6)) 
RVcolors <- RVcolors[as.vector(pam2[,34])]

heatmap.2(pam6, col=color(), scale="none", Rowv=FALSE, Colv=FALSE, dendrogram="none", ColSideColors=RVcolors, RowSideColors=RVcolors, trace="none", density.info="none")

#5 Clusters
RVclu1 <- pam(rv_matrix1, k=5)

RVcolors <- sample(rainbow(5)) 
RVcolors <- RVcolors[as.vector(RVclu1$clustering)]

pcl <- as.vector(RVclu1$clustering)

pam1 <- cbind(rv_matrix1, pcl)

pam2 <- pam1[order(pam1[,34]),]

pam3 <- t(pam2[,-34])

pam4 <- cbind(pam3, pcl)

pam5 <- pam4[order(pam4[,34]),]

pam6 <- pam5[,-34]

RVcolors <- sample(rainbow(5)) 
RVcolors <- RVcolors[as.vector(pam2[,34])]

heatmap.2(pam6, col=color(), scale="none", Rowv=FALSE, Colv=FALSE, dendrogram="none", ColSideColors=RVcolors, RowSideColors=RVcolors, trace="none", density.info="none")

#4 Clusters
RVclu1 <- pam(rv_matrix1, k=4)

RVcolors <- sample(rainbow(4)) 
RVcolors <- RVcolors[as.vector(RVclu1$clustering)]

pcl <- as.vector(RVclu1$clustering)

pam1 <- cbind(rv_matrix1, pcl)

pam2 <- pam1[order(pam1[,34]),]

pam3 <- t(pam2[,-34])

pam4 <- cbind(pam3, pcl)

pam5 <- pam4[order(pam4[,34]),]

pam6 <- pam5[,-34]

RVcolors <- sample(rainbow(4)) 
RVcolors <- RVcolors[as.vector(pam2[,34])]

heatmap.2(pam6, col=color(), scale="none", Rowv=FALSE, Colv=FALSE, dendrogram="none", ColSideColors=RVcolors, RowSideColors=RVcolors, trace="none", density.info="none")
