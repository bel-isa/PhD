#Additional computation following the random forest imputation

#Full Z transformed dataset stored in variable Z

library(cluster)
library(MASS)

row.names(Z) <- Metabs[,1]

T <- t(Z)

clu2 <- pam(T, k=2) 
#random number of clusters to get an overview
plot(clu2)


# Choose number of clusters -----------------------------------------------

for (i in 2:25){
  clu <- pam(T, k=i)
  cat(i," avg: ", round(clu$silinfo$avg.width,2))
  cat(" clus: ", round(clu$silinfo$clus.avg.widths,2),"\n")
}
#Output = silhouette values, the higher the average silhouette coefficient the better the choice

clu19 <- pam(T, k=19)
plot(clu19)
#Clusplot and Silhouette plot

eqscplot(T, pch=clu19$clustering)

