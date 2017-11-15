setwd("Y:/Isabel Orf/R-3.1.2")
# Load data ---------------------------------------------------------------

Metabs <- read.delim("Y:/Isabel Orf/Projects/Ara_Syn/Opinion/15-02-23_vergleich-arasyn.txt")

View(Metabs)


# Z-transform data --------------------------------------------------------

Metabs_z <- scale(Metabs[,2:ncol(Metabs)])
Metabs_tz <- t(Metabs_z)
Metabs_z2 <- scale(Metabs_tz[,1:ncol(Metabs_tz)])
Metabs_z2 <- t(Metabs_z2)

rownames(Metabs_z2) <- Metabs[,1]
colnames(Metabs_tz) <- Metabs[,1]

View(Metabs_z)
View(Metabs_tz)
View(Metabs_z2)
# Correlation matrix ------------------------------------------------------

library(Hmisc)

Metabs_p <- rcorr(Metabs_z, type="pearson") #can also be spearman
Metabs_s <- rcorr(Metabs_z, type="spearman")

# To file -----------------------------------------------------------------

df.Metabs_p.r = data.frame(Metabs_p$r)
df.Metabs_p.n = data.frame(Metabs_p$n)
df.Metabs_p.p = data.frame(Metabs_p$p)

write.csv(df.Metabs_p.r,'correlationmatrix_p.csv')
write.csv(df.Metabs_p.p,'correlationmatrix_sign_p.csv')

df.Metabs_s.r = data.frame(Metabs_s$r)
df.Metabs_s.n = data.frame(Metabs_s$n)
df.Metabs_s.p = data.frame(Metabs_s$p)

write.csv(df.Metabs_s.r,'correlationmatrix_s.csv')
write.csv(df.Metabs_s.p,'correlationmatrix_sign_s.csv')

write.csv(Metabs_z,'Metabs_z.csv')


# Cluster -----------------------------------------------------------------

library(pheatmap)

#Pearson correlation
#pheatmap(Metabs_z2, cluster_rows=TRUE, cluster_cols=TRUE, clustering_method="complete", legend=TRUE, show_rownames=T, show_colnames=T)
#Stimmt offensichtlich nicht... euclidean distance scheint default zu sein.


#Euclidean distance
pheatmap(Metabs_z2, cluster_rows=TRUE, cluster_cols=TRUE, clustering_distance_rows="euclidean", clustering_distance_cols="euclidean", clustering_method="complete", legend=TRUE, show_rownames=T, show_colnames=T)

#Get rid of NAs
k <- Metabs_z2[complete.cases(Metabs_z2),]

#Get rid of Metabs not existing in sample groups
Metabs_red <- Metabs_z2[-11,]
Metabs_red1 <- Metabs_red[-15,]
Metabs_red2 <- Metabs_red1[-20,]
Metabs_red3 <- Metabs_red2[-23,]
Metabs_red4 <- Metabs_red3[-21,]
Metabs_red5 <- Metabs_red4[-10,]
Metabs_red6 <- Metabs_red5[-24,]

View(Metabs_red6)

#Get rid of samples with many missing values

Metabs_red7 <- Metabs_red6[,-23:-40]

View(Metabs_red7)

red <- Metabs_red7[-2,]
red1 <- red[-15,]
red2 <- red1[-8,]
red3 <- red2[-8,]
red4 <- red3[-6,]

View(red4)

#Only columns z transformed
rownames(Metabs_z) <- Metabs[,1]

#pheatmap(Metabs_z, cluster_rows=TRUE, cluster_cols=TRUE, clustering_method="complete", legend=TRUE, show_rownames=T, show_colnames=T)
pheatmap(Metabs_z, cluster_rows=TRUE, cluster_cols=TRUE, clustering_distance_rows="euclidean", clustering_distance_cols="euclidean", clustering_method="complete", legend=TRUE, show_rownames=T, show_colnames=T)