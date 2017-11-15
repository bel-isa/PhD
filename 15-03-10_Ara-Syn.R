###############################################################
###                R Script 2015-03-10                      ###
###        Correlation analysis and Biclustering            ###
###                Comparison Ara-Syn                       ###
###                     Isabel Orf                          ###
###############################################################

# Set directory and load packages -----------------------------------------

setwd("Y:/Isabel Orf/R-3.1.2")
library(abind)
library(cluster)
library(MASS)
library(Hmisc)

# Load data ---------------------------------------------------------------

Metabs <- read.delim("Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/15-03-09_vergleich-arasyn.txt")

# Impute missing values using random forest -------------------------------

Jungle <- abind(replicate(100, missForest(Metabs[,-1])$ximp, simplify=FALSE), along=3)
### TIME CONSUMING, needs to be parallelized

Full <- apply(Jungle, c(1,2), mean)

write.table(Full,"Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-10_Results/2015-03-10_Full.txt", sep="\t")

# Scale data and check output ---------------------------------------------

Z <- scale(Full)
rownames(Z) <- Metabs[,1]

Mean <- colMeans(Z)

#Plot the column mean values, version 1
par(oma=c(1,1,1,1))
plot(Mean, main="Sample means after Z transformation", xlab="Sample", xaxt="n")
axis(1, seq(100, 700, by=100), las=2)

#Plot the column mean values, version 2
par(oma=c(2,1,1,1))
Names <- as.vector(colnames(Metabs[,-1]))
x <- seq(1,695,1)
plot(Mean, main="Sample means after Z transformation", xaxt="n", xlab="")
axis(side=1, at=x, labels= Names, las=2, cex.axis=0.1)

# Correlation matrix ------------------------------------------------------

Metabs_p <- rcorr(Z, type="pearson")
Metabs_s <- rcorr(Z, type="spearman")

#Store correlation coefficients and significance in text files
df.Metabs_p.r = data.frame(Metabs_p$r)
df.Metabs_p.n = data.frame(Metabs_p$n)
df.Metabs_p.p = data.frame(Metabs_p$P)

df.Metabs_s.r = data.frame(Metabs_s$r)
df.Metabs_s.n = data.frame(Metabs_s$n)
df.Metabs_s.p = data.frame(Metabs_s$P)

write.table(df.Metabs_p.r,"Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-10_Full_Z_pearson.txt", sep="\t")
write.table(df.Metabs_p.p,"Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-10_Full_Z_pearson_sign.txt", sep="\t")

write.table(df.Metabs_s.r,"Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-10_Full_Z_spearman.txt", sep="\t")
write.table(df.Metabs_s.p,"Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-10_Full_Z_spearman_sign.txt", sep="\t")

#Plot correlation matrices
pear <- Metabs_p$r
spear <- Metabs_s$r

myImagePlot <- function(x, ...){
  min <- min(x)
  max <- max(x)
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
      yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
      xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
      title <- Lst$title
    }
  }
  # check for null values
  if( is.null(xLabels) ){
    xLabels <- c(1:ncol(x))
  }
  if( is.null(yLabels) ){
    yLabels <- c(1:nrow(x))
  }
  
  layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
  
  # Red and green range from 0 to 1 while Blue ranges from 1 to 0
  ColorRamp <- rgb( seq(0,1,length=256),  # Red
                    seq(0,1,length=256),  # Green
                    seq(1,0,length=256))  # Blue
  ColorLevels <- seq(min, max, length=length(ColorRamp))
  
  # Reverse Y axis
  reverse <- nrow(x) : 1
  yLabels <- yLabels[reverse]
  x <- x[reverse,]
  
  # Data Map
  par(oma = c(10,10,1,1))
  image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
        ylab="", axes=FALSE, zlim=c(min,max))
  if( !is.null(title) ){
    title(main=title)
  }
  axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, las= 2, cex.axis=0.7)
  axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
       cex.axis=0.7) 
  
  # Color Scale
  par(mar = c(4,1,1,4))
  image(5, ColorLevels,
        matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
        col=ColorRamp,
        xlab="",ylab="",
        xaxt="n")
  
  layout(1)
}

myImagePlot(pear, title=c("Pearson correlation coefficient"))
myImagePlot(spear, title=c("Spearman rank coefficient"))

# PAM clustering ----------------------------------------------------------

Trans <- t(Z)

#Determine number of clusters for samples
for (i in 2:25){
  clu <- pam(Trans, k=i)
  cat(i," avg: ", round(clu$silinfo$avg.width,2))
  cat(" clus: ", round(clu$silinfo$clus.avg.widths,2),"\n")
}
###Select number of clusters with the highest average silhouette coefficient, here k=4

#Run pam for samples
clu4 <- pam(Trans, k=4)

#Plot PCA and silhouette plot
par(mar=c(5,1,3,1))
plot(clu4)
par(oma=c(1,3,1,1))
clusplot(clu4, main="Sample PCA [pam(Trans, k=4)]", color=FALSE, col.p="black", col.clus="black")

#Determine number of clusters for metabolites
for (i in 2:25){
  clu <- pam(Z, k=i)
  cat(i," avg: ", round(clu$silinfo$avg.width,2))
  cat(" clus: ", round(clu$silinfo$clus.avg.widths,2),"\n")
}
###Select number of clusters with the highest average silhouette coefficient, here k=2

#Run pam for metabolites
clu2 <- pam(Z, k=2)
clu3 <- pam(Z, k=3)

#Plot clustering for metabolites
eqscplot(Z, pch=clu2$clustering, main="Clustering metabolites [pam(Z, k=2)]")
eqscplot(Z, pch=clu3$clustering, main="Clustering metabolites [pam(Z, k=3)]")

# Create pearson correlation heatmap --------------------------------------

# Color function to generate green-red heat maps
color <- function(n = 50, low.col = 0.45, high.col=1, saturation = 1) { 
  if (n < 2) stop("n must be greater than 2")
  n1 <- n%/%2
  n2 <- n - n1
  c(hsv(low.col, saturation, seq(1,0,length=n1)), hsv(high.col, saturation, seq(0,1,length=n2))) 
}

#Cluster samples and metabolites
hmetab <- hclust(as.dist(1-cor(Trans, method="pearson")), method="complete")
hsample <- hclust(as.dist(1-cor(Z, method="pearson")), method="complete")

#Plot heatmap
par(oma=c(3,1,1,1))
heatmap(Z, Rowv=as.dendrogram(hmetab), Colv=as.dendrogram(hsample), col=color(), scale="none", main="Pearson correlation, complete clustering") 

#Cut tree according to pam clustering results
pamsample <- sample(rainbow(4)) 
pamsample <- pamsample[as.vector(clu4$clustering)]
pammetab <- sample(rainbow(4))
pammetab <- pammetab[as.vector(clu3$clustering)]
par(oma=c(5,2,2,5))
heatmap(Trans, Rowv=as.dendrogram(hsample), Colv=as.dendrogram(hmetab), col=color(), scale="none", RowSideColors=pamsample, ColSideColors=pammetab, main="Biclustering (pearson, complete, PAM ky=4 kx=3)")

# Count exp or cond Ara/Syn per cluster -----------------------------------

exp <- c(1,  1,1,  1,  1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	1,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	2,	3,	3,	3,	3,	3,	3,	3,	3,	3,	3,	3,	3,	3,	3,	3,	3,	3,	3,	3,	3,	3,	3,	3,	3,	3,	4,	4,	4,	4,	4,	4,	4,	4,	4,	4,	4,	4,	4,	4,	4,	4,	4,	4,	5,	5,	5,	5,	5,	5,	5,	5,	5,	5,	5,	5,	5,	5,	5,	5,	5,	5,	5,	5,	5,	5,	5,	5,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	6,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	7,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	8,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9,	9)

cond <- c(1,    1,  1,	1,	2,	2,	2,	2,	3,	3,	3,	3,	4,	4,	4,	4,	5,	5,	5,	6,	6,	6,	1,	1,	1,	1,	1,	1,	2,	2,	2,	2,	2,	2,	3,	3,	3,	3,	3,	3,	1,	1,	1,	1,	1,	1,	1,	1,	1,	2,	2,	2,	2,	2,	2,	2,	3,	3,	3,	3,	3,	3,	3,	3,	3,	1,	1,	1,	2,	2,	2,	3,	3,	3,	4,	4,	4,	5,	5,	5,	6,	6,	6,	7,	7,	7,	7,	7,	7,	8,	8,	8,	8,	8,	8,	9,	9,	9,	9,	9,	9,	10,	10,	10,	10,	10,	10,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	11,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	12,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	13,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	14,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	15,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	16,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	17,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	18,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	19,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	20,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	21,	22,	22,	22,	22,	22,	22,	22,	22,	22,	22,	22,	22	,22	,22	,22	,22	,22	,22	,22	,22	,22	,22	,22	,22	,22	,22	,22	,22	,22	,22	,22	,22	,22	,22	,22	,22	,22	,22	,22,	22,	22,	22)

numm <- cbind(exp, cond)

c1 <- data.frame(pam1 = pam(Trans, k=1)$clustering)
c1a <- cbind(c1, numm)
write.table(c1a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_C1.txt", sep="\t")

c2 <- data.frame(pam2 = pam(Trans, k=2)$clustering)
c2a <- cbind(c2, numm)
write.table(c2a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_C2.txt", sep="\t")

c3 <- data.frame(pam3 = pam(Trans, k=3)$clustering)
c3a <- cbind(c3, numm)
write.table(c3a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_C3.txt", sep="\t")

c4 <- data.frame(pam4 = pam(Trans, k=4)$clustering)
c4a <- cbind(c4, numm)
write.table(c4a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_C4.txt", sep="\t")

c5 <- data.frame(pam5 = pam(Trans, k=5)$clustering)
c5a <- cbind(c5, numm)
write.table(c5a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_C5.txt", sep="\t")

c6 <- data.frame(pam6 = pam(Trans, k=6)$clustering)
c6a <- cbind(c6, numm)
write.table(c6a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_C6.txt", sep="\t")

c7 <- data.frame(pam7 = pam(Trans, k=7)$clustering)
c7a <- cbind(c7, numm)
write.table(c7a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_C7.txt", sep="\t")

c8 <- data.frame(pam8 = pam(Trans, k=8)$clustering)
c8a <- cbind(c8, numm)
write.table(c8a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_C8.txt", sep="\t")

c9 <- data.frame(pam9 = pam(Trans, k=9)$clustering)
c9a <- cbind(c9, numm)
write.table(c9a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_C9.txt", sep="\t")

c10 <- data.frame(pam10 = pam(Trans, k=10)$clustering)
c10a <- cbind(c10, numm)
write.table(c10a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_C10.txt", sep="\t")

c11 <- data.frame(pam11 = pam(Trans, k=11)$clustering)
c11a <- cbind(c11, numm)
write.table(c11a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_C11.txt", sep="\t")

c12 <- data.frame(pam12 = pam(Trans, k=12)$clustering)
c12a <- cbind(c12, numm)
write.table(c12a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_C12.txt", sep="\t")

c13 <- data.frame(pam13 = pam(Trans, k=13)$clustering)
c13a <- cbind(c13, numm)
write.table(c13a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_C13.txt", sep="\t")

c14 <- data.frame(pam14 = pam(Trans, k=14)$clustering)
c14a <- cbind(c14, numm)
write.table(c14a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_C14.txt", sep="\t")

c15 <- data.frame(pam15 = pam(Trans, k=15)$clustering)
c15a <- cbind(c15, numm)
write.table(c15a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_C15.txt", sep="\t")

c16 <- data.frame(pam16 = pam(Trans, k=16)$clustering)
c16a <- cbind(c16, numm)
write.table(c16a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_C16.txt", sep="\t")

c17 <- data.frame(pam17 = pam(Trans, k=17)$clustering)
c17a <- cbind(c17, numm)
write.table(c17a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_C17.txt", sep="\t")

c18 <- data.frame(pam18 = pam(Trans, k=18)$clustering)
c18a <- cbind(c18, numm)
write.table(c18a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_C18.txt", sep="\t")

c19 <- data.frame(pam19 = pam(Trans, k=19)$clustering)
c19a <- cbind(c19, numm)
write.table(c19a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_C19.txt", sep="\t")

c20 <- data.frame(pam20 = pam(Trans, k=20)$clustering)
c20a <- cbind(c20, numm)
write.table(c20a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_C20.txt", sep="\t")

c21 <- data.frame(pam21 = pam(Trans, k=21)$clustering)
c21a <- cbind(c21, numm)
write.table(c21a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_C21.txt", sep="\t")

c22 <- data.frame(pam22 = pam(Trans, k=22)$clustering)
c22a <- cbind(c22, numm)
write.table(c22a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_C22.txt", sep="\t")

c23 <- data.frame(pam23 = pam(Trans, k=23)$clustering)
c23a <- cbind(c23, numm)
write.table(c23a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_C23.txt", sep="\t")

c24 <- data.frame(pam24 = pam(Trans, k=24)$clustering)
c24a <- cbind(c24, numm)
write.table(c24a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_C24.txt", sep="\t")

c25 <- data.frame(pam25 = pam(Trans, k=25)$clustering)
c25a <- cbind(c25, numm)
write.table(c25a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_C25.txt", sep="\t")

c26 <- data.frame(pam26 = pam(Trans, k=26)$clustering)
c26a <- cbind(c26, numm)
write.table(c26a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_C26.txt", sep="\t")

c27 <- data.frame(pam27 = pam(Trans, k=27)$clustering)
c27a <- cbind(c27, numm)
write.table(c27a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_C27.txt", sep="\t")

c28 <- data.frame(pam28 = pam(Trans, k=28)$clustering)
c28a <- cbind(c28, numm)
write.table(c28a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_C28.txt", sep="\t")

c29 <- data.frame(pam29 = pam(Trans, k=29)$clustering)
c29a <- cbind(c29, numm)
write.table(c29a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_C29.txt", sep="\t")

c30 <- data.frame(pam30 = pam(Trans, k=30)$clustering)
c30a <- cbind(c30, numm)
write.table(c30a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_C30.txt", sep="\t")

write.table(table(c1a$pam, c1a$exp), "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_exp_groups.txt", sep="\t")
write.table(table(c2a$pam, c2a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c3a$pam, c3a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c4a$pam, c4a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c5a$pam, c5a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c6a$pam, c6a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c7a$pam, c7a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c8a$pam, c8a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c9a$pam, c9a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c10a$pam, c10a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c11a$pam, c11a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c12a$pam, c12a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c13a$pam, c13a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c14a$pam, c14a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c15a$pam, c15a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c16a$pam, c16a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c17a$pam, c17a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c18a$pam, c18a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c19a$pam, c19a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c20a$pam, c20a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c21a$pam, c21a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c22a$pam, c22a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c23a$pam, c23a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c24a$pam, c24a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c25a$pam, c25a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c26a$pam, c26a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c27a$pam, c27a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c28a$pam, c28a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c29a$pam, c29a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c30a$pam, c30a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_exp_groups.txt", sep="\t", append=TRUE)

write.table(table(c1a$pam, c1a$cond), "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_cond_groups.txt", sep="\t")
write.table(table(c2a$pam, c2a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c3a$pam, c3a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c4a$pam, c4a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c5a$pam, c5a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c6a$pam, c6a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c7a$pam, c7a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c8a$pam, c8a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c9a$pam, c9a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c10a$pam, c10a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c11a$pam, c11a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c12a$pam, c12a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c13a$pam, c13a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c14a$pam, c14a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c15a$pam, c15a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c16a$pam, c16a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c17a$pam, c17a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c18a$pam, c18a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c19a$pam, c19a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c20a$pam, c20a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c21a$pam, c21a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c22a$pam, c22a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c23a$pam, c23a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c24a$pam, c24a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c25a$pam, c25a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c26a$pam, c26a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c27a$pam, c27a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c28a$pam, c28a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c29a$pam, c29a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c30a$pam, c30a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-16_Results/2015-03-16_Count_cond_groups.txt", sep="\t", append=TRUE)

#Plot graph

nclus <- c(1:30)

cond.ara.syn.nclus <- c(2.666666667,  3.75,	1.5,	2.5,	1.8,	1.333333333,	1,	0.875,	0.333333333,	0.3,	0.272727273,	0.25,	0.230769231,	0.142857143,	0.133333333,	0.125,	0,	0,	0,	0,	0.047619048,	0.045454545,	0.043478261,	0,	0,	0,	0,	0,	0,	0)

cond.syn.ara.nclus <- c(0.375,  0.416666667,	0.75,	0.15625,	0.225,	0.19047619,	0.020408163,	0.017857143,	0.037037037,	0.033333333,	0.03030303,	0.027777778,	0.025641026,	0.035714286,	0.033333333,	0.03125,	0,	0,	0,	0,	0.047619048,	0.045454545,	0.043478261,	0,	0,	0,	0,	0,	0,	0)

plot(nclus,cond.ara.syn.nclus, type="o", main="Number of experimental conditions Ara/Syn per cluster", xlab="Number of clusters", ylab="Experimental conditions Ara/Syn per cluster")
plot(nclus,cond.syn.ara.nclus, type="o", main="Number of experimental conditions Syn/Ara per cluster", xlab="Number of clusters", ylab="Experimental conditions Syn/Ara per cluster")

exp.ara.syn.nclus <- c(1.25,  0.25,	0.083333333,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0)

exp.syn.ara.nclus <- c(0.8,  1,	1.333333333,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0,	0)

plot(nclus,exp.ara.syn.nclus, type="o", main="Number of experiments Ara/Syn per cluster", xlab="Number of clusters", ylab="Experiments Ara/Syn per cluster")
plot(nclus,exp.syn.ara.nclus, type="o", main="Number of experiments Syn/Ara per cluster", xlab="Number of clusters", ylab="Experiments Syn/Ara per cluster")

# Remove all samples that belong to a specific cluster --------------------

#Here cluster 2
Trans.clu4 <- cbind(Trans, clu4$clustering)
x <- 2
Reduced <- Trans.clu4[Trans.clu4[,ncol(Trans.clu4)]!=x,]

#Remove column with cluster specifier
ncol(Reduced)
Reduced1 <- Reduced[,-31]

#Determine number of clusters for samples
for (i in 2:25){
  clu <- pam(Reduced1, k=i)
  cat(i," avg: ", round(clu$silinfo$avg.width,2))
  cat(" clus: ", round(clu$silinfo$clus.avg.widths,2),"\n")
}
###Select number of clusters with the highest average silhouette coefficient, here k=3

#Run pam for samples
Red.clu3 <- pam(Reduced1, k=3)

#Plot PCA and silhouette plot
par(mar=c(5,2,3,2))
plot(Red.clu3)
par(oma=c(1,3,1,1))
clusplot(Red.clu3, main="Sample PCA [pam(Reduced, k=3)]", color=FALSE, col.p="black", col.clus="black")

#Determine number of clusters for metabolites
Red.T <- t(Reduced1)

for (i in 2:25){
  clu <- pam(Red.T, k=i)
  cat(i," avg: ", round(clu$silinfo$avg.width,2))
  cat(" clus: ", round(clu$silinfo$clus.avg.widths,2),"\n")
}
###Select number of clusters with the highest average silhouette coefficient, here k=5 (no negative values)

#Run pam for metabolites
Red.T.clu5 <- pam(Red.T, k=5)

#Cluster samples and metabolites
h.Red <- hclust(as.dist(1-cor(Reduced1, method="pearson")), method="complete")
h.Red.T <- hclust(as.dist(1-cor(Red.T, method="pearson")), method="complete")

#Cut tree according to pam clustering results
pamsample <- sample(rainbow(3)) 
pamsample <- pamsample[as.vector(Red.clu3$clustering)]
pammetab <- sample(rainbow(6))
pammetab <- pammetab[as.vector(Red.T.clu5$clustering)]
par(oma=c(5,2,2,5))
heatmap(Reduced1, Rowv=as.dendrogram(h.Red.T), Colv=as.dendrogram(h.Red), col=color(), scale="none", RowSideColors=pamsample, ColSideColors=pammetab, main="Biclustering reduced dataset (pearson, complete, PAM ky=3 kx=5)")

#Read-out pam clusters (samples)
Red.pam <- data.frame(sort(Red.clu3$clustering))
write.table(Red.pam, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-17_Results/2015-03-17_Pam_groups_reduced.txt", sep="\t")

# Remove Ara temp experiment ----------------------------------------------

#Starts with row number 529
NoTemp <- Trans[-(529:nrow(Trans)),]

#Determine number of clusters for samples
for (i in 2:25){
  clu <- pam(NoTemp, k=i)
  cat(i," avg: ", round(clu$silinfo$avg.width,2))
  cat(" clus: ", round(clu$silinfo$clus.avg.widths,2),"\n")
}
###Select number of clusters with the highest average silhouette coefficient, here k=4

#Run pam for samples
NoT.clu4 <- pam(NoTemp, k=4)

#Plot PCA and silhouette plot
par(mar=c(5,2,3,2))
plot(NoT.clu4)
par(oma=c(1,3,1,1))
clusplot(NoT.clu4, main="Sample PCA [pam(NoTemp, k=4)]", color=FALSE, col.p="black", col.clus="black")

#Determine number of clusters for metabolites
NoT.T <- t(NoTemp)

for (i in 2:25){
  clu <- pam(NoT.T, k=i)
  cat(i," avg: ", round(clu$silinfo$avg.width,2))
  cat(" clus: ", round(clu$silinfo$clus.avg.widths,2),"\n")
}
###Select number of clusters with the highest average silhouette coefficient, here k=3

#Run pam for metabolites
NoT.T.clu3 <- pam(NoT.T, k=3)

#Cluster samples and metabolites
h.NoT <- hclust(as.dist(1-cor(NoTemp, method="pearson")), method="complete")
h.NoT.T <- hclust(as.dist(1-cor(NoT.T, method="pearson")), method="complete")

#Cut tree according to pam clustering results
pamsample <- sample(rainbow(4)) 
pamsample <- pamsample[as.vector(NoT.clu4$clustering)]
pammetab <- sample(rainbow(3))
pammetab <- pammetab[as.vector(NoT.T.clu3$clustering)]
par(oma=c(5,2,2,5))
heatmap(NoTemp, Rowv=as.dendrogram(h.NoT.T), Colv=as.dendrogram(h.NoT), col=color(), scale="none", RowSideColors=pamsample, ColSideColors=pammetab, main="Biclustering dataset without Ara temperature data (pearson, complete, PAM ky=4 kx=2)")

#Read-out pam clusters (samples)
NoT.pam <- data.frame(sort(NoT.clu4$clustering))
write.table(NoT.pam, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-17_Results/2015-03-17_Pam_groups_NoTemp.txt", sep="\t")
