library(abind)
library(cluster)
library(MASS)
library(Hmisc)

# Load data ---------------------------------------------------------------

Metabs <- read.delim("Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/15-03-09_vergleich-arasyn.txt")

# Impute missing values using random forest -------------------------------

Jungle <- abind(replicate(100, missForest(Metabs[,-1])$ximp, simplify=FALSE), along=3)
### TIME CONSUMING, needs to be parallelized!!!

Full <- apply(Jungle, c(1,2), mean)


# Scale data and check output ---------------------------------------------

Z <- scale(Full)
rownames(Z) <- Metabs[,1]
Trans <- t(Z)

# Color function to generate green-red heat maps --------------------------

color <- function(n = 50, low.col = 0.45, high.col=1, saturation = 1) { 
  if (n < 2) stop("n must be greater than 2")
  n1 <- n%/%2
  n2 <- n - n1
  c(hsv(low.col, saturation, seq(1,0,length=n1)), hsv(high.col, saturation, seq(0,1,length=n2))) 
}

# Clustering --------------------------------------------------------------

#Cluster rows by Pearson correlation

hmetab <- hclust(as.dist(1-cor(Trans, method="pearson")), method="complete") 

#Clusters columns by Pearson correlation

hsample <- hclust(as.dist(1-cor(Z, method="pearson")), method="complete") 

#Plot the data table as heatmap and the cluster results as dendrograms

heatmap(Z, Rowv=as.dendrogram(hmetab), Colv=as.dendrogram(hsample), col=my.colorFct(), scale="none") 
heatmap(Z, Rowv=as.dendrogram(hmetab), Colv=as.dendrogram(hsample), col=my.colorFct(), scale="none") 


# Cut tree ----------------------------------------------------------------

#Cut the tree at specific height and color the corresponding clusters in the heatmap color bar

mycm <- cutree(hmetab, h=max(hmetab$height)/1.5)
mycolhmetab <- sample(rainbow(256))
mycolhmetab <- mycolhmetab[as.vector(mycm)]
mycs <- cutree(hsample, h=max(hsample$height)/1.3)
mycolhsample <- sample(rainbow(256))
mycolhsample <- mycolhsample[as.vector(mycs)]
par(oma=c(5,2,2,5))
heatmap(Z, Rowv=as.dendrogram(hmetab), Colv=as.dendrogram(hsample), col=color(), scale="none", RowSideColors=mycolhmetab, ColSideColors=mycolhsample) 

# Include PAM results -----------------------------------------------------

#Compare PAM clustering results with hierarchical clustering by labeling it in heatmap color bar

clu19 <- pam(Trans, k=19)
clu3 <- pam(Z, k=3)

pamsample <- sample(rainbow(256)) 
pamsample <- pamsample[as.vector(clu19$clustering)]
pammetab <- sample(rainbow(256))
pammetab <- pammetab[as.vector(clu3$clustering)]
par(oma=c(5,2,2,5))
heatmap(Trans, Colv=as.dendrogram(hmetab), col=color(), scale="none", RowSideColors=pamsample, ColSideColors=pammetab) 

#Test different cutoffs

clu3s <- pam(Trans, k=3)
clu3m <- pam(Z, k=3)

pamsample <- sample(rainbow(256)) 
pamsample <- pamsample[as.vector(clu3s$clustering)]
pammetab <- sample(rainbow(256))
pammetab <- pammetab[as.vector(clu3m$clustering)]
par(oma=c(5,2,2,5))
heatmap(Trans, Colv=as.dendrogram(hmetab), col=color(), scale="none", RowSideColors=pamsample, ColSideColors=pammetab) 