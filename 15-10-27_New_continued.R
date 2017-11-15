setwd("Y:/Isabel Orf/Projects/01_Ara_Syn/Metabolome/2015/2015-10-27_Results")
library(pvclust)

# Repeat initial clustering with condition averages -----------------------

data <- cbind(exp1HC, exp1LC1, exp1LC2, exp2HC, exp2LC1, exp2LC2, exp3HC, exp3LC1, exp3LC2, exp4HC, exp4LC1, exp4LC2, exp4LC3, exp5HC, exp5LC1, exp5LC2, exp5LC3, exp6NC1, exp6NC2, exp6NC3, exp6NC4, exp6NC5, exp6NC6, exp6NC7, exp6VC1, exp6VC2, exp6VC3, exp6VC4, exp6VC5, exp6VC6, exp6VC7, exp7NC1, exp7NC2, exp7NC3, exp7NC4, exp7NC5, exp7NC6, exp7NC7, exp7HC1, exp7HC2, exp7HC3, exp7HC4, exp7HC5, exp7HC6, exp7HC7)
data1 <- cbind(exp1HC, exp1LC1, exp1LC2, exp2HC, exp2LC1, exp2LC2, exp3HC, exp3LC1, exp3LC2, exp4HC, exp4LC1, exp4LC2, exp4LC3, exp6NC1, exp6NC2, exp6NC3, exp6NC4, exp6NC5, exp6VC1, exp6VC2, exp6VC3, exp6VC4, exp6VC5, exp7NC1, exp7NC2, exp7NC3, exp7NC4, exp7NC5, exp7HC1, exp7HC2, exp7HC3, exp7HC4, exp7HC5)
colnames(data) <- names
names1 <- c("Syn exp1 HC", "Syn exp1 3h LC", "Syn exp1 24h LC", "Syn exp2 HC", "Syn exp2 3h LC", "Syn exp2 24h LC", "Syn exp3 HC", "Syn exp3 3h LC", "Syn exp3 24h LC", "Ara exp4 light HC", "Ara exp4 light 1d LC", "Ara exp4 light 3d LC", "Ara exp4 light 5d LC", "Ara exp6 2h NC", "Ara exp6 4h NC", "Ara exp6 6h NC", "Ara exp6 8h NC", "Ara exp6 10h NC", "Ara exp6 2h VC", "Ara exp6 4h VC", "Ara exp6 6h VC", "Ara exp6 8h VC", "Ara exp6 10h VC", "Ara exp7 2h NC", "Ara exp7 4h NC", "Ara exp7 6h NC", "Ara exp7 8h NC", "Ara exp7 10h NC", "Ara exp7 2h HC", "Ara exp7 4h HC", "Ara exp7 6h HC", "Ara exp7 8h HC", "Ara exp7 10h HC")
colnames(data1) <- names1
data2 <- log10(data1)
data3 <- scale(data1)

par(oma = c(5,1,1,1))
samp2 <- hclust(as.dist(1-cor(data, method="pearson")), method="average")
metab2 <- hclust(as.dist(1-cor(t(data), method="pearson")), method="average")
heatmap(data, main="Condition averages", Rowv=as.dendrogram(metab2), Colv=as.dendrogram(samp2), scale="none", cexCol =0.8)
plot(samp2, main="Condition averages", hang=-1, cex=0.6)

samp3 <- hclust(as.dist(1-cor(data1, method="pearson")), method="average")
metab3 <- hclust(as.dist(1-cor(t(data1), method="pearson")), method="average")
heatmap(data1, main="Condition averages without Ara dark samples", Rowv=as.dendrogram(metab3), Colv=as.dendrogram(samp3), scale="none", cexCol =0.8)
plot(samp3, main="Condition averages without Ara dark samples", hang=-1, cex=0.8)

samp4 <- hclust(as.dist(1-cor(data2, method="pearson")), method="average")
metab4 <- hclust(as.dist(1-cor(t(data2), method="pearson")), method="average")
heatmap(data2, main="Condition averages without Ara dark samples log transformed", Rowv=as.dendrogram(metab4), Colv=as.dendrogram(samp4), scale="none", cexCol =0.8)
plot(samp4, main="Condition averages without Ara dark samples log transformed", hang=-1, cex=0.8)

samp5 <- hclust(as.dist(1-cor(data3, method="pearson")), method="average")
metab5 <- hclust(as.dist(1-cor(t(data3), method="pearson")), method="average")
heatmap(data3, main="Condition averages without Ara dark samples scaled", Rowv=as.dendrogram(metab5), Colv=as.dendrogram(samp5), col=color(), scale="none", cexCol =0.8)
plot(samp5, main="Condition averages without Ara dark samples scaled", hang=-1, cex=0.8)


# Number of clusters ------------------------------------------------------
### See also File 15-10-27_Number-of_clusters.R

num.clus(data3)
num.clus(t(data3))

data3.clu.samp <- pam(t(data3), k=5)
data3.clu.metab <- pam(data3, k=6)

pamsamp <- sample(rainbow(5)) 
pamsamp <- pamsamp[as.vector(data3.clu.samp$clustering)]
pammeta <- sample(rainbow(6))
pammeta <- pammeta[as.vector(data3.clu.metab$clustering)]
par(oma=c(7,2,2,5))
heatmap(data3, Rowv=as.dendrogram(metab5), Colv=as.dendrogram(samp5), col=color(), scale="none", RowSideColors=pammeta, ColSideColors=pamsamp, main="Biclustering (pearson, average) and PAM (samp k=7, metab k=6)")

Samp.PAM <- cbind(t(data3), data3.clu.samp$clustering)
Metab.PAM <- cbind(data3, data3.clu.metab$clustering)

Samp.clu1 <- Samp.PAM[Samp.PAM[,ncol(Samp.PAM)]==1,]
rownames(Samp.clu1)
Samp.clu2 <- Samp.PAM[Samp.PAM[,ncol(Samp.PAM)]==2,]
rownames(Samp.clu2)
Samp.clu3 <- Samp.PAM[Samp.PAM[,ncol(Samp.PAM)]==3,]
rownames(Samp.clu3)
Samp.clu4 <- Samp.PAM[Samp.PAM[,ncol(Samp.PAM)]==4,]
rownames(Samp.clu4)
Samp.clu5 <- Samp.PAM[Samp.PAM[,ncol(Samp.PAM)]==5,]
rownames(Samp.clu5)
Samp.clu6 <- Samp.PAM[Samp.PAM[,ncol(Samp.PAM)]==6,]
rownames(Samp.clu6)
Samp.clu7 <- Samp.PAM[Samp.PAM[,ncol(Samp.PAM)]==7,]
rownames(Samp.clu7)

# Analyze CV correlation, optimize graphics for publication ---------------

par(mfrow = c(1,3))
plot(CV.S,CV.A, col= c("red","black")[HCLC], pch = c(15,17) [HCLC], xlab = "Synechocystis", ylab = "Arabidopsis", main= "Coefficient of variation", ylim=c(0,1.4), xlim=c(0,1.4))
metabnames <- rownames(Full)
leg.txt <- c("LC", "HC")
legend(list(x=0,y=1.3), legend=leg.txt, col= c("red","black"), pch = c(15,17))
text(CV.S, CV.A, metabnames, cex=0.5, pos=4, col="black")

reg <- lm(CV.A ~ CV.S-1)
abline(reg)
summary(reg)
impr <- rstandard(reg)
names(impr) <- c(metabnames,metabnames)
summary(impr)

data4 <- abs(impr) < 1

plot(CV.S[data4],CV.A[data4], ylim=c(0,1.4), xlim=c(0,1.4))


reg1 <- lm(CV.A[data4] ~ CV.S[data4]-1)
abline(reg1, col=c("red"))
summary(reg1)

data5 <- abs(impr) < 0.7

plot(CV.S[data5],CV.A[data5], ylim=c(0,1.4), xlim=c(0,1.4))


reg2 <- lm(CV.A[data5] ~ CV.S[data5]-1)
abline(reg2)
summary(reg2)

par(mfrow = c(1,1))
plot(CV.S,CV.A, col= c("red","black")[HCLC], pch = c(15,17) [HCLC], xlab = "Synechocystis", ylab = "Arabidopsis", main= "Coefficient of variation", ylim=c(0,1.4), xlim=c(0,1.4))
metabnames <- rownames(Full)
leg.txt <- c("LC", "HC")
legend(list(x=0,y=1.3), legend=leg.txt, col= c("red","black"), pch = c(15,17))
text(CV.S, CV.A, metabnames, cex=0.5, pos=4, col="black")
abline(reg1)
abline(reg, col=c("red"))


# Bootstrap dendrogram for RV matrix --------------------------------------

rv1_boot <- pvclust(rv_matrix1, nboot=10000)
plot(rv1_boot)


# For final figures -------------------------------------------------------

#RV plot#

library(gplots)

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