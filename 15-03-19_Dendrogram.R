#######################################
####  Analysis of reduced dataset  ####
#######################################

# Set directory and load packages -----------------------------------------

setwd("Y:/Isabel Orf/R-3.1.2")
library(cluster)
require(graphics)
require(utils)

# Load data ---------------------------------------------------------------

Full <- read.delim("Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-10_Results/2015-03-10_Full.txt")

# Scale data ---------------------------------------------

Z <- scale(Full)
rownames(Z) <- Metabs[,1]
Trans <- t(Z)

# Create reduced dataset --------------------------------------------------

Complete <- cbind(Trans, numm, clu4$clustering)
x <- 2
New <- Complete[Complete[,ncol(Complete)]!=x,]

New1 <- New[,-(31:33)]

exp <- New[,31]
cond <- New[,32]
numm <- cbind(exp, cond)

# Perform pearson clustering in PAM clusters ------------------------------

Reduced.clu <- cbind(New1, Red.clu3$clustering)

y <- 3
z <- 1

cluster1 <- Reduced.clu[Reduced.clu[,ncol(Reduced.clu)]==z,]
cluster2 <- Reduced.clu[Reduced.clu[,ncol(Reduced.clu)]==x,]
cluster3 <- Reduced.clu[Reduced.clu[,ncol(Reduced.clu)]==y,]

cluster1a <- t(cluster1[,-31])
cluster2a <- t(cluster2[,-31])
cluster3a <- t(cluster3[,-31])

baum1 <- hclust(as.dist(1-cor(cluster1a, method="pearson")), method="complete")
plot(baum1, main="Pearson complete linkage, PAM cluster 1", hang=-1, cex=0.6)

baum2 <- hclust(as.dist(1-cor(cluster2a, method="pearson")), method="complete")
plot(baum2, main="Pearson complete linkage, PAM cluster 2", hang=-1, cex=0.6)

baum3 <- hclust(as.dist(1-cor(cluster3a, method="pearson")), method="complete")
plot(baum3, main="Pearson complete linkage, PAM cluster 3", hang=-1, cex=0.6)

baum1x <- hclust(as.dist(1-cor(cluster1a, method="spearman")), method="complete")
plot(baum1x, main="Spearman complete linkage, PAM cluster 1", hang=-1, cex=0.6)

baum2x <- hclust(as.dist(1-cor(cluster2a, method="spearman")), method="complete")
plot(baum2x, main="Spearman complete linkage, PAM cluster 2", hang=-1, cex=0.6)

baum3x <- hclust(as.dist(1-cor(cluster3a, method="spearman")), method="complete")
plot(baum3x, main="Spearman complete linkage, PAM cluster 3", hang=-1, cex=0.6)

# Count exp or cond Ara/Syn per cluster -----------------------------------

c1 <- data.frame(pam1 = pam(New1, k=1)$clustering)
c1a <- cbind(c1, numm)
write.table(c1a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_C1.txt", sep="\t")

c2 <- data.frame(pam2 = pam(New1, k=2)$clustering)
c2a <- cbind(c2, numm)
write.table(c2a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_C2.txt", sep="\t")

c3 <- data.frame(pam3 = pam(New1, k=3)$clustering)
c3a <- cbind(c3, numm)
write.table(c3a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_C3.txt", sep="\t")

c4 <- data.frame(pam4 = pam(New1, k=4)$clustering)
c4a <- cbind(c4, numm)
write.table(c4a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_C4.txt", sep="\t")

c5 <- data.frame(pam5 = pam(New1, k=5)$clustering)
c5a <- cbind(c5, numm)
write.table(c5a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_C5.txt", sep="\t")

c6 <- data.frame(pam6 = pam(New1, k=6)$clustering)
c6a <- cbind(c6, numm)
write.table(c6a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_C6.txt", sep="\t")

c7 <- data.frame(pam7 = pam(New1, k=7)$clustering)
c7a <- cbind(c7, numm)
write.table(c7a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_C7.txt", sep="\t")

c8 <- data.frame(pam8 = pam(New1, k=8)$clustering)
c8a <- cbind(c8, numm)
write.table(c8a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_C8.txt", sep="\t")

c9 <- data.frame(pam9 = pam(New1, k=9)$clustering)
c9a <- cbind(c9, numm)
write.table(c9a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_C9.txt", sep="\t")

c10 <- data.frame(pam10 = pam(New1, k=10)$clustering)
c10a <- cbind(c10, numm)
write.table(c10a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_C10.txt", sep="\t")

c11 <- data.frame(pam11 = pam(New1, k=11)$clustering)
c11a <- cbind(c11, numm)
write.table(c11a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_C11.txt", sep="\t")

c12 <- data.frame(pam12 = pam(New1, k=12)$clustering)
c12a <- cbind(c12, numm)
write.table(c12a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_C12.txt", sep="\t")

c13 <- data.frame(pam13 = pam(New1, k=13)$clustering)
c13a <- cbind(c13, numm)
write.table(c13a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_C13.txt", sep="\t")

c14 <- data.frame(pam14 = pam(New1, k=14)$clustering)
c14a <- cbind(c14, numm)
write.table(c14a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_C14.txt", sep="\t")

c15 <- data.frame(pam15 = pam(New1, k=15)$clustering)
c15a <- cbind(c15, numm)
write.table(c15a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_C15.txt", sep="\t")

c16 <- data.frame(pam16 = pam(New1, k=16)$clustering)
c16a <- cbind(c16, numm)
write.table(c16a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_C16.txt", sep="\t")

c17 <- data.frame(pam17 = pam(New1, k=17)$clustering)
c17a <- cbind(c17, numm)
write.table(c17a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_C17.txt", sep="\t")

c18 <- data.frame(pam18 = pam(New1, k=18)$clustering)
c18a <- cbind(c18, numm)
write.table(c18a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_C18.txt", sep="\t")

c19 <- data.frame(pam19 = pam(New1, k=19)$clustering)
c19a <- cbind(c19, numm)
write.table(c19a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_C19.txt", sep="\t")

c20 <- data.frame(pam20 = pam(New1, k=20)$clustering)
c20a <- cbind(c20, numm)
write.table(c20a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_C20.txt", sep="\t")

c21 <- data.frame(pam21 = pam(New1, k=21)$clustering)
c21a <- cbind(c21, numm)
write.table(c21a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_C21.txt", sep="\t")

c22 <- data.frame(pam22 = pam(New1, k=22)$clustering)
c22a <- cbind(c22, numm)
write.table(c22a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_C22.txt", sep="\t")

c23 <- data.frame(pam23 = pam(New1, k=23)$clustering)
c23a <- cbind(c23, numm)
write.table(c23a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_C23.txt", sep="\t")

c24 <- data.frame(pam24 = pam(New1, k=24)$clustering)
c24a <- cbind(c24, numm)
write.table(c24a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_C24.txt", sep="\t")

c25 <- data.frame(pam25 = pam(New1, k=25)$clustering)
c25a <- cbind(c25, numm)
write.table(c25a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_C25.txt", sep="\t")

c26 <- data.frame(pam26 = pam(New1, k=26)$clustering)
c26a <- cbind(c26, numm)
write.table(c26a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_C26.txt", sep="\t")

c27 <- data.frame(pam27 = pam(New1, k=27)$clustering)
c27a <- cbind(c27, numm)
write.table(c27a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_C27.txt", sep="\t")

c28 <- data.frame(pam28 = pam(New1, k=28)$clustering)
c28a <- cbind(c28, numm)
write.table(c28a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_C28.txt", sep="\t")

c29 <- data.frame(pam29 = pam(New1, k=29)$clustering)
c29a <- cbind(c29, numm)
write.table(c29a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_C29.txt", sep="\t")

c30 <- data.frame(pam30 = pam(New1, k=30)$clustering)
c30a <- cbind(c30, numm)
write.table(c30a, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_C30.txt", sep="\t")

write.table(table(c1a$pam, c1a$exp), "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_exp_groups.txt", sep="\t")
write.table(table(c2a$pam, c2a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c3a$pam, c3a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c4a$pam, c4a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c5a$pam, c5a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c6a$pam, c6a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c7a$pam, c7a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c8a$pam, c8a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c9a$pam, c9a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c10a$pam, c10a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c11a$pam, c11a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c12a$pam, c12a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c13a$pam, c13a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c14a$pam, c14a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c15a$pam, c15a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c16a$pam, c16a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c17a$pam, c17a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c18a$pam, c18a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c19a$pam, c19a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c20a$pam, c20a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c21a$pam, c21a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c22a$pam, c22a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c23a$pam, c23a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c24a$pam, c24a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c25a$pam, c25a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c26a$pam, c26a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c27a$pam, c27a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c28a$pam, c28a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c29a$pam, c29a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_exp_groups.txt", sep="\t", append=TRUE)
write.table(table(c30a$pam, c30a$exp), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_exp_groups.txt", sep="\t", append=TRUE)

write.table(table(c1a$pam, c1a$cond), "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_cond_groups.txt", sep="\t")
write.table(table(c2a$pam, c2a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c3a$pam, c3a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c4a$pam, c4a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c5a$pam, c5a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c6a$pam, c6a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c7a$pam, c7a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c8a$pam, c8a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c9a$pam, c9a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c10a$pam, c10a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c11a$pam, c11a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c12a$pam, c12a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c13a$pam, c13a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c14a$pam, c14a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c15a$pam, c15a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c16a$pam, c16a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c17a$pam, c17a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c18a$pam, c18a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c19a$pam, c19a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c20a$pam, c20a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c21a$pam, c21a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c22a$pam, c22a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c23a$pam, c23a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c24a$pam, c24a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c25a$pam, c25a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c26a$pam, c26a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c27a$pam, c27a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c28a$pam, c28a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c29a$pam, c29a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_cond_groups.txt", sep="\t", append=TRUE)
write.table(table(c30a$pam, c30a$cond), col.names=FALSE, "Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-03-19_Results/2015-03-19_Count_cond_groups.txt", sep="\t", append=TRUE)

# Plot graph

nclus <- c(1:10)

cond.ara.syn.nclus <- c(2.666666667,  2.875,	3,	2,	0.2,	0.333333333,	0,	0,	0,	0)

cond.syn.ara.nclus <- c(0.375,  0.365079365,	0.037037037,	0.03125,	0.2,	0.083333333,	0,	0,	0,	0)

plot(nclus,cond.ara.syn.nclus, type="o", main="Number of experimental conditions Ara/Syn per cluster", xlab="Number of clusters", ylab="Experimental conditions Ara/Syn per cluster")
plot(nclus,cond.syn.ara.nclus, type="o", main="Number of experimental conditions Syn/Ara per cluster", xlab="Number of clusters", ylab="Experimental conditions Syn/Ara per cluster")

exp.ara.syn.nclus <- c(1.25,  0.75,	0,	0,	0,	0,	0,	0,	0,	0)

exp.syn.ara.nclus <- c(0.8,  1.333333333,	0,	0,	0,	0,	0,	0,	0,	0)

plot(nclus,exp.ara.syn.nclus, type="o", main="Number of experiments Ara/Syn per cluster", xlab="Number of clusters", ylab="Experiments Ara/Syn per cluster")
plot(nclus,exp.syn.ara.nclus, type="o", main="Number of experiments Syn/Ara per cluster", xlab="Number of clusters", ylab="Experiments Syn/Ara per cluster")
