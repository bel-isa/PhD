
setwd("C:/Users/orf.MPIMP-GOLM/ownCloud/Project/002_Pseudomonas_Omics/Isabel/LC-MS")

library(missForest)
library(abind)
library(ggplot2)
library(mixOmics)

# Data

Metabs <- read.delim("Secondary_all1.txt", sep="\t", header=TRUE, na.strings="NA")
tmp <- Metabs[,-1]
tmp1 <- as.matrix(tmp)
row.names(tmp1) <- Metabs[,1]
Metabs <- tmp1

# Impute missing values using random forest function # Not required here, because Hezi replaced missing values by zero... to be discussed 

#jungle <- function(x){
#  tmp <- abind(replicate(100,missForest(x)$ximp,simplify=F),along=3)
#  tmp1 <- apply(tmp, c(1,2), mean)
#  return(tmp1)
#}

# The replacement is performed 100 times and the average over all replacements is used for further calculations

#Metabs <- jungle(Metabs)

Metabs_scale <- scale(t(Metabs), center=TRUE, scale=TRUE)

genotype <- rep("darkgreen", nrow(Metabs_scale))
genotype[grep("C24", fixed = T, rownames(Metabs_scale))] <- "red"

# Perform the PLS-DA
plsda <- plsda(Metabs_scale, genotype)

pdf("PLS-DA.pdf", width = 10, height = 7)
plotIndiv(plsda, col = genotype, ind.names = T, cex = 1, geom = "point", 
          X.label = "PLS-DA dim1", Y.label = "PLS-DA dim2", ellipse=T, title="", style ="graphics")
abline(h = 0, v = 0, lty = 2, col = "grey")
dev.off()
pdf("Loadings_comp1.pdf", width = 60, height = 60)
plotLoadings(plsda, comp = 1)
dev.off()
pdf("Loadings_comp2.pdf", width = 60, height = 60)
plotLoadings(plsda, comp = 2)
dev.off()

svg("PLS-DA.svg", width = 10, height = 7)
plotIndiv(plsda, col = genotype, ind.names = T, cex = 1, geom = "point", 
          X.label = "PLS-DA dim1", Y.label = "PLS-DA dim2", ellipse=T, title="", style ="graphics")
abline(h = 0, v = 0, lty = 2, col = "grey")
dev.off()
svg("Loadings_comp1.svg", width = 60, height = 60)
plotLoadings(plsda, comp = 1)
dev.off()
svg("Loadings_comp2.svg", width = 60, height = 60)
plotLoadings(plsda, comp = 2)
dev.off()

write.table(plsda$loadings$X,file= 'Loadings.txt', sep="\t")

# Group data into genotypes
Col0 <- Metabs_scale[c(grep("Col.0", rownames(Metabs_scale))),]
C24 <- Metabs_scale[c(grep("C24", rownames(Metabs_scale))),]

treatment <- rep("darkgreen", nrow(Col0))
treatment[grep("Pst", fixed = T, rownames(Col0))] <- "red"

plsda1 <- plsda(Col0, treatment)

pdf("PLS-DA_Col0.pdf", width = 10, height = 7)
plotIndiv(plsda1, col = treatment, ind.names = T, cex = 1, X.label = "PLS-DA dim1", Y.label = "PLS-DA dim2", ellipse=T, 
          title="", style ="graphics")
abline(h = 0, v = 0, lty = 2, col = "grey")
dev.off()
pdf("Loadings_comp1_Col0.pdf", width = 60, height = 60)
plotLoadings(plsda1, comp = 1)
dev.off()
pdf("Loadings_comp2_Col0.pdf", width = 60, height = 60)
plotLoadings(plsda1, comp = 2)
dev.off()

svg("PLS-DA_Col0.svg", width = 10, height = 7)
plotIndiv(plsda1, col = treatment, ind.names = T, cex = 1, X.label = "PLS-DA dim1", Y.label = "PLS-DA dim2", ellipse=T, 
          title="", style ="graphics")
abline(h = 0, v = 0, lty = 2, col = "grey")
dev.off()
svg("Loadings_comp1_Col0.svg", width = 60, height = 60)
plotLoadings(plsda1, comp = 1)
dev.off()
svg("Loadings_comp2_Col0.svg", width = 60, height = 60)
plotLoadings(plsda1, comp = 2)
dev.off()

write.table(plsda1$loadings$X,file= 'Loadings_Col0.txt', sep="\t")

treatment <- rep("darkgreen", nrow(C24))
treatment[grep("Pst", fixed = T, rownames(C24))] <- "red"

plsda2 <- plsda(C24, treatment)

pdf("PLS-DA_C24.pdf", width = 10, height = 7)
plotIndiv(plsda2, col = treatment, ind.names = T, cex = 1, X.label = "PLS-DA dim1", Y.label = "PLS-DA dim2", ellipse=T, 
          title="", style ="graphics")
abline(h = 0, v = 0, lty = 2, col = "grey")
dev.off()
pdf("Loadings_comp1_C24.pdf", width = 60, height = 60)
plotLoadings(plsda2, comp = 1)
dev.off()
pdf("Loadings_comp2_C24.pdf", width = 60, height = 60)
plotLoadings(plsda2, comp = 2)
dev.off()

svg("PLS-DA_C24.svg", width = 10, height = 7)
plotIndiv(plsda2, col = treatment, ind.names = T, cex = 1, X.label = "PLS-DA dim1", Y.label = "PLS-DA dim2", ellipse=T, 
          title="", style ="graphics")
abline(h = 0, v = 0, lty = 2, col = "grey")
dev.off()
svg("Loadings_comp1_C24.svg", width = 60, height = 60)
plotLoadings(plsda2, comp = 1)
svg("Loadings_comp2_C24.svg", width = 60, height = 60)
plotLoadings(plsda2, comp = 2)
dev.off()

write.table(plsda2$loadings$X,file= 'Loadings_C24.txt', sep="\t")

# Group data into treatments
Col0_M <- Col0[c(grep("mock", rownames(Col0))),]
Col0_P <- Col0[c(grep("Pst", rownames(Col0))),]

C24_M <- C24[c(grep("mock", rownames(C24))),]
C24_P <- C24[c(grep("Pst", rownames(C24))),]

Pst <- rbind(Col0_P, C24_P)

genotype <- rep("darkgreen", nrow(Pst))
genotype[grep("C24", fixed = T, rownames(Pst))] <- "red"

plsda3 <- plsda(Pst, genotype)

pdf("PLS-DA_Pst.pdf", width = 10, height = 7)
plotIndiv(plsda3, col = genotype, ind.names = T, cex = 1, X.label = "PLS-DA dim1", Y.label = "PLS-DA dim2", ellipse=T, 
          title="", style ="graphics")
abline(h = 0, v = 0, lty = 2, col = "grey")
dev.off()
pdf("Loadings_comp1_Pst.pdf", width = 60, height = 60)
plotLoadings(plsda3, comp = 1)
dev.off()
pdf("Loadings_comp2_Pst.pdf", width = 60, height = 60)
plotLoadings(plsda3, comp = 2)
dev.off()

svg("PLS-DA_Pst.svg", width = 10, height = 7)
plotIndiv(plsda3, col = genotype, ind.names = T, cex = 1, X.label = "PLS-DA dim1", Y.label = "PLS-DA dim2", ellipse=T, 
          title="", style ="graphics")
abline(h = 0, v = 0, lty = 2, col = "grey")
dev.off()
svg("Loadings_comp1_Pst.svg", width = 60, height = 60)
plotLoadings(plsda3, comp = 1)
dev.off()
svg("Loadings_comp2_Pst.svg", width = 60, height = 60)
plotLoadings(plsda3, comp = 2)
dev.off()

write.table(plsda3$loadings$X,file= 'Loadings_Pst.txt', sep="\t")

Mock <- rbind(Col0_M, C24_M)

genotype <- rep("darkgreen", nrow(Mock))
genotype[grep("C24", fixed = T, rownames(Mock))] <- "red"

plsda4 <- plsda(Mock, genotype)

pdf("PLS-DA_Mock.pdf", width = 10, height = 7)
plotIndiv(plsda4, col = genotype, ind.names = T, cex = 1, X.label = "PLS-DA dim1", Y.label = "PLS-DA dim2", ellipse=T, 
          title="", style ="graphics")
abline(h = 0, v = 0, lty = 2, col = "grey")
dev.off()
pdf("Loadings_comp1_Mock.pdf", width = 60, height = 60)
plotLoadings(plsda4, comp = 1)
dev.off()
pdf("Loadings_comp2_Mock.pdf", width = 60, height = 60)
plotLoadings(plsda4, comp = 2)
dev.off()

svg("PLS-DA_Mock.svg", width = 10, height = 7)
plotIndiv(plsda4, col = genotype, ind.names = T, cex = 1, X.label = "PLS-DA dim1", Y.label = "PLS-DA dim2", ellipse=T, 
          title="", style ="graphics")
abline(h = 0, v = 0, lty = 2, col = "grey")
dev.off()
svg("Loadings_comp1_Mock.svg", width = 60, height = 60)
plotLoadings(plsda4, comp = 1)
dev.off()
svg("Loadings_comp2_Mock.svg", width = 60, height = 60)
plotLoadings(plsda4, comp = 2)
dev.off()

write.table(plsda4$loadings$X,file= 'Loadings_Mock.txt', sep="\t")
