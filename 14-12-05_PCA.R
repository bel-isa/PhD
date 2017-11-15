setwd("Y:/Isabel Orf/R-3.1.2")


# Load data matrix --------------------------------------------------------

data <- read.delim("Y:/Isabel Orf/Projects/Delta4/Manuscript/PCA_d4_ndhR/PCA_ndhR_d4_ann.txt", na.strings= "NA")


# Remove missing values/optimize ------------------------------------------

data <- na.omit(data)

rownames(data) <- data[,1]

colnames(data) <- sub("^X.", "", colnames(data))


# Calculate PCs -----------------------------------------------------------

library(pcaMethods)

pca <- pca((data[,2:10]),
           method="ppca", center=TRUE, scale='none', nPcs=9, seed=3456)



# Show results ------------------------------------------------------------

print(pca)

summary(pca)


# Create output files -----------------------------------------------------
#Can also be used to plot loadings c=(T,F)
cls <- c(rep("red"), rep("blue"), rep("orange"))

pdf('PCA_ndhR_d4_ann.pdf',7,7)
slplot(pca, pcs=c(1,2), scoresLoadings = c(F, T), scol=cls, lcex=0.8) 
slplot(pca, pcs=c(1,3), scoresLoadings = c(F, T), scol=cls, lcex=0.8)
slplot(pca, pcs=c(2,3), scoresLoadings = c(F, T), scol=cls, lcex=0.8)
slplot(pca, pcs=c(1,4), scoresLoadings = c(F, T), scol=cls, lcex=0.8)
slplot(pca, pcs=c(2,4), scoresLoadings = c(F, T), scol=cls, lcex=0.8)
slplot(pca, pcs=c(3,4), scoresLoadings = c(F, T), scol=cls, lcex=0.8)
slplot(pca, pcs=c(1,5), scoresLoadings = c(F, T), scol=cls, lcex=0.8)
slplot(pca, pcs=c(2,5), scoresLoadings = c(F, T), scol=cls, lcex=0.8)
slplot(pca, pcs=c(3,5), scoresLoadings = c(F, T), scol=cls, lcex=0.8)
slplot(pca, pcs=c(4,5), scoresLoadings = c(F, T), scol=cls, lcex=0.8)
slplot(pca, pcs=c(1,6), scoresLoadings = c(F, T), scol=cls, lcex=0.8)
slplot(pca, pcs=c(2,6), scoresLoadings = c(F, T), scol=cls, lcex=0.8)
slplot(pca, pcs=c(3,6), scoresLoadings = c(F, T), scol=cls, lcex=0.8)
slplot(pca, pcs=c(4,6), scoresLoadings = c(F, T), scol=cls, lcex=0.8)
slplot(pca, pcs=c(5,6), scoresLoadings = c(F, T), scol=cls, lcex=0.8)
slplot(pca, pcs=c(1,7), scoresLoadings = c(F, T), scol=cls, lcex=0.8)
slplot(pca, pcs=c(2,7), scoresLoadings = c(F, T), scol=cls, lcex=0.8)
slplot(pca, pcs=c(3,7), scoresLoadings = c(F, T), scol=cls, lcex=0.8)
slplot(pca, pcs=c(4,7), scoresLoadings = c(F, T), scol=cls, lcex=0.8)
slplot(pca, pcs=c(5,7), scoresLoadings = c(F, T), scol=cls, lcex=0.8)
slplot(pca, pcs=c(6,7), scoresLoadings = c(F, T), scol=cls, lcex=0.8)
slplot(pca, pcs=c(1,8), scoresLoadings = c(F, T), scol=cls, lcex=0.8)
slplot(pca, pcs=c(2,8), scoresLoadings = c(F, T), scol=cls, lcex=0.8)
slplot(pca, pcs=c(3,8), scoresLoadings = c(F, T), scol=cls, lcex=0.8)
slplot(pca, pcs=c(4,8), scoresLoadings = c(F, T), scol=cls, lcex=0.8)
slplot(pca, pcs=c(5,8), scoresLoadings = c(F, T), scol=cls, lcex=0.8)
slplot(pca, pcs=c(6,8), scoresLoadings = c(F, T), scol=cls, lcex=0.8)
slplot(pca, pcs=c(7,8), scoresLoadings = c(F, T), scol=cls, lcex=0.8)
slplot(pca, pcs=c(1,9), scoresLoadings = c(F, T), scol=cls, lcex=0.8)
slplot(pca, pcs=c(2,9), scoresLoadings = c(F, T), scol=cls, lcex=0.8)
slplot(pca, pcs=c(3,9), scoresLoadings = c(F, T), scol=cls, lcex=0.8)
slplot(pca, pcs=c(4,9), scoresLoadings = c(F, T), scol=cls, lcex=0.8)
slplot(pca, pcs=c(5,9), scoresLoadings = c(F, T), scol=cls, lcex=0.8)
slplot(pca, pcs=c(6,9), scoresLoadings = c(F, T), scol=cls, lcex=0.8)
slplot(pca, pcs=c(7,9), scoresLoadings = c(F, T), scol=cls, lcex=0.8)
slplot(pca, pcs=c(8,9), scoresLoadings = c(F, T), scol=cls, lcex=0.8)
dev.off()

# Additional plots --------------------------------------------------------

library(reshape2)
library(devtools)
install_github("ggbiplot", "vqv")

library(ggbiplot)

heatmap <- qplot(x=Var1, y=Var2, data=melt(cor(data[,2:10])), geom="tile", fill=value)

pdf('Heatmap_corr.pdf',7,5)
print(heatmap)
dev.off()

library(reshape2)
library(ggplot2)

#prcomp cannot be used if one column is constant or zero

pca.1 <- prcomp(data[,2:10], scale=T)
melted <- melt(pca.1$rotation[,1:9])

barplot <- ggplot(data=melted) +
  geom_bar(aes(x=Var1, y=value, fill=Var1), stat="identity") +
  facet_wrap(~Var2)

pdf('Barplot_PCs_ann.pdf',7,5)
print(barplot)
dev.off()

# Colorful loadings plots -------------------------------------------------

loadings <- data.frame(colnames(data[,2:10]), pca.1$r[,1:9])

pc1.2 <- qplot(x=PC1, y=PC2, data=loadings, colour=factor(colnames(data[,2:10])))
pc1.3 <- qplot(x=PC1, y=PC3, data=loadings, colour=factor(colnames(data[,2:10])))
pc2.3 <- qplot(x=PC2, y=PC3, data=loadings, colour=factor(colnames(data[,2:10])))
pc1.4 <- qplot(x=PC1, y=PC4, data=loadings, colour=factor(colnames(data[,2:10])))
pc2.4 <- qplot(x=PC2, y=PC4, data=loadings, colour=factor(colnames(data[,2:10])))
pc3.4 <- qplot(x=PC3, y=PC4, data=loadings, colour=factor(colnames(data[,2:10])))
pc1.5 <- qplot(x=PC1, y=PC5, data=loadings, colour=factor(colnames(data[,2:10])))
pc2.5 <- qplot(x=PC2, y=PC5, data=loadings, colour=factor(colnames(data[,2:10])))
pc3.5 <- qplot(x=PC3, y=PC5, data=loadings, colour=factor(colnames(data[,2:10])))
pc4.5 <- qplot(x=PC4, y=PC5, data=loadings, colour=factor(colnames(data[,2:10])))
pc1.6 <- qplot(x=PC1, y=PC6, data=loadings, colour=factor(colnames(data[,2:10])))
pc2.6 <- qplot(x=PC2, y=PC6, data=loadings, colour=factor(colnames(data[,2:10])))
pc3.6 <- qplot(x=PC3, y=PC6, data=loadings, colour=factor(colnames(data[,2:10])))
pc4.6 <- qplot(x=PC4, y=PC6, data=loadings, colour=factor(colnames(data[,2:10])))
pc5.6 <- qplot(x=PC5, y=PC6, data=loadings, colour=factor(colnames(data[,2:10])))
pc1.7 <- qplot(x=PC1, y=PC7, data=loadings, colour=factor(colnames(data[,2:10])))
pc2.7 <- qplot(x=PC2, y=PC7, data=loadings, colour=factor(colnames(data[,2:10])))
pc3.7 <- qplot(x=PC3, y=PC7, data=loadings, colour=factor(colnames(data[,2:10])))
pc4.7 <- qplot(x=PC4, y=PC7, data=loadings, colour=factor(colnames(data[,2:10])))
pc5.7 <- qplot(x=PC5, y=PC7, data=loadings, colour=factor(colnames(data[,2:10])))
pc6.7 <- qplot(x=PC6, y=PC7, data=loadings, colour=factor(colnames(data[,2:10])))
pc1.8 <- qplot(x=PC1, y=PC8, data=loadings, colour=factor(colnames(data[,2:10])))
pc2.8 <- qplot(x=PC2, y=PC8, data=loadings, colour=factor(colnames(data[,2:10])))
pc3.8 <- qplot(x=PC3, y=PC8, data=loadings, colour=factor(colnames(data[,2:10])))
pc4.8 <- qplot(x=PC4, y=PC8, data=loadings, colour=factor(colnames(data[,2:10])))
pc5.8 <- qplot(x=PC5, y=PC8, data=loadings, colour=factor(colnames(data[,2:10])))
pc6.8 <- qplot(x=PC6, y=PC8, data=loadings, colour=factor(colnames(data[,2:10])))
pc7.8 <- qplot(x=PC7, y=PC8, data=loadings, colour=factor(colnames(data[,2:10])))
pc1.9 <- qplot(x=PC1, y=PC9, data=loadings, colour=factor(colnames(data[,2:10])))
pc2.9 <- qplot(x=PC2, y=PC9, data=loadings, colour=factor(colnames(data[,2:10])))
pc3.9 <- qplot(x=PC3, y=PC9, data=loadings, colour=factor(colnames(data[,2:10])))
pc4.9 <- qplot(x=PC4, y=PC9, data=loadings, colour=factor(colnames(data[,2:10])))
pc5.9 <- qplot(x=PC5, y=PC9, data=loadings, colour=factor(colnames(data[,2:10])))
pc6.9 <- qplot(x=PC6, y=PC9, data=loadings, colour=factor(colnames(data[,2:10])))
pc7.9 <- qplot(x=PC7, y=PC9, data=loadings, colour=factor(colnames(data[,2:10])))
pc8.9 <- qplot(x=PC8, y=PC9, data=loadings, colour=factor(colnames(data[,2:10])))

pdf('PCA_ndhR_d4_ann_color.pdf',7,5)
pc1.2
pc1.3
pc2.3
pc1.4
pc2.4
pc3.4
pc1.5
pc2.5
pc3.5
pc4.5
pc1.6
pc2.6
pc3.6
pc4.6
pc5.6
pc1.7
pc2.7
pc3.7
pc4.7
pc5.7
pc6.7
pc1.8
pc2.8
pc3.8
pc4.8
pc5.8
pc6.8
pc7.8
pc1.9
pc2.9
pc3.9
pc4.9
pc5.9
pc6.9
pc7.9
pc8.9
dev.off()