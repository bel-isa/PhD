setwd("Y:/Isabel Orf/R-3.1.2")
# Load data ---------------------------------------------------------------

Metabs <- read.delim("Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/15-03-09_vergleich-arasyn.txt")

View(Metabs)


# Random forest x1 --------------------------------------------------------

library(missForest)

Wald <- missForest(Metabs, maxiter=100, ntree=100)

Full <- Wald$ximp[,-1]
Full <- data.matrix(Full)
class(Full) <- "numeric"
row.names(Full) <- Metabs[,1]

View(Full)
dim(Full)
class(Full[,1])


# Random forest x100 ------------------------------------------------------

library(abind)

Jungle <- abind(replicate(100, missForest(Metabs[,-1])$ximp, simplify=FALSE), along=3)

Full <- apply(Jungle, c(1,2), mean)


# Z transform data --------------------------------------------------------

Z <- scale(Full)

test <- apply(Full, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X))) #manually transform into a range between 0 and 1
testz <- scale(test, scale=F)

View(Z)

Mean <- colMeans(Z)

View(Mean)

boxplot(Z, las=2, xaxt="n")
abline(h = 0, col = "red")

boxplot(test, las=2, xaxt="n")
abline(h = 0, col = "red")
boxplot(testz, las=2, xaxt="n")
abline(h = 0, col = "red")


# log transform and median center data (Heike Sprengers function) ------------------------------------------------------

func_log_transform <- function(samples) {
  
  median_samples <- apply(samples, 1, median, na.rm=TRUE) # calculate median rowwise --> per analyte over all samples
  
  transformed <- data.frame(matrix(rep(NA, nrow(samples)*ncol(samples)), nrow=nrow(samples)))
  
  for (i  in 1:length(samples[1,])) {
    for (j in 1:length(samples[ ,1])) {
      transformed[j,i] <- log10(samples[j,i]/median_samples[j]) # calculate log10 of ratio of value/median
    }
  }
  colnames(transformed) <- colnames(samples)
  rownames(transformed) <- rownames(samples)
  return(transformed)
}

L <- func_log_transform(Full)

boxplot(L, las=2, xaxt="n")
abline(h = 0, col = "red")
