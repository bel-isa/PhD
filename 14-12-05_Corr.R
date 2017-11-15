setwd("Y:/Isabel Orf/R-3.1.2")
# Load data ---------------------------------------------------------------

Metabs_Syn <- read.delim("Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/14-12-08_Data_Syn.txt")

View(Metabs_Syn)

# Remove additional information and store data in new variable

Metabs <- cbind(Metabs_Syn[,1], Metabs_Syn[,3:ncol(Metabs_Syn)])

View(Metabs)


# Z-transform data --------------------------------------------------------

Metabs_z <- scale(Metabs[,2:ncol(Metabs)])


# Correlation matrix ------------------------------------------------------

library(Hmisc)

Metabs_cor <- rcorr(Metabs_z, type="pearson") #can also be spearman


# To file -----------------------------------------------------------------

df.Metabs_cor.r = data.frame(Metabs_cor$r)
df.Metabs_cor.n = data.frame(Metabs_cor$n)
df.Metabs_cor.p = data.frame(Metabs_cor$p)

write.csv(df.Metabs_cor.r,'correlationmatrix.csv')
write.csv(df.Metabs_cor.p,'correlationmatrix_sign.csv')
