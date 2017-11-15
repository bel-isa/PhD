setwd("Y:/Isabel Orf/R-3.1.2")
# Load data ---------------------------------------------------------------

Metabs <- read.delim("Y:/Isabel Orf/Projects/Ara_Syn/Opinion/15-02-20_vergleich-arasyn.txt")

View(Metabs)


# Z-transform data --------------------------------------------------------

Metabs_z <- scale(Metabs[,2:ncol(Metabs)])

View(Metabs_z)

# Correlation matrix ------------------------------------------------------

library(Hmisc)

Metabs_p <- rcorr(Metabs_z, type="pearson") #can also be spearman
Metabs_s <- rcorr(Metabs_z, type="spearman")

# To file -----------------------------------------------------------------

df.Metabs_p.r = data.frame(Metabs_p$r)
df.Metabs_p.n = data.frame(Metabs_p$n)
df.Metabs_p.p = data.frame(Metabs_p$p)

write.csv(df.Metabs_p.r,'correlationmatrix.csv')
write.csv(df.Metabs_p.p,'correlationmatrix_sign.csv')
