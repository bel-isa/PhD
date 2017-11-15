(setwd("Y:/Isabel Orf/R-3.1.2")

library(FactoMineR)

# Load data ---------------------------------------------------------------

Ara <- read.delim("Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-09-01_Results/Ara.txt", sep="\t")
Syn <- read.delim("Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-09-01_Results/Syn.txt", sep="\t")

Ara <- Ara[,-1]
Syn <- Syn[,-1]

Metabs <- read.delim("Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/15-03-09_vergleich-arasyn.txt")
rownames(Ara) <- Metabs[,1]
rownames(Syn) <- Metabs[,1]

# CV function -----------------------------------------------------

calc.CV <- function(x) {
  tmp <- as.numeric(apply(x, 1, function(z) sd(z)/mean(z)))
  return(tmp)
}

# Calculate CV matrix -----------------------------------------------------

Syn_HC <- as.matrix(Syn[,1:22])
Syn_LC <- as.matrix(Syn[,23:ncol(Syn)])

CV.Syn_HC <- calc.CV(Syn_HC)
CV.Syn_LC <- calc.CV(Syn_LC)

Ara_HC <- as.matrix(cbind(Ara[,1:6], Ara[,109:ncol(Ara)]))
Ara_LC <- as.matrix(Ara[,7:108])

CV.Ara_HC <- calc.CV(Ara_HC)
CV.Ara_LC <- calc.CV(Ara_LC)

# Plot --------------------------------------------------------------------

CV.Syn <- c(CV.Syn_LC,CV.Syn_HC)
CV.Ara <- c(CV.Ara_LC,CV.Ara_HC)

a <- rep(1,30)
b <- rep(2,30)
x <- c(a,b)

CV.Syn <- cbind(CV.Syn, x)
CV.Ara <- cbind(CV.Ara, x)

plot(CV.Syn,CV.Ara, col= c("red","black")[x])

#nice plot
plot(CV.Syn,CV.Ara, col= c("red","black")[x], pch = c(15,17) [x], xlab = "Synechocystis", ylab = "Arabidopsis", main= "Coefficient of variation")

leg.txt <- c("LC", "HC")
legend(list(x=1,y=3.7), legend=leg.txt, col= c("red","black"), pch = c(15,17))
text(CV.Syn, CV.Ara, row.names(CV.Syn), cex=0.5, pos=4, col="black") 


# Repeat with sample means ------------------------------------------------

a <- apply(Syn_HC[,1:4],1,function(x) mean(x))
b <- apply(Syn_HC[,5:7],1,function(x) mean(x))
c <- apply(Syn_HC[,8:10],1,function(x) mean(x))
d <- apply(Syn_HC[,11:13],1,function(x) mean(x))
e <- apply(Syn_HC[,14:16],1,function(x) mean(x))
f <- apply(Syn_HC[,17:19],1,function(x) mean(x))
g <- apply(Syn_HC[,20:22],1,function(x) mean(x))

Av_Syn_HC <- cbind(a,b,c,d,e,f,g)

a <- apply(Syn_LC[,1:4],1,function(x) mean(x))
b <- apply(Syn_LC[,5:7],1,function(x) mean(x))
c <- apply(Syn_LC[,8:10],1,function(x) mean(x))
d <- apply(Syn_LC[,11:13],1,function(x) mean(x))
e <- apply(Syn_LC[,14:15],1,function(x) mean(x))
f <- apply(Syn_LC[,16:17],1,function(x) mean(x))
g <- apply(Syn_LC[,18:20],1,function(x) mean(x))
h <- apply(Syn_LC[,21:24],1,function(x) mean(x))
i <- apply(Syn_LC[,25:27],1,function(x) mean(x))
j <- apply(Syn_LC[,28:30],1,function(x) mean(x))
k <- apply(Syn_LC[,31:33],1,function(x) mean(x))
l <- apply(Syn_LC[,34:36],1,function(x) mean(x))
m <- apply(Syn_LC[,37:39],1,function(x) mean(x))
n <- apply(Syn_LC[,40:42],1,function(x) mean(x))

Av_Syn_LC <- cbind(a,b,c,d,e,f,g,h,i,j,k,l,m,n)

a <- apply(Ara_HC[,1:3],1,function(x) mean(x))
b <- apply(Ara_HC[,4:6],1,function(x) mean(x))
c <- apply(Ara_HC[,7:12],1,function(x) mean(x))
d <- apply(Ara_HC[,13:18],1,function(x) mean(x))
e <- apply(Ara_HC[,19:24],1,function(x) mean(x))
f <- apply(Ara_HC[,25:30],1,function(x) mean(x))
g <- apply(Ara_HC[,31:36],1,function(x) mean(x))
h <- apply(Ara_HC[,37:42],1,function(x) mean(x))
i <- apply(Ara_HC[,43:48],1,function(x) mean(x))
j <- apply(Ara_HC[,49:54],1,function(x) mean(x))
k <- apply(Ara_HC[,55:60],1,function(x) mean(x))
l <- apply(Ara_HC[,61:68],1,function(x) mean(x))
m <- apply(Ara_HC[,69:74],1,function(x) mean(x))
n <- apply(Ara_HC[,75:80],1,function(x) mean(x))
o <- apply(Ara_HC[,81:86],1,function(x) mean(x))
p <- apply(Ara_HC[,87:92],1,function(x) mean(x))

Av_Ara_HC <- cbind(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p)

a <- apply(Ara_LC[,1:3],1,function(x) mean(x))
b <- apply(Ara_LC[,4:6],1,function(x) mean(x))
c <- apply(Ara_LC[,7:12],1,function(x) mean(x))
d <- apply(Ara_LC[,13:18],1,function(x) mean(x))
e <- apply(Ara_LC[,19:30],1,function(x) mean(x))
f <- apply(Ara_LC[,31:42],1,function(x) mean(x))
g <- apply(Ara_LC[,43:54],1,function(x) mean(x))
h <- apply(Ara_LC[,55:66],1,function(x) mean(x))
i <- apply(Ara_LC[,67:78],1,function(x) mean(x))
j <- apply(Ara_LC[,79:90],1,function(x) mean(x))
k <- apply(Ara_LC[,91:102],1,function(x) mean(x))

Av_Ara_LC <- cbind(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p)

CV.Av_Syn_HC <- apply(Av_Syn_HC,1,function(x) sd(x)/mean(x))
CV.Av_Syn_LC <- apply(Av_Syn_LC,1,function(x) sd(x)/mean(x))

CV.Av_Syn <- c(CV.Av_Syn_LC,CV.Av_Syn_HC)
CV.Av_Syn <- cbind(CV.Av_Syn, x)

CV.Av_Ara_HC <- apply(Av_Ara_HC,1,function(x) sd(x)/mean(x))
CV.Av_Ara_LC <- apply(Av_Ara_LC,1,function(x) sd(x)/mean(x))

CV.Av_Ara <- c(CV.Av_Ara_LC,CV.Av_Ara_HC)
CV.Av_Ara <- cbind(CV.Av_Ara, x)

# Plot --------------------------------------------------------------------

CV.Av_Syn <- c(CV.Av_Syn_LC,CV.Av_Syn_HC)
CV.Av_Ara <- c(CV.Av_Ara_LC,CV.Av_Ara_HC)

a <- rep(1,30)
b <- rep(2,30)
x <- c(a,b)

CV.Av_Syn <- cbind(CV.Av_Syn, x)
CV.Av_Ara <- cbind(CV.Av_Ara, x)

plot(CV.Av_Syn,CV.Av_Ara, col= c("red","black")[x])

#nice plot
plot(CV.Av_Syn,CV.Av_Ara, col= c("red","black")[x], pch = c(15,17) [x], xlab = "Synechocystis", ylab = "Arabidopsis", main= "Coefficient of variation")

leg.txt <- c("LC", "HC")
legend(list(x=1,y=2.7), legend=leg.txt, col= c("red","black"), pch = c(15,17))
text(CV.Av_Syn, CV.Av_Ara, row.names(CV.Av_Syn), cex=0.5, pos=4, col="black") 