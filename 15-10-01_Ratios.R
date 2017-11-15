setwd("Y:/Isabel Orf/R-3.1.2")
 
library(FactoMineR)

 
# Load data ---------------------------------------------------------------
 
Ara <- read.delim("Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-09-01_Results/Ara.txt", sep="\t")
Syn <- read.delim("Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-09-01_Results/Syn.txt", sep="\t")
 
Ara <- Ara[,-1]
Syn <- Syn[,-1]
 
Metabs <- read.delim("Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/15-03-09_vergleich-arasyn.txt")
 
rownames(Ara) <- Metabs[,1]
rownames(Syn) <- Metabs[,1]

# Split Syn data into smaller sets --------------------------------------------

Syn_HC <- as.matrix(Syn[,1:22])
Syn_LC <- as.matrix(Syn[,23:ncol(Syn)])

SHC1 <- apply(Syn_HC[,1:4],1, function(x) mean(x)) #10076if
SHC2 <- apply(Syn_HC[,5:10],1, function(x) mean(x)) #12246if
SHC3 <- apply(Syn_HC[,11:19],1, function(x) mean(x)) #13122mo
SHC4 <- apply(Syn_HC[,20:22],1, function(x) mean(x)) #14217if
S3C1 <- apply(Syn_LC[,1:4],1, function(x) mean(x)) #10076if
S3C2 <- apply(Syn_LC[,5:10],1, function(x) mean(x)) #12246if
S3C3 <- apply(Syn_LC[,11:17],1, function(x) mean(x)) #13122mo
S3C4 <- apply(Syn_LC[,18:20],1, function(x) mean(x)) #14217if
SLC1 <- apply(Syn_LC[,21:24],1, function(x) mean(x)) #10076if
SLC2 <- apply(Syn_LC[,25:30],1, function(x) mean(x)) #12246if
SLC3 <- apply(Syn_LC[,31:39],1, function(x) mean(x)) #13122mo
SLC4 <- apply(Syn_LC[,40:42],1, function(x) mean(x)) #14217if

# Split Ara data into smaller sets --------------------------------------------

Ara_HC <- as.matrix(cbind(Ara[,1:6], Ara[,109:ncol(Ara)]))
Ara_LC <- as.matrix(Ara[,7:108])

AHC1 <- apply(Ara_HC[,1:3],1, mean) #Exp 1
AHC2 <- apply(Ara_HC[,7:48],1, mean #Exp 2
AVC <- apply(Ara_HC[,49:ncol(Ara_HC)],1, mean #Exp 2
ALC1 <- apply(Ara_LC[,1:3],1, mean #Exp 1
ALC3 <- apply(Ara_LC[,7:9],1, mean #Exp 1
ALC5 <- apply(Ara_LC[,13:15],1, mean) #Exp 1
ANC <- apply(Ara_LC[,19:ncol(Ara_LC)],1, mean) #Exp 2

# Calculate ratio matrices --------------------------------------------

ratioMat <- function(x) {
  tmp <- matrix(nrow = length(x), ncol = length(x))
  for (i in 1:length(x)){
    for (j in 1:length(x)){
      tmp[i,j] <- (x[i]/x[j])
    }
  }
  return(tmp)
}

# Calculate ratio matrices Syn --------------------------------------------

SHC1Mat <- ratioMat(SHC1) #10076if
SHC2Mat <- ratioMat(SHC2) #12246if
SHC3Mat <- ratioMat(SHC3) #13122mo
SHC4Mat <- ratioMat(SHC4) #14217if
S3C1Mat <- ratioMat(S3C1) #10076if
S3C2Mat <- ratioMat(S3C2) #12246if
S3C3Mat <- ratioMat(S3C3) #13122mo
S3C4Mat <- ratioMat(S3C4) #14217if
SLC1Mat <- ratioMat(SLC1) #10076if
SLC2Mat <- ratioMat(SLC2) #12246if
SLC3Mat <- ratioMat(SLC3) #13122mo
SLC4Mat <- ratioMat(SLC4) #14217if

# Calculate ratio matrices Ara --------------------------------------------

AHC1Mat <- ratioMat(AHC1Mat) #Exp1
AHC2Mat <- ratioMat(AHC2Mat) #Exp2
AVCMat <- ratioMat(AVCMat) #Exp2
ALC1Mat <- ratioMat(ALC1Mat) #Exp1
ALC3Mat <- ratioMat(ALC3Mat) #Exp1
ALC5Mat <- ratioMat(ALC5Mat) #Exp1
ANCMat <- ratioMat(ANCMat) #Exp2

# Create list of ratio matrices --------------------------------------------

RMlist <- list(AHC1Mat, ALC1Mat, ALC3Mat, ALC5Mat, AVCMat, AHC2Mat, ANCMat, SHC1Mat, S3C1Mat, SLC1Mat, SHC2Mat, S3C2Mat, SLC2Mat, SHC3Mat, S3C3Mat, SLC3Mat, SHC4Mat, S3C4Mat, SLC4Mat)

# Calculate RV coefficients and create new matrix --------------------------------------------

names <- c("Ara_Exp1_HC", "Ara_Exp1_1d_LC", "Ara_Exp1_3d_LC", "Ara_Exp1_5d_LC", "Ara_Exp2_VC", "Ara_Exp2_HC", "Ara_Exp2_NC", "Syn_Exp1_HC", "Syn_Exp1_3h_LC", "Syn_Exp1_24h_LC", "Syn_Exp2_HC", "Syn_Exp2_3h_LC", "Syn_Exp2_24h_LC", "Syn_Exp3_HC", "Syn_Exp3_3h_LC", "Syn_Exp3_24h_LC", "Syn_Exp4_HC", "Syn_Exp4_3h_LC", "Syn_Exp4_24h_LC")

rv_matrix <- matrix(nrow = length(RMlist), ncol = length(RMlist))

colnames(rv_matrix) <- names
rownames(rv_matrix) <- names

for(i in 1:length(RMlist)) {
  tmp1 <- vector()
  for(j in 1:length(RMlist)) {
    tmp <- coeffRV(RMlist[[i]], RMlist[[j]])$rv
    tmp1 <- c(tmp1, tmp)
  }
  rv_matrix[i,] <- tmp1
}

fix(rv_matrix)

write.csv(rv_matrix,'rv_matrix.csv')

# Cluster matrix --------------------------------------------

par(oma=c(5,1,1,1))
heatmap(rv_matrix)

for (i in 2:10){
  clu <- pam(rv_matrix, k=i)
  cat(i," avg: ", round(clu$silinfo$avg.width,2))
  cat(" clus: ", round(clu$silinfo$clus.avg.widths,2),"\n")
}

#2  avg:  0.38 clus:  0.78 0.28 
#3  avg:  0.41 clus:  0.77 0.38 0.29 
#4  avg:  0.42 clus:  0.77 0.66 0.16 0.29 
####5  avg:  0.53 clus:  0.77 0.64 0.35 0.38 0.59 
#6  avg:  0.46 clus:  0.77 0.64 0.34 0.05 0.27 0.57 
#7  avg:  0.48 clus:  0.76 0.61 0 0.59 0.36 0.12 0.55 
#8  avg:  0.47 clus:  0.76 0.61 0 0.59 0.28 0 0.35 0.45 
#9  avg:  0.45 clus:  0.76 0 0.77 0 0.59 0.28 0 0.35 0.45 
#10  avg:  0.42 clus:  0 0.76 0 0.77 0 0.59 0.28 0 0.35 0.45

#Run pam
clu <- pam(rv_matrix, k=5)

color <- function(n = 50, low.col = 0.45, high.col=1, saturation = 1) { 
  if (n < 2) stop("n must be greater than 2")
  n1 <- n%/%2
  n2 <- n - n1
  c(hsv(low.col, saturation, seq(1,0,length=n1)), hsv(high.col, saturation, seq(0,1,length=n2))) 
}

pamcolors <- sample(rainbow(5)) 
pamcolors <- pamcolors[as.vector(clu$clustering)]

par(oma=c(5,2,2,5))
heatmap(rv_matrix, col=color(), scale="none", Rowv= NA, Colv = NA, RowSideColors=pamcolors, ColSideColors=pamcolors, main="PAM clustering of RV coefficients between metabolite ratio matrices (k=5)")
