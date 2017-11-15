#######################################################################################
# Script to calculate averages, standard error and perform t-tests between conditions #
#######################################################################################

setwd("C:/Users/orf.MPIMP-GOLM/ownCloud/Project/002_Pseudomonas_Omics/Isabel/LC-MS")

# Data

Metabs <- read.delim("Secondary_metabs_all.txt", sep="\t", header=TRUE, na.strings="NA")
tmp <- Metabs[,-1]
tmp1 <- as.matrix(tmp)
row.names(tmp1) <- Metabs[,1]
Metabs <- tmp1

Metabs_t <- t(Metabs)

# Group data into genotypes

Col0 <- Metabs_t[c(grep("X0", rownames(Metabs_t))),]
C24 <- Metabs_t[c(grep("X4", rownames(Metabs_t))),]

# Group data into treatments

Col0_M <- Col0[c(grep("m", rownames(Col0))),]
Col0_P <- Col0[c(grep("p", rownames(Col0))),]

C24_M <- C24[c(grep("m", rownames(C24))),]
C24_P <- C24[c(grep("p", rownames(C24))),]

# Group data into time points

A <- Col0_M[c(grep("m2_", rownames(Col0_M))),]
B <- Col0_M[c(grep("m4_", rownames(Col0_M))),]
C <- Col0_M[c(grep("m8_", rownames(Col0_M))),]
D <- Col0_M[c(grep("m16_", rownames(Col0_M))),]
E <- Col0_M[c(grep("m24_", rownames(Col0_M))),]
G <- Col0_M[c(grep("m48_", rownames(Col0_M))),]

H <- Col0_P[c(grep("p2_", rownames(Col0_P))),]
I <- Col0_P[c(grep("p4_", rownames(Col0_P))),]
J <- Col0_P[c(grep("p8_", rownames(Col0_P))),]
K <- Col0_P[c(grep("p16_", rownames(Col0_P))),]
L <- Col0_P[c(grep("p24_", rownames(Col0_P))),]
M <- Col0_P[c(grep("p48_", rownames(Col0_P))),]

N <- C24_M[c(grep("m2_", rownames(C24_M))),]
O <- C24_M[c(grep("m4_", rownames(C24_M))),]
P <- C24_M[c(grep("m8_", rownames(C24_M))),]
Q <- C24_M[c(grep("m16_", rownames(C24_M))),]
R <- C24_M[c(grep("m24_", rownames(C24_M))),]
S <- C24_M[c(grep("m48_", rownames(C24_M))),]

U <- C24_P[c(grep("p2_", rownames(C24_P))),]
V <- C24_P[c(grep("p4_", rownames(C24_P))),]
W <- C24_P[c(grep("p8_", rownames(C24_P))),]
X <- C24_P[c(grep("p16_", rownames(C24_P))),]
Y <- C24_P[c(grep("p24_", rownames(C24_P))),]
Z <- C24_P[c(grep("p48_", rownames(C24_P))),]

# Calculate averages
myav <- function(x) {if(length(x[!is.na(x)])>1){mean(x, na.rm =T)} 
  else(return("NA"))
}

names <- c("Av 2h", "Av 4h", "Av 8h", "Av 16h", "Av 24h", "Av 48h")

AvCol0_M <- cbind(apply(A, 2, myav), apply(B, 2, myav), apply(C, 2, myav), apply(D, 2, myav), apply(E, 2, myav), apply(G, 2, myav))
colnames(AvCol0_M) <- names

AvCol0_P <- cbind(apply(H, 2, myav), apply(I, 2, myav), apply(J, 2, myav), apply(K, 2, myav), apply(L, 2, myav), apply(M, 2, myav))
colnames(AvCol0_P) <- names

AvC24_M <- cbind(apply(N, 2, myav), apply(O, 2, myav), apply(P, 2, myav), apply(Q, 2, myav), apply(R, 2, myav), apply(S, 2, myav))
colnames(AvC24_M) <- names

AvC24_P <- cbind(apply(U, 2, myav), apply(V, 2, myav), apply(W, 2, myav), apply(X, 2, myav), apply(Y, 2, myav), apply(Z, 2, myav))
colnames(AvC24_P) <- names

# Calculate standard error

serr <- function(x) {if(length(x[!is.na(x)])>1){sd(x, na.rm =T)/sqrt(length(x))}
  else(return("NA"))
}

A1 <- apply(A, 2, serr)
B1 <- apply(B, 2, serr)
C1 <- apply(C, 2, serr)
D1 <- apply(D, 2, serr)
E1 <- apply(E, 2, serr)
G1 <- apply(G, 2, serr)

H1 <- apply(H, 2, serr)
I1 <- apply(I, 2, serr)
J1 <- apply(J, 2, serr)
K1 <- apply(K, 2, serr)
L1 <- apply(L, 2, serr)
M1 <- apply(M, 2, serr)

N1 <- apply(N, 2, serr)
O1 <- apply(O, 2, serr)
P1 <- apply(P, 2, serr)
Q1 <- apply(Q, 2, serr)
R1 <- apply(R, 2, serr)
S1 <- apply(S, 2, serr)

U1 <- apply(U, 2, serr)
V1 <- apply(V, 2, serr)
W1 <- apply(W, 2, serr)
X1 <- apply(X, 2, serr)
Y1 <- apply(Y, 2, serr)
Z1 <- apply(Z, 2, serr)

names1 <- c("SE 2h", "SE 4h", "SE 8h", "SE 16h", "SE 24h", "SE 48h")

SECol0_M <- cbind(A1, B1, C1, D1, E1, G1)
colnames(SECol0_M) <- names1

SECol0_P <- cbind(H1, I1, J1, K1, L1, M1)
colnames(SECol0_P) <- names1

SEC24_M <- cbind(N1, O1, P1, Q1, R1, S1)
colnames(SEC24_M) <- names1

SEC24_P <- cbind(U1, V1, W1, X1, Y1, Z1)
colnames(SEC24_P) <- names1

# Combine table and save file

write.table(cbind(AvCol0_M, SECol0_M), "Col0_M.txt", sep="\t")
write.table(cbind(AvCol0_P, SECol0_P), "Col0_P.txt", sep="\t")

write.table(cbind(AvC24_M, SEC24_M), "C24_M.txt", sep="\t")
write.table(cbind(AvC24_P, SEC24_P), "C24_P.txt", sep="\t")

# T-Test

mytest <- function(x ,y, n){ #Function on two matrices x and y, comparing n
  metabX <- x[,n] #n is a column in x
  metabY <- y[,n] #n is a column in y
  if(sum(!is.na(metabX)) > 1 & sum(!is.na(metabY)) > 1){ #Only run the function if metabX and metabY have more than 1 value
    P <- t.test(metabX, metabY)$p.value #This is the function
    }else{ #Otherwise return NA
    P <- NA
  }
  return(P)
}

AH <- vector() #Create vector to be filled with result of mytest
for(i in 1:ncol(A)){ #For-loop... for every column in A
  AH <- c(AH, mytest(A, H, i)) #Write result of function mytest in vector
}
#AH <- p.adjust(AH, method= c("fdr")) #Benjamini-Hochberg correction for multiple testing

BI <- vector()
for(i in 1:ncol(B)){
  BI <- c(BI, mytest(B, I, i))
}
#BI <- p.adjust(BI, method= c("fdr"))

CJ <- vector()
for(i in 1:ncol(C)){
  CJ <- c(CJ, mytest(C, J, i))
}
#CJ <- p.adjust(CJ, method= c("fdr"))
               
DK <- vector()
for(i in 1:ncol(D)){
  DK <- c(DK, mytest(D, K, i))
}
#DK <- p.adjust(DK, method= c("fdr"))

EL <- vector()
for(i in 1:ncol(E)){
  EL <- c(EL, mytest(E, L, i))
}
#EL <- p.adjust(EL, method= c("fdr"))

GM <- vector()
for(i in 1:ncol(G)){
  GM <- c(GM, mytest(G, M, i))
}
#GM <- p.adjust(GM, method= c("fdr"))

NU <- vector()
for(i in 1:ncol(N)){
  NU <- c(NU, mytest(N, U, i))
}
#NU <- p.adjust(NU, method= c("fdr"))

OV <- vector()
for(i in 1:ncol(O)){
  OV <- c(OV, mytest(O, V, i))
}
#OV <- p.adjust(OV, method= c("fdr"))

PW <- vector()
for(i in 1:ncol(P)){
  PW <- c(PW, mytest(P, W, i))
}
#PW <- p.adjust(PW, method= c("fdr"))

QX <- vector()
for(i in 1:ncol(Q)){
  QX <- c(QX, mytest(Q, X, i))
}
#QX <- p.adjust(QX, method= c("fdr"))

RY <- vector()
for(i in 1:ncol(R)){
  RY <- c(RY, mytest(R, Y, i))
}
#RY <- p.adjust(RY, method= c("fdr"))

SZ <- vector()
for(i in 1:ncol(S)){
  SZ <- c(SZ, mytest(S, Z, i))
}
#SZ <- p.adjust(SZ, method= c("fdr"))

# For Col-0

names2 <- c("p 2h", "p 4h", "p 8h", "p 16h", "p 24h", "p 48h")
testCol0 <- cbind(AH, BI, CJ, DK, EL, GM)
colnames(testCol0) <- names2
rownames(testCol0) <- colnames(A)

write.table(testCol0, "T-Test_Col-0.txt", sep="\t")

# For C24

testC24 <- cbind(NU, OV, PW, QX, RY, SZ)
colnames(testC24) <- names2
rownames(testC24) <- colnames(A)

write.table(testC24, "T-Test_C24.txt", sep="\t")

# Create table of all results

write.table(cbind(AvCol0_M, SECol0_M, AvCol0_P, SECol0_P, AvC24_M, SEC24_M, AvC24_P, SEC24_P, testCol0, testC24), "Secondary_metabolites_R_analysis.txt", sep="\t")
