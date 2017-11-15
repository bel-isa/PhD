(setwd("Y:/Isabel Orf/R-3.1.2")

library(FactoMineR)
library(ape)
 
# Load data ---------------------------------------------------------------
 
Ara <- read.delim("Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-09-01_Results/Ara.txt", sep="\t")
Syn <- read.delim("Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/2015-09-01_Results/Syn.txt", sep="\t")
 
Ara <- Ara[,-1]
Syn <- Syn[,-1]
 
Metabs <- read.delim("Y:/Isabel Orf/Projects/Ara_Syn/Metabolome/2015/15-03-09_vergleich-arasyn.txt")

rownames(Ara) <- Metabs[,1]
rownames(Syn) <- Metabs[,1]


# Split data into smaller sets --------------------------------------------

Syn_HC <- as.matrix(Syn[,1:22])
Syn_LC <- as.matrix(Syn[,23:ncol(Syn)])
 
Ara_HC <- as.matrix(cbind(Ara[,1:6], Ara[,109:ncol(Ara)]))
Ara_LC <- as.matrix(Ara[,7:108])
 

# Calculate ratio matrices ------------------------------------------------
 
SHC <- apply(Syn_HC,1,function(x) mean(x))

SHCMat <- matrix(nrow=30,ncol=30)

for (i in 1:length(SHC)){
  for (j in 1:length(SHC)){
    SHCMat[i,j] <- (SHC[i]/SHC[j])
  }
} 

SLC3 <- apply(Syn_LC[,1:20],1,function(x) mean(x))

SLC3Mat <- matrix(nrow=30,ncol=30)

for (i in 1:length(SLC3)){
  for (j in 1:length(SLC3)){
    SLC3Mat[i,j] <- (SLC3[i]/SLC3[j])
  }
}

SLC <- apply(Syn_LC[,21:ncol(Syn_LC)],1,function(x) mean(x))

SLCMat <- matrix(nrow=30,ncol=30)

for (i in 1:length(SLC)){
  for (j in 1:length(SLC)){
    SLCMat[i,j] <- (SLC[i]/SLC[j])
  }
}

AHC1 <- apply(Ara_HC[,1:3],1,function(x) mean(x))

AHC1Mat <- matrix(nrow=30,ncol=30)

for (i in 1:length(AHC1)){
  for (j in 1:length(AHC1)){
    AHC1Mat[i,j] <- (AHC1[i]/AHC1[j])
  }
} 

AHC2 <- apply(Ara_HC[,7:48],1,function(x) mean(x))

AHC2Mat <- matrix(nrow=30,ncol=30)

for (i in 1:length(AHC2)){
  for (j in 1:length(AHC2)){
    AHC2Mat[i,j] <- (AHC2[i]/AHC2[j])
  }
}

AVC <- apply(Ara_HC[,49:ncol(Ara_HC)],1,function(x) mean(x))

AVCMat <- matrix(nrow=30,ncol=30)

for (i in 1:length(AVC)){
  for (j in 1:length(AVC)){
    AVCMat[i,j] <- (AVC[i]/AVC[j])
  }
}

ALC1 <- apply(Ara_LC[,1:3],1,function(x) mean(x))

ALC1Mat <- matrix(nrow=30,ncol=30)

for (i in 1:length(ALC1)){
  for (j in 1:length(ALC1)){
    ALC1Mat[i,j] <- (ALC1[i]/ALC1[j])
  }
}

ALC3 <- apply(Ara_LC[,7:9],1,function(x) mean(x))

ALC3Mat <- matrix(nrow=30,ncol=30)

for (i in 1:length(ALC3)){
  for (j in 1:length(ALC3)){
    ALC3Mat[i,j] <- (ALC3[i]/ALC3[j])
  }
}

ALC5 <- apply(Ara_LC[,13:15],1,function(x) mean(x))

ALC5Mat <- matrix(nrow=30,ncol=30)

for (i in 1:length(ALC5)){
  for (j in 1:length(ALC5)){
    ALC5Mat[i,j] <- (ALC5[i]/ALC5[j])
  }
}

ANC <- apply(Ara_LC[,19:ncol(Ara_LC)],1,function(x) mean(x))

ANCMat <- matrix(nrow=30,ncol=30)

for (i in 1:length(ANC)){
  for (j in 1:length(ANC)){
    ANCMat[i,j] <- (ANC[i]/ANC[j])
  }
}

# RV coefficient ----------------------------------------------------------

#HC over LC Synechocystis

RV1 <- coeffRV(SHCMat,SLCMat) #0.4955128 pvalue 5.530661e-05
RV2 <- coeffRV(SHCMat,SLC3Mat) #0.7632532 pvalue 4.605497e-07

#HC over LC Arabidopsis

#Exp1
RV3 <- coeffRV(AHC1Mat,ALC1Mat) #0.5024054 pvalue 0.001666713
RV4 <- coeffRV(AHC1Mat,ALC3Mat) #0.6282526 pvalue 0.0007472327
RV5 <- coeffRV(AHC1Mat,ALC5Mat) #0.635488 pvalue 0.0007914657

#Exp2
RV6 <- coeffRV(AHC2Mat,ANCMat) #0.9551975 pvalue 0.003350622
RV7 <- coeffRV(AVCMat,ANCMat) #0.4205908 pvalue 0.0005411318

#Syn over Ara HC

#Ara Exp1
RV8 <- coeffRV(SHCMat,AHC1Mat) #0.03250871 pvalue 0.349052

#Ara Exp2
RV9 <- coeffRV(SHCMat,AHC2Mat) #0.003883954 pvalue 0.7534436
RV10 <- coeffRV(SHCMat,AVCMat) #0.002129556 pvalue 0.7790611

#Syn over Ara LC

#Syn 3h LC, Ara Exp1
RV11 <- coeffRV(SLC3Mat,ALC1Mat) #0.0140444 pvalue 0.5653941
RV12 <- coeffRV(SLC3Mat,ALC3Mat) #0.01690825 pvalue 0.535832
RV13 <- coeffRV(SLC3Mat,ALC5Mat) #0.03415384 pvalue 0.3523066

#Syn 3h LC, Ara Exp2
RV14 <- coeffRV(SLC3Mat,ANCMat) #0.002236374 pvalue 0.8357884

#Syn 24h LC, Ara Exp1
RV15 <- coeffRV(SLCMat,ALC1Mat) #0.01400645 pvalue 0.5670864
RV16 <- coeffRV(SLCMat,ALC3Mat) #0.01182259 pvalue 0.6129212
RV17 <- coeffRV(SLCMat,ALC5Mat) #0.02281387 pvalue 0.4698386

#Syn 24h LC, Ara Exp2
RV18 <- coeffRV(SLCMat,ANCMat) #0.01240679 pvalue 0.6577077

# Mantel Test -------------------------------------------------------------

#HC over LC Synechocystis
mantel.test(SHCMat,SLCMat, graph = TRUE,
            main = "Mantel test",
            xlab = "z-statistic", ylab = "Density",
            sub = "Synechocystis HC vs. Synechocystis 24h LC")
#z statistic: 935.752

mantel.test(SHCMat,SLC3Mat, graph = TRUE,
            main = "Mantel test",
            xlab = "z-statistic", ylab = "Density",
            sub = "Synechocystis HC vs. Synechocystis 3h LC")
#z statistic: 753.2473

#HC over LC Arabidopsis
mantel.test(AHC1Mat,ALC1Mat, graph = TRUE,
            main = "Mantel test",
            xlab = "z-statistic", ylab = "Density",
            sub = "Arabidopsis HC vs. Arabidopsis 1d LC (Exp1)")
#z statistic: 19972.25

mantel.test(AHC1Mat,ALC3Mat, graph = TRUE,
            main = "Mantel test",
            xlab = "z-statistic", ylab = "Density",
            sub = "Arabidopsis HC vs. Arabidopsis 3d LC (Exp1)")
#z statistic: 20664.59

mantel.test(AHC1Mat,ALC5Mat, graph = TRUE,
            main = "Mantel test",
            xlab = "z-statistic", ylab = "Density",
            sub = "Arabidopsis HC vs. Arabidopsis 5d LC (Exp1)")
#z statistic: 17936.55

mantel.test(AHC2Mat,ANCMat, graph = TRUE,
            main = "Mantel test",
            xlab = "z-statistic", ylab = "Density",
            sub = "Arabidopsis HC vs. Arabidopsis LC (Exp2)")
#z statistic: 16045

mantel.test(AVCMat,ANCMat, graph = TRUE,
            main = "Mantel test",
            xlab = "z-statistic", ylab = "Density",
            sub = "Arabidopsis VC vs. Arabidopsis LC (Exp2)")
#z statistic: 12439.9

#Syn over Ara HC
mantel.test(SHCMat,AHC1Mat, graph = TRUE,
            main = "Mantel test",
            xlab = "z-statistic", ylab = "Density",
            sub = "Synechocystis HC vs. Arabidopsis HC (Exp1)")
#z statistic: 1298.52

mantel.test(SHCMat,AHC2Mat, graph = TRUE,
            main = "Mantel test",
            xlab = "z-statistic", ylab = "Density",
            sub = "Synechocystis HC vs. Arabidopsis HC (Exp2)")
#z statistic: 800.8158

mantel.test(SHCMat,AVCMat, graph = TRUE,
            main = "Mantel test",
            xlab = "z-statistic", ylab = "Density",
            sub = "Synechocystis HC vs. Arabidopsis VC (Exp2)")
#z statistic: 762.4464

#Syn over Ara LC
mantel.test(SLC3Mat,ALC1Mat, graph = TRUE,
            main = "Mantel test",
            xlab = "z-statistic", ylab = "Density",
            sub = "Synechocystis 3h LC vs. Arabidopsis 1d LC (Exp1)")
#z statistic: 1223.841

mantel.test(SLC3Mat,ALC3Mat, graph = TRUE,
            main = "Mantel test",
            xlab = "z-statistic", ylab = "Density",
            sub = "Synechocystis 3h LC vs. Arabidopsis 3d LC (Exp1)")
#z statistic: 1160.111

mantel.test(SLC3Mat,ALC5Mat, graph = TRUE,
            main = "Mantel test",
            xlab = "z-statistic", ylab = "Density",
            sub = "Synechocystis 3h LC vs. Arabidopsis 5d LC (Exp1)")
#z statistic: 1058.23

mantel.test(SLC3Mat,ANCMat, graph = TRUE,
            main = "Mantel test",
            xlab = "z-statistic", ylab = "Density",
            sub = "Synechocystis 3h LC vs. Arabidopsis LC (Exp2)")
#z statistic: 823.5141

mantel.test(SLCMat,ALC1Mat, graph = TRUE,
            main = "Mantel test",
            xlab = "z-statistic", ylab = "Density",
            sub = "Synechocystis 24h LC vs. Arabidopsis 1d LC (Exp1)")
#z statistic: 1393.808

mantel.test(SLCMat,ALC3Mat, graph = TRUE,
            main = "Mantel test",
            xlab = "z-statistic", ylab = "Density",
            sub = "Synechocystis 24h LC vs. Arabidopsis 3d LC (Exp1)")
#z statistic: 1274.243

mantel.test(SLCMat,ALC5Mat, graph = TRUE,
            main = "Mantel test",
            xlab = "z-statistic", ylab = "Density",
            sub = "Synechocystis 24h LC vs. Arabidopsis 5d LC (Exp1)")
#z statistic: 1169.332

mantel.test(SLCMat,ANCMat, graph = TRUE,
            main = "Mantel test",
            xlab = "z-statistic", ylab = "Density",
            sub = "Synechocystis 24h LC vs. Arabidopsis LC (Exp2)")
#z statistic: 896.4015