# Set directory and load packages -----------------------------------------

setwd("Y:/Isabel Orf/Projects/01_Ara_Syn/Metabolome/2015/2015-10-05_Results")

library(abind)
library(missForest)
library(cluster)
library(MASS)
library(Hmisc)
library(FactoMineR)
library(gplots)

# Read datasets -----------------------------------------

SynExp1 <- read.delim("SynExp1.txt", header=T)
rownames(SynExp1) <- SynExp1[,1]
SynExp1 <- SynExp1[,-1]
SynExp2 <- read.delim("SynExp2.txt", header=T)
rownames(SynExp2) <- SynExp2[,1]
SynExp2 <- SynExp2[,-1]
SynExp3 <- read.delim("SynExp3.txt", header=T)
rownames(SynExp3) <- SynExp3[,1]
SynExp3 <- SynExp3[,-1]
SynExp4 <- read.delim("SynExp4.txt", header=T)
rownames(SynExp4) <- SynExp4[,1]
SynExp4 <- SynExp4[,-1]

#Exp1 and Exp4 contain samples for the ccmM mutant, in case the mutant should be included in the comparison as well
SynExp1ccmM <- SynExp1[,13:ncol(SynExp1)]
SynExp1 <- SynExp1[,1:12]

SynExp4ccmM <- SynExp4[,10:ncol(SynExp4)]
SynExp4 <- SynExp4[,1:9]

AraExp1day <- read.delim("AraExp1day.txt", header=T)
rownames(AraExp1day) <- AraExp1day[,1]
AraExp1day <- AraExp1day[,-1]
AraExp1night <- read.delim("AraExp1night.txt", header=T)
rownames(AraExp1night) <- AraExp1night[,1]
AraExp1night <- AraExp1night[,-1]
AraExp2VC <- read.delim("AraExp2_VC.txt", header=T)
rownames(AraExp2VC) <- AraExp2VC[,1]
AraExp2VC <- AraExp2VC[,-1]
AraExp2HC <- read.delim("AraExp2_HC.txt", header=T)
rownames(AraExp2HC) <- AraExp2HC[,1]
AraExp2HC <- AraExp2HC[,-1]

# Count NA function -----------------------------------------

metab.nas <- function(x) {
  tmp <- as.numeric(apply(x, 1, function(z) sum(is.na(z))))
  whichmetabs <- rownames(x)[which(tmp > ncol(x)/2)]
  return(whichmetabs)
}

sum.nas <- function(x) {
  tmp <- as.numeric(apply(x, 1, function(z) sum(is.na(z))))
  whichmetabs <- which(tmp > ncol(x)/2)
  return(whichmetabs)
}

# Output = Metabolites that contain more than 50% NAs per Experiment 

# Count NAs -----------------------------------------

metab.nas(SynExp1)
metab.nas(SynExp2)
metab.nas(SynExp3)
metab.nas(SynExp4)
metab.nas(AraExp1day)
metab.nas(AraExp1night)
metab.nas(AraExp2VC)
metab.nas(AraExp2HC)

# Remove Metabolites with too many NAs function -----------------------------------------

remove <- function(x) {
  tmp <- as.matrix(x[-c(1,2,4,5,7,9,10:13,16,19,22,24,26,30),])
  return(tmp)
}

# Remove Metabolites with too many NAs -----------------------------------------

SynExp1.1 <- remove(SynExp1)
SynExp3.1 <- remove(SynExp3)
SynExp4.1 <- remove(SynExp4)

AraExp1day.1 <- remove(AraExp1day)
AraExp1night.1 <- remove(AraExp1night)
AraExp2VC.1 <- remove(AraExp2VC)
AraExp2HC.1 <- remove(AraExp2HC)

# Percentage of missing values funtion -------------------------------

percent.nas <- function(x) {
  tmp <- sum(is.na(x))
  percentvalues <- tmp/length(x)*100
  return(percentvalues)
}

# How many missing values left? -------------------------------

percent.nas(SynExp1.1)
percent.nas(SynExp3.1)
percent.nas(SynExp4.1)
percent.nas(AraExp1day.1)
percent.nas(AraExp1night.1)
percent.nas(AraExp2VC.1)
percent.nas(AraExp2HC.1)

# Impute missing values using random forest function -------------------------------

jungle <- function(x){
  tmp <- abind(replicate(100,missForest(x)$ximp,simplify=F),along=3)
  tmp1 <- apply(tmp, c(1,2), mean)
  return(tmp1)
}

#The replacement is performed 100 times and the average over all replacements is used for further calculations

# Impute missing values -------------------------------

SynExp1.2 <- jungle(SynExp1.1)
SynExp3.2 <- jungle(SynExp3.1)
SynExp4.2 <- jungle(SynExp4.1)

AraExp1day.2 <- jungle(AraExp1day.1)
AraExp1night.2 <- jungle(AraExp1night.1)
AraExp2VC.2 <- jungle(AraExp2VC.1)
AraExp2HC.2 <- jungle(AraExp2HC.1)

# Median center function -----------------------------------------

median.center <- function(x) {
  t(apply(x, 1, function(x) x/median(x)))
}

#Function to normalize the datasets

# Median center datasets -----------------------------------------

SynExp1.3 <- median.center(SynExp1.2)
SynExp3.3 <- median.center(SynExp3.2)
SynExp4.3 <- median.center(SynExp4.2)

AraExp1day.3 <- median.center(AraExp1day.2)
AraExp1night.3 <- median.center(AraExp1night.2)
AraExp2VC.3 <- median.center(AraExp2VC.2)
AraExp2HC.3 <- median.center(AraExp2HC.2)

#Add experimental identifiers to the dataset for later separation

a <- rep(1, ncol(SynExp1.3))
b <- rep(2, ncol(SynExp3.3))
c <- rep(3, ncol(SynExp4.3))
d <- rep(4, ncol(AraExp1day.3))
e <- rep(5, ncol(AraExp1night.3))
f <- rep(6, ncol(AraExp2VC.3))
g <- rep(7, ncol(AraExp2HC.3))

exp.ident <- c(a,b,c,d,e,f,g)

# Combine datasets -----------------------------------------

Full <- cbind(SynExp1.3, SynExp3.3, SynExp4.3, AraExp1day.3, AraExp1night.3, AraExp2VC.3, AraExp2HC.3)

# Classify samples (split dataset again) --------------------------------------------------------

Full.1 <- rbind(Full, exp.ident)

exp1 <- subset(t(Full.1), exp.ident==1)
exp2 <- subset(t(Full.1), exp.ident==2)
exp3 <- subset(t(Full.1), exp.ident==3)
exp4 <- subset(t(Full.1), exp.ident==4)
exp5 <- subset(t(Full.1), exp.ident==5)
exp6 <- subset(t(Full.1), exp.ident==6)
exp7 <- subset(t(Full.1), exp.ident==7)

# Create condition vectors for the experiments and combine with data

cond <- c(rep("SHC",4), rep("S3LC",4), rep("S24LC",4))
cond.ident <- c(rep(1,4), rep(2,4), rep(3,4))
exp1 <- cbind(as.data.frame(exp1), cond, cond.ident)

cond <- c(rep("SHC",9), rep("S3LC",9), rep("S24LC",7))
cond.ident <- c(rep(1,9), rep(2,9), rep(3,7))
exp2 <- cbind(as.data.frame(exp2), cond, cond.ident)

cond <- c(rep("SHC",3), rep("S3LC",3), rep("S24LC",3))
cond.ident <- c(rep(1,3), rep(2,3), rep(3,3))
exp3 <- cbind(as.data.frame(exp3), cond, cond.ident)

cond <- c(rep("AHCd",3), rep("A1LCd",3), rep("A3LCd",3), rep("A5LCd",3))
cond.ident <- c(rep(4,3), rep(5,3), rep(6,3), rep(7,3))
exp4 <- cbind(as.data.frame(exp4), cond, cond.ident)

cond <- c(rep("AHCn",3), rep("A1LCn",3), rep("A3LCn",3), rep("A5LCn",3))
cond.ident <- c(rep(8,3), rep(9,3), rep(10,3), rep(11,3))
exp5 <- cbind(as.data.frame(exp5), cond, cond.ident)

exp6NC1 <- exp6[c(grep("NC1", rownames(exp6))),]
exp6NC2 <- exp6[c(grep("NC2", rownames(exp6))),]
exp6NC3 <- exp6[c(grep("NC3", rownames(exp6))),]
exp6NC4 <- exp6[c(grep("NC4", rownames(exp6))),]
exp6NC5 <- exp6[c(grep("NC5", rownames(exp6))),]
exp6NC6 <- exp6[c(grep("NC6", rownames(exp6))),]
exp6NC7 <- exp6[c(grep("NC7", rownames(exp6))),]
exp6NC <- rbind(exp6NC1, exp6NC2, exp6NC3, exp6NC4, exp6NC5, exp6NC6, exp6NC7)
cond <- c(rep("A1NC", 6), rep("A2NC", 6), rep("A3NC", 6), rep("A4NC", 6), rep("A5NC", 6), rep("A6NC", 6), rep("A7NC", 6))
cond.ident <- c(rep(12, 6), rep(13, 6), rep(14, 6), rep(15, 6), rep(16, 6), rep(17, 6), rep(18, 6))
exp6NC <- cbind(as.data.frame(exp6NC), cond, cond.ident)

exp7NC1 <- exp7[c(grep("NC1", rownames(exp7))),]
exp7NC2 <- exp7[c(grep("NC2", rownames(exp7))),]
exp7NC3 <- exp7[c(grep("NC3", rownames(exp7))),]
exp7NC4 <- exp7[c(grep("NC4", rownames(exp7))),]
exp7NC5 <- exp7[c(grep("NC5", rownames(exp7))),]
exp7NC6 <- exp7[c(grep("NC6", rownames(exp7))),]
exp7NC7 <- exp7[c(grep("NC7", rownames(exp7))),]
exp7NC <- rbind(exp7NC1, exp7NC2, exp7NC3, exp7NC4, exp7NC5, exp7NC6, exp7NC7)
cond <- c(rep("A1NC", 6), rep("A2NC", 6), rep("A3NC", 6), rep("A4NC", 6), rep("A5NC", 6), rep("A6NC", 6), rep("A7NC", 6))
cond.ident <- c(rep(12, 6), rep(13, 6), rep(14, 6), rep(15, 6), rep(16, 6), rep(17, 6), rep(18, 6))
exp7NC <- cbind(as.data.frame(exp7NC), cond, cond.ident)

exp7HC1 <- exp7[c(grep("HC1", rownames(exp7))),]
exp7HC2 <- exp7[c(grep("HC2", rownames(exp7))),]
exp7HC3 <- exp7[c(grep("HC3", rownames(exp7))),]
exp7HC4 <- exp7[c(grep("HC4", rownames(exp7))),]
exp7HC5 <- exp7[c(grep("HC5", rownames(exp7))),]
exp7HC6 <- exp7[c(grep("HC6", rownames(exp7))),]
exp7HC7 <- exp7[c(grep("HC7", rownames(exp7))),]
exp7HC <- rbind(exp7HC1, exp7HC2, exp7HC3, exp7HC4, exp7HC5, exp7HC6, exp7HC7)
cond <- c(rep("A1HC", 6), rep("A2HC", 6), rep("A3HC", 6), rep("A4HC", 6), rep("A5HC", 6), rep("A6HC", 6), rep("A7HC", 6))
cond.ident <- c(rep(19, 6), rep(20, 6), rep(21, 6), rep(22, 6), rep(23, 6), rep(24, 6), rep(25, 6))
exp7HC <- cbind(as.data.frame(exp7HC), cond, cond.ident)

exp6VC1 <- exp6[c(grep("VC1", rownames(exp6))),]
exp6VC2 <- exp6[c(grep("VC2", rownames(exp6))),]
exp6VC3 <- exp6[c(grep("VC3", rownames(exp6))),]
exp6VC4 <- exp6[c(grep("VC4", rownames(exp6))),]
exp6VC5 <- exp6[c(grep("VC5", rownames(exp6))),]
exp6VC6 <- exp6[c(grep("VC6", rownames(exp6))),]
exp6VC7 <- exp6[c(grep("VC7", rownames(exp6))),]
exp6VC <- rbind(exp6VC1, exp6VC2, exp6VC3, exp6VC4, exp6VC5, exp6VC6, exp6VC7)
cond <- c(rep("A1VC", 6), rep("A2VC", 6), rep("A3VC", 8), rep("A4VC", 6), rep("A5VC", 6), rep("A6VC", 6), rep("A7VC", 6))
cond.ident <- c(rep(26, 6), rep(27, 6), rep(28, 8), rep(29, 6), rep(30, 6), rep(31, 6), rep(32, 6))
exp6VC <- cbind(as.data.frame(exp6VC), cond, cond.ident)

Full.2 <- rbind(exp1, exp2, exp3, exp4, exp5, exp6VC, exp6NC, exp7HC, exp7NC)

exp1HC <- subset(exp1, exp1[,"cond.ident"]==1)
exp1HC <- exp1HC[,-c(15:ncol(exp1HC))]
exp1HC <- apply(exp1HC, 2, mean)
exp1LC1 <- subset(exp1, exp1[,"cond.ident"]==2)
exp1LC1 <- exp1LC1[,-c(15:ncol(exp1LC1))]
exp1LC1 <- apply(exp1LC1, 2, mean)
exp1LC2 <- subset(exp1, exp1[,"cond.ident"]==3)
exp1LC2 <- exp1LC2[,-c(15:ncol(exp1LC2))]
exp1LC2 <- apply(exp1LC2, 2, mean)

exp2HC <- subset(exp2, exp2[,"cond.ident"]==1)
exp2HC <- exp2HC[,-c(15:ncol(exp2HC))]
exp2HC <- apply(exp2HC, 2, mean)
exp2LC1 <- subset(exp2, exp2[,"cond.ident"]==2)
exp2LC1 <- exp2LC1[,-c(15:ncol(exp2LC1))]
exp2LC1 <- apply(exp2LC1, 2, mean)
exp2LC2 <- subset(exp2, exp2[,"cond.ident"]==3)
exp2LC2 <- exp2LC2[,-c(15:ncol(exp2LC2))]
exp2LC2 <- apply(exp2LC2, 2, mean)

exp3HC <- subset(exp3, exp3[,"cond.ident"]==1)
exp3HC <- exp3HC[,-c(15:ncol(exp3HC))]
exp3HC <- apply(exp3HC, 2, mean)
exp3LC1 <- subset(exp3, exp3[,"cond.ident"]==2)
exp3LC1 <- exp3LC1[,-c(15:ncol(exp3LC1))]
exp3LC1 <- apply(exp3LC1, 2, mean)
exp3LC2 <- subset(exp3, exp3[,"cond.ident"]==3)
exp3LC2 <- exp3LC2[,-c(15:ncol(exp3LC2))]
exp3LC2 <- apply(exp3LC2, 2, mean)

exp4HC <- subset(exp4, exp4[,"cond.ident"]==4)
exp4HC <- exp4HC[,-c(15:ncol(exp4HC))]
exp4HC <- apply(exp4HC, 2, mean)
exp4LC1 <- subset(exp4, exp4[,"cond.ident"]==5)
exp4LC1 <- exp4LC1[,-c(15:ncol(exp4LC1))]
exp4LC1 <- apply(exp4LC1, 2, mean)
exp4LC2 <- subset(exp4, exp4[,"cond.ident"]==6)
exp4LC2 <- exp4LC2[,-c(15:ncol(exp4LC2))]
exp4LC2 <- apply(exp4LC2, 2, mean)
exp4LC3 <- subset(exp4, exp4[,"cond.ident"]==7)
exp4LC3 <- exp4LC3[,-c(15:ncol(exp4LC3))]
exp4LC3 <- apply(exp4LC3, 2, mean)

exp5HC <- subset(exp5, exp5[,"cond.ident"]==8)
exp5HC <- exp5HC[,-c(15:ncol(exp5HC))]
exp5HC <- apply(exp5HC, 2, mean)
exp5LC1 <- subset(exp5, exp5[,"cond.ident"]==9)
exp5LC1 <- exp5LC1[,-c(15:ncol(exp5LC1))]
exp5LC1 <- apply(exp5LC1, 2, mean)
exp5LC2 <- subset(exp5, exp5[,"cond.ident"]==10)
exp5LC2 <- exp5LC2[,-c(15:ncol(exp5LC2))]
exp5LC2 <- apply(exp5LC2, 2, mean)
exp5LC3 <- subset(exp5, exp5[,"cond.ident"]==11)
exp5LC3 <- exp5LC3[,-c(15:ncol(exp5LC3))]
exp5LC3 <- apply(exp5LC3, 2, mean)

exp6VC1 <- exp6VC1[,-c(15:ncol(exp6VC1))]
exp6VC1 <- apply(exp6VC1, 2, mean)
exp6VC2 <- exp6VC2[,-c(15:ncol(exp6VC2))]
exp6VC2 <- apply(exp6VC2, 2, mean)
exp6VC3 <- exp6VC3[,-c(15:ncol(exp6VC3))]
exp6VC3 <- apply(exp6VC3, 2, mean)
exp6VC4 <- exp6VC4[,-c(15:ncol(exp6VC4))]
exp6VC4 <- apply(exp6VC4, 2, mean)
exp6VC5 <- exp6VC5[,-c(15:ncol(exp6VC5))]
exp6VC5 <- apply(exp6VC5, 2, mean)
exp6VC6 <- exp6VC6[,-c(15:ncol(exp6VC6))]
exp6VC6 <- apply(exp6VC6, 2, mean)
exp6VC7 <- exp6VC7[,-c(15:ncol(exp6VC7))]
exp6VC7 <- apply(exp6VC7, 2, mean)

exp6NC1 <- exp6NC1[,-c(15:ncol(exp6NC1))]
exp6NC1 <- apply(exp6NC1, 2, mean)
exp6NC2 <- exp6NC2[,-c(15:ncol(exp6NC2))]
exp6NC2 <- apply(exp6NC2, 2, mean)
exp6NC3 <- exp6NC3[,-c(15:ncol(exp6NC3))]
exp6NC3 <- apply(exp6NC3, 2, mean)
exp6NC4 <- exp6NC4[,-c(15:ncol(exp6NC4))]
exp6NC4 <- apply(exp6NC4, 2, mean)
exp6NC5 <- exp6NC5[,-c(15:ncol(exp6NC5))]
exp6NC5 <- apply(exp6NC5, 2, mean)
exp6NC6 <- exp6NC6[,-c(15:ncol(exp6NC6))]
exp6NC6 <- apply(exp6NC6, 2, mean)
exp6NC7 <- exp6NC7[,-c(15:ncol(exp6NC7))]
exp6NC7 <- apply(exp6NC7, 2, mean)

exp7HC1 <- exp7HC1[,-c(15:ncol(exp7HC1))]
exp7HC1 <- apply(exp7HC1, 2, mean)
exp7HC2 <- exp7HC2[,-c(15:ncol(exp7HC2))]
exp7HC2 <- apply(exp7HC2, 2, mean)
exp7HC3 <- exp7HC3[,-c(15:ncol(exp7HC3))]
exp7HC3 <- apply(exp7HC3, 2, mean)
exp7HC4 <- exp7HC4[,-c(15:ncol(exp7HC4))]
exp7HC4 <- apply(exp7HC4, 2, mean)
exp7HC5 <- exp7HC5[,-c(15:ncol(exp7HC5))]
exp7HC5 <- apply(exp7HC5, 2, mean)
exp7HC6 <- exp7HC6[,-c(15:ncol(exp7HC6))]
exp7HC6 <- apply(exp7HC6, 2, mean)
exp7HC7 <- exp7HC7[,-c(15:ncol(exp7HC7))]
exp7HC7 <- apply(exp7HC7, 2, mean)

exp7NC1 <- exp7NC1[,-c(15:ncol(exp7NC1))]
exp7NC1 <- apply(exp7NC1, 2, mean)
exp7NC2 <- exp7NC2[,-c(15:ncol(exp7NC2))]
exp7NC2 <- apply(exp7NC2, 2, mean)
exp7NC3 <- exp7NC3[,-c(15:ncol(exp7NC3))]
exp7NC3 <- apply(exp7NC3, 2, mean)
exp7NC4 <- exp7NC4[,-c(15:ncol(exp7NC4))]
exp7NC4 <- apply(exp7NC4, 2, mean)
exp7NC5 <- exp7NC5[,-c(15:ncol(exp7NC5))]
exp7NC5 <- apply(exp7NC5, 2, mean)
exp7NC6 <- exp7NC6[,-c(15:ncol(exp7NC6))]
exp7NC6 <- apply(exp7NC6, 2, mean)
exp7NC7 <- exp7NC7[,-c(15:ncol(exp7NC7))]
exp7NC7 <- apply(exp7NC7, 2, mean)


# CV function -----------------------------------------------------

calc.CV <- function(x) {
  tmp <- as.numeric(apply(x, 2, function(z) sd(z)/mean(z)))
  return(tmp)
}

#Calculate coefficient of variation

# CV HC and LC -----------------------------------------------------

call.cond <- Full.2[,"cond.ident"]

#Group samples into HC anc LC

SHC <- subset(Full.2, call.cond==1)
SLC <- subset(Full.2, call.cond==3)

AHC <- rbind(subset(Full.2, call.cond==4), subset(Full.2, call.cond==22), subset(Full.2, call.cond==23), subset(Full.2, call.cond==29), subset(Full.2, call.cond==30))
ALC <- rbind(subset(Full.2, call.cond==5), subset(Full.2, call.cond==6), subset(Full.2, call.cond==7), subset(Full.2, call.cond==15), subset(Full.2, call.cond==16))

#Calculate CV

CV.SHC <- calc.CV(SHC[,-c(15:17)])
CV.SLC <- calc.CV(SLC[,-c(15:17)])

CV.AHC <- calc.CV(AHC[,-c(15:17)])
CV.ALC <- calc.CV(ALC[,-c(15:17)])

CV.S <- c(CV.SLC,CV.SHC)
CV.A <- c(CV.ALC,CV.AHC)

HCLC <- c(rep(1,14),rep(2,14))

plot(CV.S, CV.A, col= c("red","black")[HCLC])

# Analyze CV correlation, optimize graphics for publication ---------------

par(mfrow = c(1,3))
plot(CV.S,CV.A, col= c("red","black")[HCLC], pch = c(15,17) [HCLC], xlab = "Synechocystis", ylab = "Arabidopsis", main= "Coefficient of variation", ylim=c(0,1.4), xlim=c(0,1.4))
metabnames <- rownames(Full)
leg.txt <- c("LC", "HC")
legend(list(x=0,y=1.3), legend=leg.txt, col= c("red","black"), pch = c(15,17))
text(CV.S, CV.A, metabnames, cex=0.5, pos=4, col="black")

reg <- lm(CV.A ~ CV.S-1)
abline(reg)
summary(reg)
impr <- rstandard(reg)
names(impr) <- c(metabnames,metabnames)
summary(impr)

#Test two different residual error cut-offs

data4 <- abs(impr) < 1

plot(CV.S[data4],CV.A[data4], ylim=c(0,1.4), xlim=c(0,1.4))


reg1 <- lm(CV.A[data4] ~ CV.S[data4]-1)
abline(reg1, col=c("red"))
summary(reg1)

data5 <- abs(impr) < 0.7

plot(CV.S[data5],CV.A[data5], ylim=c(0,1.4), xlim=c(0,1.4))


reg2 <- lm(CV.A[data5] ~ CV.S[data5]-1)
abline(reg2)
summary(reg2)

#Final plot

par(mfrow = c(1,1))
plot(CV.S,CV.A, col= c("red","black")[HCLC], pch = c(15,17) [HCLC], xlab = "Synechocystis", ylab = "Arabidopsis", main= "Coefficient of variation", ylim=c(0,1.4), xlim=c(0,1.4))
metabnames <- rownames(Full)
leg.txt <- c("LC", "HC")
legend(list(x=0,y=1.3), legend=leg.txt, col= c("red","black"), pch = c(15,17))
text(CV.S, CV.A, metabnames, cex=0.5, pos=4, col="black")
abline(0,1)
abline(reg)

# Sample and Metabolite clustering -----------------------

data1 <- cbind(exp1HC, exp1LC1, exp1LC2, exp2HC, exp2LC1, exp2LC2, exp3HC, exp3LC1, exp3LC2, exp4HC, exp4LC1, exp4LC2, exp4LC3, exp6NC1, exp6NC2, exp6NC3, exp6NC4, exp6NC5, exp6VC1, exp6VC2, exp6VC3, exp6VC4, exp6VC5, exp7NC1, exp7NC2, exp7NC3, exp7NC4, exp7NC5, exp7HC1, exp7HC2, exp7HC3, exp7HC4, exp7HC5)
names1 <- c("Syn exp1 HC", "Syn exp1 3h LC", "Syn exp1 24h LC", "Syn exp2 HC", "Syn exp2 3h LC", "Syn exp2 24h LC", "Syn exp3 HC", "Syn exp3 3h LC", "Syn exp3 24h LC", "Ara exp4 light HC", "Ara exp4 light 1d LC", "Ara exp4 light 3d LC", "Ara exp4 light 5d LC", "Ara exp6 2h NC", "Ara exp6 4h NC", "Ara exp6 6h NC", "Ara exp6 8h NC", "Ara exp6 10h NC", "Ara exp6 2h VC", "Ara exp6 4h VC", "Ara exp6 6h VC", "Ara exp6 8h VC", "Ara exp6 10h VC", "Ara exp7 2h NC", "Ara exp7 4h NC", "Ara exp7 6h NC", "Ara exp7 8h NC", "Ara exp7 10h NC", "Ara exp7 2h HC", "Ara exp7 4h HC", "Ara exp7 6h HC", "Ara exp7 8h HC", "Ara exp7 10h HC")
colnames(data1) <- names1
data3 <- scale(data1)

#First visualization, to check outcome

# Color function to generate green-red heat maps
color <- function(n = 50, low.col = 0.45, high.col=1, saturation = 1) { 
  if (n < 2) stop("n must be greater than 2")
  n1 <- n%/%2
  n2 <- n - n1
  c(hsv(low.col, saturation, seq(1,0,length=n1)), hsv(high.col, saturation, seq(0,1,length=n2))) 
}

par(oma = c(5,1,1,1))
samp5 <- hclust(as.dist(1-cor(data3, method="pearson")), method="average")
metab5 <- hclust(as.dist(1-cor(t(data3), method="pearson")), method="average")
heatmap(data3, Rowv=as.dendrogram(metab5), Colv=as.dendrogram(samp5), col=color(), scale="none", cexCol =0.8)

# Number of clusters ------------------------------------------------------
### See also File 15-10-27_Number-of_clusters.R

num.clus <- function (x){
  for (i in 2:10){
    clu <- pam(x, k=i)
    cat(i," avg: ", round(clu$silinfo$avg.width,2))
    cat(" clus: ", round(clu$silinfo$clus.avg.widths,2),"\n")
  }
}

num.clus(data3)
#Metabolites, 2 clusters (since not really used in publication, this is fine)

num.clus(t(data3))
#Samples, 6 clusters chosen

data3.clu.samp <- pam(t(data3), k=6)
data3.clu.metab <- pam(data3, k=2)

#Create Figure

pamsamp <- sample(rainbow(6)) 
pamsamp <- pamsamp[as.vector(data3.clu.samp$clustering)]
pammeta <- sample(rainbow(2))
pammeta <- pammeta[as.vector(data3.clu.metab$clustering)]

par(oma=c(7,2,2,5))
heatmap.2(data3, Rowv=as.dendrogram(metab5), Colv=as.dendrogram(samp5), col=color(), scale="none", RowSideColors=pammeta, ColSideColors=pamsamp, trace="none", density.info="none")

#Read-out clusters

Samp.PAM <- cbind(t(data3), data3.clu.samp$clustering)
Metab.PAM <- cbind(data3, data3.clu.metab$clustering)

Samp.clu1 <- Samp.PAM[Samp.PAM[,ncol(Samp.PAM)]==1,]
rownames(Samp.clu1)
Samp.clu2 <- Samp.PAM[Samp.PAM[,ncol(Samp.PAM)]==2,]
rownames(Samp.clu2)
Samp.clu3 <- Samp.PAM[Samp.PAM[,ncol(Samp.PAM)]==3,]
rownames(Samp.clu3)
Samp.clu4 <- Samp.PAM[Samp.PAM[,ncol(Samp.PAM)]==4,]
rownames(Samp.clu4)
Samp.clu5 <- Samp.PAM[Samp.PAM[,ncol(Samp.PAM)]==5,]
rownames(Samp.clu5)
Samp.clu6 <- Samp.PAM[Samp.PAM[,ncol(Samp.PAM)]==6,]
rownames(Samp.clu6)

# Calculate ratio matrices function --------------------------------------------

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

exp1HCMat <- ratioMat(exp1HC)
exp1LC1Mat <- ratioMat(exp1LC1)
exp1LC2Mat <- ratioMat(exp1LC2)
exp2HCMat <- ratioMat(exp2HC)
exp2LC1Mat <- ratioMat(exp2LC1)
exp2LC2Mat <- ratioMat(exp2LC2)
exp3HCMat <- ratioMat(exp3HC)
exp3LC1Mat <- ratioMat(exp3LC1)
exp3LC2Mat <- ratioMat(exp3LC2)

# Calculate ratio matrices Ara --------------------------------------------

exp4HCMat <- ratioMat(exp4HC)
exp4LC1Mat <- ratioMat(exp4LC1)
exp4LC2Mat <- ratioMat(exp4LC2)
exp4LC3Mat <- ratioMat(exp4LC3)
exp5HCMat <- ratioMat(exp5HC)
exp5LC1Mat <- ratioMat(exp5LC1)
exp5LC2Mat <- ratioMat(exp5LC2)
exp5LC3Mat <- ratioMat(exp5LC3)
exp6NC1Mat <- ratioMat(exp6NC1)
exp6NC2Mat <- ratioMat(exp6NC2)
exp6NC3Mat <- ratioMat(exp6NC3)
exp6NC4Mat <- ratioMat(exp6NC4)
exp6NC5Mat <- ratioMat(exp6NC5)
exp6NC6Mat <- ratioMat(exp6NC6)
exp6NC7Mat <- ratioMat(exp6NC7)
exp6VC1Mat <- ratioMat(exp6VC1)
exp6VC2Mat <- ratioMat(exp6VC2)
exp6VC3Mat <- ratioMat(exp6VC3)
exp6VC4Mat <- ratioMat(exp6VC4)
exp6VC5Mat <- ratioMat(exp6VC5)
exp6VC6Mat <- ratioMat(exp6VC6)
exp6VC7Mat <- ratioMat(exp6VC7)
exp7NC1Mat <- ratioMat(exp7NC1)
exp7NC2Mat <- ratioMat(exp7NC2)
exp7NC3Mat <- ratioMat(exp7NC3)
exp7NC4Mat <- ratioMat(exp7NC4)
exp7NC5Mat <- ratioMat(exp7NC5)
exp7NC6Mat <- ratioMat(exp7NC6)
exp7NC7Mat <- ratioMat(exp7NC7)
exp7HC1Mat <- ratioMat(exp7HC1)
exp7HC2Mat <- ratioMat(exp7HC2)
exp7HC3Mat <- ratioMat(exp7HC3)
exp7HC4Mat <- ratioMat(exp7HC4)
exp7HC5Mat <- ratioMat(exp7HC5)
exp7HC6Mat <- ratioMat(exp7HC6)
exp7HC7Mat <- ratioMat(exp7HC7)

# Create list of ratio matrices --------------------------------------------

RMlist1 <- list(exp1HCMat, exp1LC1Mat, exp1LC2Mat, exp2HCMat, exp2LC1Mat, exp2LC2Mat, exp3HCMat, exp3LC1Mat, exp3LC2Mat, exp4HCMat, exp4LC1Mat, exp4LC2Mat, exp4LC3Mat, exp6NC1Mat, exp6NC2Mat, exp6NC3Mat, exp6NC4Mat, exp6NC5Mat, exp6VC1Mat, exp6VC2Mat, exp6VC3Mat, exp6VC4Mat, exp6VC5Mat, exp7NC1Mat, exp7NC2Mat, exp7NC3Mat, exp7NC4Mat, exp7NC5Mat, exp7HC1Mat, exp7HC2Mat, exp7HC3Mat, exp7HC4Mat, exp7HC5Mat)

# Calculate RV coefficients and create new matrix --------------------------------------------

names1 <- c("Syn exp1 HC", "Syn exp1 3h LC", "Syn exp1 24h LC", "Syn exp2 HC", "Syn exp2 3h LC", "Syn exp2 24h LC", "Syn exp3 HC", "Syn exp3 3h LC", "Syn exp3 24h LC", "Ara exp4 light HC", "Ara exp4 light 1d LC", "Ara exp4 light 3d LC", "Ara exp4 light 5d LC", "Ara exp6 2h NC", "Ara exp6 4h NC", "Ara exp6 6h NC", "Ara exp6 8h NC", "Ara exp6 10h NC", "Ara exp6 2h VC", "Ara exp6 4h VC", "Ara exp6 6h VC", "Ara exp6 8h VC", "Ara exp6 10h VC", "Ara exp7 2h NC", "Ara exp7 4h NC", "Ara exp7 6h NC", "Ara exp7 8h NC", "Ara exp7 10h NC", "Ara exp7 2h HC", "Ara exp7 4h HC", "Ara exp7 6h HC", "Ara exp7 8h HC", "Ara exp7 10h HC")

rv_matrix1 <- matrix(nrow = length(RMlist1), ncol = length(RMlist1))

colnames(rv_matrix1) <- names1
rownames(rv_matrix1) <- names1

for(i in 1:length(RMlist1)) {
  tmp1 <- vector()
  for(j in 1:length(RMlist1)) {
    tmp <- coeffRV(RMlist1[[i]], RMlist1[[j]])$rv
    tmp1 <- c(tmp1, tmp)
  }
  rv_matrix1[i,] <- tmp1
}

fix(rv_matrix)

write.csv(rv_matrix1,'rv_matrix1.csv')

# Cluster matrix --------------------------------------------

num.clus(rv_matrix1)
#Number of clusters=10

#Run pam
RVclu1 <- pam(rv_matrix1, k=6)

RVcolors <- sample(rainbow(6)) 
RVcolors <- RVcolors[as.vector(RVclu1$clustering)]

#Plot clusters (simple)

pear <- hclust(as.dist(1-cor(rv_matrix1, method="pearson")), method="average")
par(oma = c(1,1,1,1))
plot(pear, main="RV coefficients between metabolite ratio matrices", hang=-1, cex=1)

par(oma = c(5,1,1,5))
heatmap(rv_matrix1, col=color(), scale="none", Rowv=as.dendrogram(pear), Colv="Rowv", ColSideColors=RVcolors, RowSideColors=RVcolors, main="PAM clustering of RV coefficients between metabolite ratio matrices (k=7)")

#Plot clusters and sort by cluster number

pcl <- as.vector(RVclu1$clustering)

pam1 <- cbind(rv_matrix1, pcl)

pam2 <- pam1[order(pam1[,34]),]

pam3 <- t(pam2[,-34])

pam4 <- cbind(pam3, pcl)

pam5 <- pam4[order(pam4[,34]),]

pam6 <- pam5[,-34]

RVcolors <- sample(rainbow(6)) 
RVcolors <- RVcolors[as.vector(pam2[,34])]

heatmap.2(pam6, col=color(), scale="none", Rowv=FALSE, Colv=FALSE, dendrogram="none", ColSideColors=RVcolors, RowSideColors=RVcolors, trace="none", density.info="none")

#Read-out clusters

PAM1 <- cbind(rv_matrix1, RVclu1$clustering)

cl1 <- PAM1[PAM1[,ncol(PAM1)]==1,]
rownames(cl1)
cl2 <- PAM1[PAM1[,ncol(PAM1)]==2,]
rownames(cl2)
cl3 <- PAM1[PAM1[,ncol(PAM1)]==3,]
rownames(cl3)
cl4 <- PAM1[PAM1[,ncol(PAM1)]==4,]
rownames(cl4)
cl5 <- PAM1[PAM1[,ncol(PAM1)]==5,]
rownames(cl5)
cl6 <- PAM1[PAM1[,ncol(PAM1)]==6,]
rownames(cl6)
