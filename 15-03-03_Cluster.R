#Full Z transformed dataset stored in variable Z

row.names(Z) <- Metabs[,1]

T <- t(Z)

attach(T)

# Calculating pairwise Euclidean distances between the (standardiz --------

dist.T <- dist(T)

#Single linkage

T.single.link <- hclust(dist.T, method='single')

#Plotting the single linkage dendrogram

plclust(T.single.link, labels=rownames(T), ylab="Distance")

#complete linkage

T.complete.link <- hclust(dist.T, method='complete')

#Plotting the complete linkage dendrogram

plclust(T.complete.link, labels=rownames(T), ylab="Distance")

#Average linkage

T.avg.link <- hclust(dist.T, method='average')

#Plotting the average linkage dendrogram

plclust(T.avg.link, labels=rownames(T), ylab="Distance")


# Cutting the complete-linkage dendrogram to form k=19 clusters he --------

cut.19 <- cutree(T.complete.link, k=19)
cut.19     
#printing the "clustering vector"

T.19.clust <- lapply(1:19, function(nc) rownames(T)[cut.19==nc])  
T.19.clust   
#printing the clusters in terms of the labels


# Visualization of Clusters ----------------------------------------------

#scatterplot matrix

pairs(T[,-1], panel=function(x,y) text(x,y,cut.19))

#plot of the scores on the first 2 principal components, with the clusters separated by color

T.pc <- princomp(T[,-1],cor=TRUE)

# Setting up the colors for the 5 clusters on the plot
my.color.vector <- rep(c(139)
my.color.vector[cut.19==2] <- c(26)
my.color.vector[cut.19==3] <- c(135)
my.color.vector[cut.19==4] <- c(142)
my.color.vector[cut.19==5] <- c(451)
my.color.vector[cut.19==6] <- c(572)
my.color.vector[cut.19==7] <- c(635)
my.color.vector[cut.19==8] <- c(603)
my.color.vector[cut.19==9] <- c(77)
my.color.vector[cut.19==10] <- c(12)
my.color.vector[cut.19==11] <- c(137)
my.color.vector[cut.19==12] <- c(417)
my.color.vector[cut.19==13] <- c(261)
my.color.vector[cut.19==14] <- c(99)
my.color.vector[cut.19==15] <- c(259)
my.color.vector[cut.19==16] <- c(557)
my.color.vector[cut.19==17] <- c(491)
my.color.vector[cut.19==18] <- c(411)
my.color.vector[cut.19==19] <- c(81)

# Plotting the PC scores

par(pty="s")
plot(T.pc$scores[,1], T.pc$scores[,2], ylim=range(T.pc$scores[,1]), 
     xlab="PC 1", ylab="PC 2", type ='n', lwd=2)
text(T.pc$scores[,1], T.pc$scores[,2], labels=rownames(T), cex=0.7, lwd=2,
     col=my.color.vector)



# #### Everything for metabolites again -----------------------------------

# Calculating pairwise Euclidean distances between the (standardiz --------

dist.Z <- dist(Z)

#Single linkage

Z.single.link <- hclust(dist.Z, method='single')

#Plotting the single linkage dendrogram

plclust(Z.single.link, labels=rownames(Z), ylab="Distance")

#complete linkage

Z.complete.link <- hclust(dist.Z, method='complete')

#Plotting the complete linkage dendrogram

plclust(Z.complete.link, labels=rownames(Z), ylab="Distance")

#Average linkage

Z.avg.link <- hclust(dist.Z, method='average')

#Plotting the average linkage dendrogram

plclust(Z.avg.link, labels=rownames(Z), ylab="Distance")


# Cutting the complete-linkage dendrogram to form k=3 clusters he --------

cut.3 <- cutree(Z.complete.link, k=3)
cut.3     
#printing the "clustering vector"

Z.3.clust <- lapply(1:3, function(nc) rownames(Z)[cut.3==nc])  
Z.3.clust   
#printing the clusters in terms of the labels


# Visualization of Clusters ----------------------------------------------

#scatterplot matrix

par(mar = rep(2, 4))
pairs(Z[,-1], panel=function(x,y) text(x,y,cut.3))

#plot of the scores on the first 2 principal components, with the clusters separated by color

Z.pc <- prcomp(Z[,-1],cor=TRUE)

# Setting up the colors for the 3 clusters on the plot
my.color.vector1 <- rep(c(139), times=nrow(Z))
my.color.vector1[cut.3==2] <- c(26)
my.color.vector1[cut.3==3] <- c(135)


#Plotting the PC scores

par(pty="s")
plot(Z.pc$scores[,1], Z.pc$scores[,2], ylim=range(Z.pc$scores[,1]), 
     xlab="PC 1", ylab="PC 2", type ='n', lwd=2)
text(Z.pc$scores[,1], Z.pc$scores[,2], labels=rownames(Z), cex=0.7, lwd=2,
     col=my.color.vector1)