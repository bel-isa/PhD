#http://stackoverflow.com/questions/15376075/cluster-analysis-in-r-determine-the-optimal-number-of-clusters

mydata <- Reduced1
d <- Reduced1

# One: Look for a bend or elbow in the sum of squared error (SSE)  --------

wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
for (i in 2:30) wss[i] <- sum(kmeans(mydata, centers=i)$withinss)

par(oma=c(5,2,5,2))
plot(1:30, wss, type="o", ylab="Within groups sum of squares", xlab="Number of Clusters", main="Sum of squares error")


# Two: partitioning around medoids using the pamk function in the  --------

library(fpc)

pamk.best <- pamk(d)
cat("number of clusters estimated by optimum average silhouette width:", pamk.best$nc, "\n")
plot(pam(d, pamk.best$nc))

# Three: Calinsky criterion -----------------------------------------------

require(vegan)

fit <- cascadeKM(d, 1, 30, iter = 1000)

plot(fit, sortg = TRUE, grpmts.plot = TRUE)
calinski.best <- as.numeric(which.max(fit$results[2,]))
cat("Calinski criterion optimal number of clusters:", calinski.best, "\n")

# Four: Bayesian Information Criterion for expectation-maximizatio --------

# See http://www.jstatsoft.org/v18/i06/paper http://www.stat.washington.edu/fraley/mclust/tr504.pdf
library(mclust)

# Run the function to see how many clusters it finds to be optimal, set it to search for at least 1 model and up 30
d_clust <- Mclust(as.matrix(d), G=1:30)
m.best <- dim(d_clust$z)[2]
cat("model-based optimal number of clusters:", m.best, "\n")

plot(d_clust)

# Five: Affinity propagation (AP) clustering ------------------------------

library(apcluster)

d.apclus <- apcluster(negDistMat(r=2), d)
cat("affinity propogation optimal number of clusters:", length(d.apclus@clusters), "\n")

heatmap(d.apclus, main="Affinity propagation clustering")
plot(d.apclus, d)

# Six: Gap Statistic ------------------------------------------------------

library(cluster)

luecke <- clusGap(d, pam, 30, B = 100, verbose = interactive())

plot(1:30, luecke$Tab[,3], type="o", xlab="Number of clusters", ylab="Gap", main="Gap Statistics")
