#7 Clusters

data3.clu.samp <- pam(t(data3), k=7)
data3.clu.metab <- pam(data3, k=2)

pamsamp <- sample(rainbow(7)) 
pamsamp <- pamsamp[as.vector(data3.clu.samp$clustering)]
pammeta <- sample(rainbow(2))
pammeta <- pammeta[as.vector(data3.clu.metab$clustering)]

par(oma=c(7,2,2,5))
heatmap.2(data3, Rowv=as.dendrogram(metab5), Colv=as.dendrogram(samp5), col=color(), scale="none", RowSideColors=pammeta, ColSideColors=pamsamp, trace="none", density.info="none")

#6 Clusters

data3.clu.samp <- pam(t(data3), k=6)
data3.clu.metab <- pam(data3, k=2)

pamsamp <- sample(rainbow(6)) 
pamsamp <- pamsamp[as.vector(data3.clu.samp$clustering)]
pammeta <- sample(rainbow(2))
pammeta <- pammeta[as.vector(data3.clu.metab$clustering)]

par(oma=c(7,2,2,5))
heatmap.2(data3, Rowv=as.dendrogram(metab5), Colv=as.dendrogram(samp5), col=color(), scale="none", RowSideColors=pammeta, ColSideColors=pamsamp, trace="none", density.info="none")

#5 Clusters

data3.clu.samp <- pam(t(data3), k=5)
data3.clu.metab <- pam(data3, k=2)

pamsamp <- sample(rainbow(5)) 
pamsamp <- pamsamp[as.vector(data3.clu.samp$clustering)]
pammeta <- sample(rainbow(2))
pammeta <- pammeta[as.vector(data3.clu.metab$clustering)]

par(oma=c(7,2,2,5))
heatmap.2(data3, Rowv=as.dendrogram(metab5), Colv=as.dendrogram(samp5), col=color(), scale="none", RowSideColors=pammeta, ColSideColors=pamsamp, trace="none", density.info="none")

#4 Clusters

data3.clu.samp <- pam(t(data3), k=4)
data3.clu.metab <- pam(data3, k=2)

pamsamp <- sample(rainbow(4)) 
pamsamp <- pamsamp[as.vector(data3.clu.samp$clustering)]
pammeta <- sample(rainbow(2))
pammeta <- pammeta[as.vector(data3.clu.metab$clustering)]

par(oma=c(7,2,2,5))
heatmap.2(data3, Rowv=as.dendrogram(metab5), Colv=as.dendrogram(samp5), col=color(), scale="none", RowSideColors=pammeta, ColSideColors=pamsamp, trace="none", density.info="none")
