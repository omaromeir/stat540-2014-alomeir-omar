library(RColorBrewer)
library(cluster)
library(pvclust)
library(xtable)
library(limma)
library(plyr)
library(lattice)

prDat <- read.table("../data/GSE4051_data.tsv", header = TRUE, row.names = 1) # the whole enchilada
str(prDat, max.level = 0)
prDes <- readRDS("../data/GSE4051_design.rds")
str(prDes)

sprDat <- t(scale(t(prDat)))
str(sprDat, max.level = 0, give.attr = FALSE)

round(data.frame(avgBefore = rowMeans(head(prDat)), avgAfter = rowMeans(head(sprDat)), varBefore = apply(head(prDat), 1, var), varAfter = apply(head(sprDat), 1, var)), 2)

########################## Sample clustering

########################## Heirarchical

pr.dis <- dist(t(sprDat), method = 'euclidean')

prDes$grp <- with(prDes, interaction(gType, devStage))
summary(prDes$grp)

pr.hc.s <- hclust(pr.dis, method = 'single')
pr.hc.c <- hclust(pr.dis, method = 'complete')
pr.hc.a <- hclust(pr.dis, method = 'average')
pr.hc.w <- hclust(pr.dis, method = 'ward')

op <- par(mar = c(0,4,4,2), mfrow = c(2,2))

plot(pr.hc.s, labels = FALSE, main = "Single", xlab = "")
plot(pr.hc.c, labels = FALSE, main = "Complete", xlab = "")
plot(pr.hc.a, labels = FALSE, main = "Average", xlab = "")
plot(pr.hc.w, labels = FALSE, main = "Ward", xlab = "")

par(op)

op <- par(mar = c(1,4,4,1))
plot(pr.hc.w, labels = prDes$grp, cex = 0.6, main = "Ward showing 10 clusters")
rect.hclust(pr.hc.w, k = 10)

par(op)

jGraysFun <- colorRampPalette(brewer.pal(n = 9, "Greys"))
gTypeCols <- brewer.pal(11, "RdGy")[c(4,7)]
heatmap(as.matrix(sprDat), Rowv = NA, col = jGraysFun(256), hclustfun = function(x) hclust(x, method = 'ward'), scale = "none", labCol = prDes$grp, labRow = NA, margins = c(8,1), ColSideColor = gTypeCols[unclass(prDes$gType)])
legend("topright", legend = levels(prDes$gType), col = gTypeCols, lty = 1, lwd = 5, cex = 0.5)

########################## Exercise

# Trying a different method
heatmap(as.matrix(sprDat), Rowv = NA, col = jGraysFun(256), hclustfun = function(x) hclust(x, method = 'complete'), scale = "none", labCol = prDes$grp, labRow = NA, margins = c(8,1), ColSideColor = gTypeCols[unclass(prDes$gType)])

# Clustering the original data
heatmap(as.matrix(prDat), Rowv = NA, col = jGraysFun(256), hclustfun = function(x) hclust(x, method = 'average'), scale = "none", labCol = prDes$grp, labRow = NA, margins = c(8,1), ColSideColor = gTypeCols[unclass(prDes$gType)])

########################## Partitioning

set.seed(31)
k <- 5
pr.km <- kmeans(t(sprDat), centers = k, nstart =  50)

pr.km$withinss

pr.kmTable <- data.frame(devStage = prDes$devStage, cluster = pr.km$cluster)
prTable  <-  xtable(with(pr.kmTable, table(devStage,cluster)), caption='Number of samples from each develomental stage within each k-means cluster')

align(prTable) <- "lccccc"
print(prTable, type = 'html', caption.placement = 'top')

########################## Repeat the analysis using a different seed and check if you get the same clusters.

set.seed(15)
k <- 5
pr.km <- kmeans(t(sprDat), centers = k, nstart =  50)

pr.km$withinss

pr.kmTable <- data.frame(devStage = prDes$devStage, cluster = pr.km$cluster)
prTable  <-  xtable(with(pr.kmTable, table(devStage,cluster)), caption='Number of samples from each develomental stage within each k-means cluster')

align(prTable) <- "lccccc"
print(prTable, type = 'html', caption.placement = 'top')

########################## Pam

pr.pam <- pam(pr.dis, k = k)
pr.pamTable <- data.frame(devStage = prDes$devStage, cluster = pr.pam$clustering)
pamTable  <-  xtable(with(pr.pamTable, table(devStage, cluster)), caption='Number of samples from each develomental stage within each PAM cluster')

align(pamTable) <- "lccccc"
print(pamTable, type = 'html', caption.placement = 'top')

op <- par(mar = c(5,1,4,4))
plot(pr.pam, main = "Silhouette Plot for 5 clusters")

par(op)

########################## Exercise
k = 1:5
plot(k, pr.pam$silinfo$clus.avg.widths)

########################## Gene clustering

desMat <- model.matrix(~devStage, prDes)
prFit <- lmFit(prDat, desMat)
ebFit <- eBayes(prFit)

hits <- topTable(ebFit, coef=grep("devStage", colnames(coef(ebFit))), p.value=1e-5, number=972)

topDat <- sprDat[rownames(hits), ]

geneC.dis <- dist(topDat, method = 'euclidean')

geneC.hc.a <- hclust(geneC.dis, method = 'average')

plot(geneC.hc.a, labels = FALSE, main = "Hierarchical with Average Linkage", xlab = "")

########################## Partitioning

set.seed(1234)
k <- 5
kmeans.genes <- kmeans(topDat, centers = k)

clusterNum <- 1 

plot(kmeans.genes$centers[clusterNum, ], ylim = c(-4, 4), type = 'n', xlab = "Samples", ylab = "Relative expression" ) 

matlines(y = t(topDat[kmeans.genes$cluster == clusterNum, ]), col = 'grey') 

points(kmeans.genes$centers[clusterNum, ], type = 'l') 

points(kmeans.genes$centers[clusterNum, ],  col = prDes$devStage, pch = 20) 

devStageCols <- brewer.pal(11, "RdGy")[c(2,4,7,9,11)]
heatmap(as.matrix(topDat), col = jGraysFun(256),
        hclustfun = function(x) hclust(x, method = 'average'),
        labCol = prDes$grp, labRow = NA, margin = c(8,1), scale = "none",
        ColSideColor = devStageCols[unclass(prDes$devStage)])
legend("topleft", levels(prDes$devStage), col = devStageCols,
       lty = 1, lwd = 5, cex = 0.5)

annoTopDat <- stack(as.data.frame(topDat)) # stack probe data tall and skinny
annoTopDat$probeset <- rownames(topDat) # add probeset ID as variable
## get info on gType and devStage, then average over reps within devStage
annoTopDat <- merge(annoTopDat, prDes, by.x = "ind", by.y = "sidChar")
devStageAvg <- ddply(annoTopDat, ~ probeset, function(x) {
  avgByDevStage <- aggregate(values ~ devStage, x, mean)$values
  names(avgByDevStage) <- levels(x$devStage)
  avgByDevStage
})
## put probset info back into rownames
rownames(devStageAvg) <- devStageAvg$probeset
devStageAvg$probeset <- NULL
str(devStageAvg)

heatmap(as.matrix(devStageAvg), Colv = NA, col = jGraysFun(256),
        hclustfun = function(x) hclust(x,method = 'average'),
        labCol = colnames(devStageAvg), labRow = NA, margin = c(8,1))

k <- 4
geneDS.km <- kmeans(devStageAvg, centers = k, nstart = 50)
clust.centers <- geneDS.km$centers

#Look at all clusters
op <- par(mfrow = c(2, 2))
for(clusterNum in 1:4) {
  # Set up the axes without plotting; ylim set based on trial run.
  plot(clust.centers[clusterNum,], ylim = c(-4,4), type='n',
       xlab = "Develomental Stage", ylab = "Relative expression",
       axes = F, main = paste("Cluster", clusterNum, sep = " ")) 
  axis(2)
  axis(1, 1:5, c(colnames(clust.centers)[1:4],"4W"), cex.axis = 0.9)
  
  # Plot the expression of all the genes in the selected cluster in grey.
  matlines(y = t(devStageAvg[geneDS.km$cluster == clusterNum, ]),
           col = 'grey') 
  
  # Add the cluster center. This is last so it isn't underneath the members
  points(clust.centers[clusterNum, ] , type = 'l') 
  
  # Optional: points to show development stages.
  points(clust.centers[clusterNum, ],  pch = 20)
} 

par(op)

plot(clust.centers[clusterNum, ], ylim = c(-4, 4), type = 'n',
     xlab = "Develomental Stage", ylab = "Average expression",
     axes = FALSE, main = "Clusters centers") 
axis(2)
axis(1, 1:5, c(colnames(clust.centers)[1:4],"4W"), cex.axis = 0.9)

for(clusterNum in 1:4) {
  points(clust.centers[clusterNum,], type = 'l', col = clusterNum, lwd=2) 
  points(clust.centers[clusterNum,] , col = clusterNum, pch = 20)
}

cloud(devStageAvg[ ,"E16"] ~ devStageAvg[ ,"P6"] * devStageAvg[ ,"4_weeks"], col = geneDS.km$clust, xlab = "E16", ylab = "P6", zlab = "4_weeks")

pvc <- pvclust(topDat, nboot = 100)

plot(pvc, labels = prDes$grp, cex = 0.6)
pvrect(pvc, alpha = 0.95) 

########################## PCA
pcs <- prcomp(sprDat, center = F, scale = F) 

# scree plot
plot(pcs) 

# append the rotations for the first 10 PCs to the phenodata
prinComp <- cbind(prDes, pcs$rotation[prDes$sidNum, 1:10]) 

# scatter plot showing us how the first few PCs relate to covariates
plot(prinComp[ ,c("sidNum", "devStage", "gType", "PC1", "PC2", "PC3")], pch = 19, cex = 0.8) 

plot(prinComp[ ,c("PC1","PC2")], bg = prDes$devStage, pch = 21, cex = 1.5)
legend(list(x = 0.2, y = 0.3), as.character(levels(prDes$devStage)), pch = 21, pt.bg = c(1,2,3,4,5))

######################## Exercise

#one plot in lattice
xyplot(PC1 ~ PC2, prinComp, groups=devStage, auto.key=TRUE)