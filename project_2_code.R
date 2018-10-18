#### K-means clustering of Darlingtonia ####
library(dplyr)
library(cluster)
library(vegan)
library(pvclust)

# Load darlingtonia
darling <- read.csv("project-1/Darlingtonia.csv", header = T,
                    stringsAsFactors = F, row.names =1 )

darling <- darling %>%
  mutate(tube_diam = log1p(tube_diam),
         keel_diam = log1p(keel_diam),
         wing2_length = log1p(wing2_length),
         hoodarea = log1p(hoodarea),
         wingarea = log1p(wingarea),
         tubearea = log1p(tubearea)) %>%
  scale()

max_clusters = nrow(darling) - 1

wss = rep(0, max_clusters)

for(i in 1:max_clusters){
 
  wss[i] = sum(kmeans(darling, centers = i,nstart=25)$withinss)
   
}

plot(1:10, wss[1:10], type = "b", xlab = "Number of groups", ylab = "Within groups sum of squares")

# 2 seems to be the inflection point for the scree plot.

sil <- rep(0, max_clusters)
for (i in 2:max_clusters) sil[i] <- summary(silhouette(kmeans(darling, centers = i, iter.max = 100, 
                                                   nstart = 25)$cluster, dist(darling)))$avg.width
plot(2:10, sil[2:10], type = "b", xlab = "Number of groups", ylab = "average silhouette width ")

# 2 has the highest silhouette width.

darling.kop <- kmeans(darling, centers = 2, iter.max = 100, nstart= 25)

pairs(summary(darling.pc), panel = function(x, y, z) text(x, y, darling.kop$cluster))

summary(darling.pc)

darling.pc <- princomp(darling, cor = F)
summary(darling.pc)
darling.pc$loadings

color.vector <- rep("green", times = nrow(darling))
color.vector[darling.kop$cluster == 1] <- "blue"
color.vector[darling.kop$cluster == 2] <- "green"


plot(darling.pc$scores[, 1], darling.pc$scores[, 2], ylim = range(darling.pc$scores[,1]), 
     xlim = range(darling.pc$scores[, 1] * 1.25), xlab = "PC 1", ylab = "PC 2",
     type = "n", lwd = 2)
points(darling.pc$scores[, 1], darling.pc$scores[, 2],
     cex = 1.25, lwd = 2, col = color.vector)


biplot(darling.pc, xlabs = rep("", nrow(darling)), xlim = range(-0.55, 0.55))
points(darling.pc$scores[, 1], darling.pc$scores[, 2], 
     cex = 1.25, lwd = 2, col = color.vector)
title(main = 'Darlingtonia biplot')


darling.pc$loadings

# Comp 1 best separates the clusters
# It is kind of a composite of hood area, wing length, and hood area.

darling.kop$centers
# center for wing area is .599 for cluster 1 and -.84 for cluster 2.
mean(darling[which(darling.kop$cluster == 1),'wingarea'])
mean(darling[which(darling.kop$cluster == 2),'wingarea'])

# MRPP
groups <- darling.kop$cluster

distDarling <- vegdist(darling, method = 'euclidean')
darlingMRPP <- mrpp(distDarling, groups, permutations = 1000)

hist(darlingMRPP$boot.deltas, main = 'Histogram of darling MRPP deltas')
points(darlingMRPP$delta, 0, pch = 19, col = 'red', bg= 'red', cex = 2)
darlingMRPP$E.delta # 4.143135
darlingMRPP$delta # actual delta 3.485
darlingMRPP$Pvalue # 0.00099

# Conclude: significantly different groups

darlingAnosim <- anosim(distDarling, groups, permutations = 1000)
darlingAnosim 
# R = 0.5159
# Significance: 0.000999

hist(darlingAnosim$perm, main = "Histogram of R statistics for Darlingtonia", 
     xlim = c(-0.5, 1))
points(darlingAnosim$statistic, 0, pch = 19, col = "red", bg = "red", cex = 2)

# Clusters are significantly different.

rm(list=ls())

#### Polythetic agglomerative hierarchical clustering ####

dunes <- read.csv('project-2/dune_data.csv', stringsAsFactors = F, row.names = 1)

distDunes <- vegdist(dunes, 'jaccard')


# Run various clustering algorithms
singleTree <- hclust(distDunes, method = "single")
completeTree <- hclust(distDunes, method = "complete")
centroidTree <- hclust(distDunes, method = "centroid")
medianTree <- hclust(distDunes, method = "median")
averageTree <- hclust(distDunes, method = "average")
wardTree <- hclust(distDunes, method = "ward.D2")

plot(singleTree)
plot(completeTree)
plot(centroidTree)
plot(medianTree)
plot(averageTree)
plot(wardTree)

ag1 <- coef.hclust(singleTree)
ag2 <- coef.hclust(completeTree)
ag3 <- NA
ag4 <- NA
ag5 <- coef.hclust(averageTree)
ag6 <- coef.hclust(wardTree)
methods <- c("single", "complete", "centroid", "median", "average", "ward")
agc <- round(c(ag1, ag2, ag3, ag4, ag5, ag6), 2)
agcTable <- data.frame(methods, agc)
agcTable

cc1 <- cor(distDunes, cophenetic(singleTree))
cc2 <- cor(distDunes, cophenetic(completeTree))
cc3 <- cor(distDunes, cophenetic(centroidTree))
cc4 <- cor(distDunes, cophenetic(medianTree))
cc5 <- cor(distDunes, cophenetic(averageTree))
cc6 <- cor(distDunes, cophenetic(wardTree))
cophCor <- round(c(cc1, cc2, cc3, cc4, cc5, cc6), 2)
methods <- c("single", "complete", "centroid", "median", "average", "ward")
dendrogramTable <- data.frame(methods, cophCor, agc)
dendrogramTable

# The average linkage method has the highest cophenetic clustering
# 0.88

# The ward has the most cluster structure (.73)

# Bootstrapping

boot1 <- pvclust(t(dunes), method.hclust = "average", method.dist = "binary", 
                 nboot = 1000)
plot(boot1)
pvrect(boot1, alpha = 0.95, pv = "au")

# it looks like only 1 cluster emerges - but doesn't that mean 2 clusters?
# ASK BEN

#### Polythetic Divisive Hierarchical Clustering ####

diTree <- diana(distDunes)
plot(diTree, which.plots = 2)

# Divisive coeffiicent
diTree$dc
## [1] 0.5588

# Cophenetic correlation coefficient:

d.coph <- cor(distDunes, cophenetic(diTree))
d.coph
# 0.804

# The PAHC does slightly better (higher cophenetic coeff)
# The PDHC has less cluster structure. 

