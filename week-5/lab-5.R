library(raster)
library(cluster)


snails <- read.csv("week-5/snail_data.csv", row = 1, header = TRUE)[, 1:3]

snail.tot <- apply(snails, 2, sum)
cv(snail.tot)

#Set a vector for the loop to fill. 
wss <- rep(0, 8)

#Run a loop for 1 to 8 clusters:
for (i in 1:8) {
  
  wss[i] <- sum(kmeans(snails, centers = i,nstart=25)$withinss) # run the kmeans function for each number of clusters (i) and extract the within sum of squares for each.

}
#Check out you vector of within group sum of squares for one to eight groups:
wss 
plot(1:8, wss, type = "b", xlab = "Number of groups", ylab = "Within groups sum of squares")

sil <- rep(0, 8)
for (i in 2:8) sil[i] <- summary(silhouette(kmeans(snails, centers = i, iter.max = 100, 
                                                   nstart = 25)$cluster, dist(snails)))$avg.width
plot(2:8, sil[2:8], type = "b", xlab = "Number of groups", ylab = "average silhouette width ")

# 2 clusters seems to be good - you don't get a better
# silhoutte width with additional clusters.

snails.kop <- kmeans(snails, centers = 2, iter.max = 10, nstart = 25)

pairs(snails, panel = function(x, y, z) text(x, y, snails.kop$cluster))

# The groups differentiate along spireheight.

snail.pc <- princomp(snails, cor = F)
summary(snail.pc)
snail.pc$loadings

my.color.vector <- rep("green", times = nrow(snails))
my.color.vector[snails.kop$cluster == 1] <- "blue"
my.color.vector[snails.kop$cluster == 2] <- "green"


plot(snail.pc$scores[, 1], snail.pc$scores[, 2], ylim = range(snail.pc$scores[, 
                                                                              1]), xlim = range(snail.pc$scores[, 1] * 1.25), xlab = "PC 1", ylab = "PC 2", 
     type = "n", lwd = 2)
text(snail.pc$scores[, 1], snail.pc$scores[, 2], labels = rownames(snails), 
     cex = 1.25, lwd = 2, col = my.color.vector)


biplot(snail.pc, xlabs = rep("", 9), xlim = range(-0.55, 0.55))
text(snail.pc$scores[, 1], snail.pc$scores[, 2], labels = rownames(snails), 
     cex = 1.25, lwd = 2, col = my.color.vector)


