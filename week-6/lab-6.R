library(raster)
library(cluster)
library(vegan)
library(pvclust)

birds <- read.csv('week-6/Caribbean_birds.csv', row = 1, header = T, 
                  stringsAsFactors = F)
# Make a dissimilarity matrix
distBirds <- vegdist(birds, 'jaccard')

# Run various clustering algorithms
singleTree <- hclust(distBirds, method = "single")
completeTree <- hclust(distBirds, method = "complete")
centroidTree <- hclust(distBirds, method = "centroid")
medianTree <- hclust(distBirds, method = "median")
averageTree <- hclust(distBirds, method = "average")
wardTree <- hclust(distBirds, method = "ward.D2")

plot(singleTree)
plot(completeTree)
plot(centroidTree)
plot(medianTree)
plot(averageTree)
plot(wardTree)

par(mfrow = c(2, 3))
plot(singleTree)
plot(completeTree)
plot(centroidTree)
plot(medianTree)
plot(averageTree)
plot(wardTree)

# Jamaica, Cuba, Hispaniola, Puerto Rico consistently cluster together as
# a relatively-out group
# Generally, the clusters seem similar but their relative distances seem
# to change? it's very hard to look at these and extract comparative meaning. 


# Evaluating the clusters:

# agglomerative coefficient:
# how agglomerative is the data?
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

# cophenetic correlation coefficient
# how well do the dendrograms approximate the original distance matrix?
cc1 <- cor(distBirds, cophenetic(singleTree))
cc2 <- cor(distBirds, cophenetic(completeTree))
cc3 <- cor(distBirds, cophenetic(centroidTree))
cc4 <- cor(distBirds, cophenetic(medianTree))
cc5 <- cor(distBirds, cophenetic(averageTree))
cc6 <- cor(distBirds, cophenetic(wardTree))
cophCor <- round(c(cc1, cc2, cc3, cc4, cc5, cc6), 2)
methods <- c("single", "complete", "centroid", "median", "average", "ward")
dendrogramTable <- data.frame(methods, cophCor, agc)
dendrogramTable

# function telling the loop to read the input as text
e = function(expr) eval(parse(text = expr))

# sets up a variable to fill with the output of the loop
cc <- NULL

# list of names for the loop
methodList <- c("singleTree", "completeTree", "centroidTree", "medianTree", 
                "averageTree", "wardTree")

# run the loop
for (i in methodList) {
  cc[i] <- round(cor(distBirds, cophenetic(e(i))), 2)
}

cc

# Bootstrapping
# how many/which clusters are recovered?
# usually you should use more like 1000 bootstraps

boot1 <- pvclust(t(birds), method.hclust = "single", method.dist = "binary", 
                 nboot = 100)
boot2 <- pvclust(t(birds), method.hclust = "complete", method.dist = "binary", 
                 nboot = 100)
boot3 <- pvclust(t(birds), method.hclust = "average", method.dist = "binary", 
                 nboot = 100)

# Here, “binary” is Jaccard distance
par(mfrow = c(2, 3))

plot(boot1)
pvrect(boot1, alpha = 0.95, pv = "au")
plot(boot2)
pvrect(boot2, alpha = 0.95, pv = "au")
plot(boot3)
pvrect(boot3, alpha = 0.95, pv = "au")

plot(boot1)
pvrect(boot1, alpha = 0.95, pv = "bp")
plot(boot2)
pvrect(boot2, alpha = 0.95, pv = "bp")
plot(boot3)
pvrect(boot3, alpha = 0.95, pv = "bp")


# Polythetic divisive hierarchical clustering
# use diana()

diTree <- diana(distBirds)
plot(diTree, which.plots = 2)

diTree$dc
## [1] 0.6052537
# and calculate the cophenetic correlation coefficient:

d.coph <- cor(distBirds, cophenetic(diTree))


# now try with portal data?

rm(list=ls())

library(dplyr)

portal <- read.csv("week-6/portal_dat.csv")

tobin <- function(x) {
  if(x > 0) return(1)
  else return(x)
}

portal <- portal %>%
  mutate(dates = row_number()) %>%
  mutate(dates = as.character(dates))

portal[,2:22] <- apply(portal[,2:22], c(1,2), tobin)


