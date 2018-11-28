#### Data processing ####
source('final-project/R/store_plant_data.R')
source('final-project/R/store_rodent_data.R')
source('final-project/R/store_adjusted_data.R')

library(dplyr)

store_plant_data()
store_rodent_data()

rodents <- read.csv('final-project/data/summer-rodents-raw.csv', 
                    stringsAsFactors = F)
summer_plants <- read.csv('final-project/data/summer-plants-raw.csv', 
                          stringsAsFactors = F)
winter_plants <- read.csv('final-project/data/winter-plants-raw.csv', 
                          stringsAsFactors = F)
all_plants <- rbind(summer_plants, winter_plants)

store_adjusted_data(plant_data = all_plants, rodent_data = rodents,
                    focal_type = 'rodent')
store_adjusted_data(plant_data = summer_plants, focal_type = 'plant',
                    season = 'summer')
store_adjusted_data(plant_data = winter_plants, focal_type = 'plant',
                    season = 'winter')

rm(list=ls())

#### Data exploration ####
summer_plants <- read.csv('final-project/data/summer-plants-adjusted.csv', 
                          stringsAsFactors = F)

 # logColSums <- log(colSums(summer_plants))
 # hist(logColSums)
# 
# veryrare <- which(logColSums <= 2)
# rare <- which(logColSums <= 4)
# common <- which(logColSums >= 8)
# verycommon <- which(logColSums >= 10)

summer_plants <- select(summer_plants, -year)
species <- colnames(summer_plants)
species <- as.data.frame(species)
species$occurrences <- 0

for(i in 1:nrow(species)) {
  not0 <- as.integer(length(which(summer_plants[,i] >= 1)))
  species$occurrences[i] <- not0
}

hist(species$occurrences)


length(which(species$occurrences >= 5))
filter(species, occurrences >= 30)

species$maximum <- 0

for(i in 1:nrow(species)){
  species$maximum[i] <- max(summer_plants[,i], na.rm = T)
}

hist(species$maximum[ which(species$maximum < 50)])

filter(species, maximum < 25)

# OK so we're going to go with removing species whose maximum abundance
# in any year is less than 50. 

rm(list=ls())

#### PCoA ####
summer_plants <- read.csv('final-project/data/summer-plants-adjusted.csv', 
                          stringsAsFactors = F)

summer_plants <- summer_plants %>%
  select(-year)

species <- colnames(summer_plants)
species <- as.data.frame(species)
species$maximum <- 0
for(i in 1:nrow(species)){
  species$maximum[i] <- max(summer_plants[,i], na.rm = T)
}

summer_plants <- summer_plants %>%
  select(which(species$maximum >= 50))

birds <- read.csv('week-4/Current_Hawaiian_Birds.csv', row = 1,
                  header = T)
birds2 <- read.csv('week-4/combined_birds.csv', row = 1,
                   header = T)
tree <- read.csv('week-4/tree.csv', row = 1, header = T)

library(vegan)
library(ca)

#### PCoordinatesA ####

# Bray-curtis index on current birds data
# this is the same as sorenson's, if used with 1 and 0
jbirds <- vegdist(birds, 'bray')

cmd <- cmdscale(jbirds, k = 5, eig = T)

str(cmd)
cmd$points

# PCoordinatesA table to look at the eigenvalues
# and the proportion of variance they capture

eigenvalues <- cmd$eig[1:5]
propVar <- eigenvalues/sum(eigenvalues)
cumVar <- cumsum(propVar)
PCoA_Table <- cbind(eigenvalues, propVar, cumVar)
PCoA_Table


# Scree plot:
plot(eigenvalues)
lines(lowess(eigenvalues))

# if it's not totally intractable, I'd keep 3 axes
# to explain a cumulative 91% of variance

# Plotting the first two PCoA axes

x <- cmd$points[,1]
y <- cmd$points[,2]
plot(x, y, xlab = 'Coord 1', ylab = 'Coord 2',
     xlim = range(x) * 1.2, 
     ylim = range(y) * 1.2, 
     type = 'n')
text(x,y, labels = rownames(cmd$points), cex = 0.9)


# Another (kind of terrifyingly messy) way to plot


ordiplot(scores(cmd)[, c(1, 2)], type = "t", cex = 1, main = "Hawaiian Bird PCoA")
## species scores not available
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)

# Add species
`?`(wascores)

species <- wascores(cmd$points[, 1:2], birds)
text(species, rownames(species), cex = 0.7, col = "red")

