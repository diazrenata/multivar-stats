install.packages('MVA')
install.packages('psych')
install.packages('Hmisc')
install.packages('vegan')
install.packages('StatMatch')
install.packages('MASS')
install.packages('raster')
install.packages('cluster')

library(MVA)
library(psych)
library(Hmisc)
library(vegan)
library(StatMatch)
library(MASS)
library(raster)
library(cluster)

# load air pollution from MVA, and modified air
usAir <- USairpollution

usAir_mod <- read.csv('week-1/data/usAir_mod.csv', row = 1, header = T)


# check for unrealistic values
describeBy(usAir)

# I don't see anything obviously odd, but I don't know all these variables

describeBy(usAir_mod)

# Max temp of 700 anythings is indeed odd.

## Check for missing data (NAs)
describe(usAir_mod)

# Missing values for wind (2) and precip (1)

# complete case - remove rows with missing data

usAir_mod[complete.cases(usAir_mod), ]

# alternatively, you could replace the missing values with means for the variable
# 'imputation'

meanz <- colMeans(usAir_mod, na.rm = T)
naFunc <- function(column) {
  column[is.na(column)] = round(mean(column, na.rm = T), 2)
  return(column)
}

Impute <- apply(usAir_mod, 2, naFunc)

describe(Impute)
# no more missing variables

## Transformations

# log transform

# start with SO2

hist(usAir$SO2)
usAirlog <- log1p(usAir)
hist(usAirlog$SO2)
par(mfrow = c(1,2))
hist(usAir[,1])
hist(usAirlog[,1])

for (i in 2:ncol(usAir)){
  par(mfrow = c(1,2))
  hist(usAir[,i])
  hist(usAirlog[,i])
}

# cols 2, 5, 7 look plausibly normal even when not log transformed

# sqrt
usAirsqrt <- sqrt(usAir)
hist(usAirsqrt$SO2)
par(mfrow = c(1, 2))

hist(usAir[, 1])
hist(usAirsqrt[, 1])


for (i in 2:ncol(usAir)){
  par(mfrow = c(1,2))
  hist(usAir[,i])
  hist(usAirsqrt[,i])
}


# arcsine - need to generate toy data

newData <- runif(100, 0, 1)
asin(sqrt(newData))
par(mfrow = c(1,2))
hist(newData)
hist(asin(sqrt(newData)))


# Data standardization

cSums <- colSums(usAir)
Sdev <- sd(cSums)
M <- mean(cSums)
Cv <- Sdev/M * 100

# welp, cv is 129, so... that's > 50

scaledData <- scale(usAir)
par(mfrow = c(1, 2))

hist(usAir[, 1], main = colnames(usAir)[1], xlab = " ")
hist(scaledData[, 1], main = colnames(usAir)[1], xlab = " ")

for (i in 1:ncol(scaledData)) {
  par(mfrow = c(1, 2))
  
  hist(usAir[, i], main = colnames(usAir)[i], xlab = " ")
  hist(scaledData[, i], main = colnames(usAir)[i], xlab = " ")
}



cSums <- round(colSums(scaledData), digits = 7)
Sdev <- sd(cSums)
M <- mean(cSums)
Cv <- Sdev/M * 100


# outliers

# looking at each variable seperately, z standardize and then look for 
# values outside 3 sds for each variable

scaledData <- scale(usAir)

par(mfrow = c(2, 4))
mapply(hist, as.data.frame(usAir), main = colnames(usAir), xlab = " ")

out <- function(x) {
  lier <- x[abs(x) > 3]
  return(lier)
}

apply(scaledData, 2, out)
# chicago is a weirdo

# skip multivar outliers for now 
# i.e. looking at all the variables together


# Distance and dissimilarity

# euclidian dist

scaledData <- scale(usAir)

eucDist <- vegdist(scaledData, 'euclidean')
hist(eucDist)
# so many of them are about equally distant from each other...but there are some weirdos


# some toy fruit data
Fruit <- rbind(c(1, 0, 1, 1), c(2, 1, 0, 0), c(3, 0, 4, 4))
colnames(Fruit) <- c("Farm", "Strawaberry", "Peach", "Rasberry")
Fruit

eucDist<- vegdist(Fruit[,-1], 'euclidean')

eucDist

# note that eucDist gets weird for 0s - if you have 0 overlap,
# you can still be closer than two that have v little overlap

# city block distance
# doesn't this depend strongly on the units?
cbDist <- vegdist(scaledData, "manhattan")

hist(cbDist)

# both are unimodal, but cbDist is generally higher?


# bray-curtis
brayDist <- vegdist(usAir, 'bray')
hist(brayDist)

# less pronuncedly unimodal

brayFruit <- vegdist(Fruit[, -1], 'bray')
brayFruit

# gower dissimilarity
# allows for mixed variables and missing values
# gower lives in the raster package, and it's called daisy.

daisy(usAir_mod, metric = 'gower')

# Multivariate outliers
brayDist <- vegdist(usAir, 'bray')
multOut <- scale(colMeans(as.matrix(brayDist)))
hist(multOut)
multOut[multOut > 3, ]

# alternatively find observations >3 sds from the mean
colBray <- colMeans(as.matrix(brayDist))
mBray <- mean(colBray)
stdBray <- sd(colBray)
threeSD <- stdBray * 3 + mBray
hist(colBray)
colBray[colBray>threeSD]
