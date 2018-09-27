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


# Running a PCoA by hand (step by step as in class)

jbirds <- vegdist(birds, "bray")
CORD <- -1/2 * jbirds^2
C <- as.matrix(CORD)
cs <- colMeans(C)
rs <- rowMeans(C)
C1 <- sweep(C, MARGIN = 2, cs, FUN = "-")
C2 <- sweep(C1, MARGIN = 1, rs, FUN = "-")
delta <- mean(C) + C2

# Next, run an eigen analysis:
EG <- eigen(delta)
eigenvalues2 <- EG$values[1:5]

# And make our PCoA table:
propVar2 <- eigenvalues2/sum(eigenvalues2)
cumVar2 <- cumsum(propVar2)
PCoA_Table2 <- cbind(eigenvalues2, propVar2, cumVar2)
PCoA_Table2

# You scale the eigenvectors by the square root of their eigenvalues to get
# the coordinates (points):
points2 <- sweep(EG$vectors[, 1:5], MARGIN = 2, sqrt(eigenvalues2), FUN = "*")
points2

# The coordinates:
plot(x, y, xlab = "Coordinate 1", ylab = "Coordinate 2", xlim = range(x) * 1.2, 
     ylim = range(y) * 1.2, type = "n")
text(x, y, labels = rownames(birds), cex = 0.9)


# Calculate weighted species scores:
scores1 <- sweep(birds, MARGIN = 1, x, FUN = "*")
species1 <- colSums(scores1)/colSums(birds)
scores2 <- sweep(birds, MARGIN = 1, y, FUN = "*")
species2 <- colSums(scores2)/colSums(birds)

# Add to the plot:

text(cbind(species1, species2), colnames(birds), cex = 0.7, col = "red")


#### NMDS ####

head(birds2)

# Calculate Sorenson's dissimilarity for the second
# birds dataset

jbirds2 <- vegdist(birds2, 'bray')


# Use metaMDS to get the NMDS 
?metaMDS

nmdsBird <- metaMDS(jbirds2, k = 2, trace = T, trymax = 100000000000)
stressplot(nmdsBird)

nmdsBird

# High r2 (.9, .98) for the stress plot tells us that 
# the axes accurately preserve the rank-dissimilarlities
# in the original data.
# Min stress of 0.126; this is 'fair'. 
# Stress is more informative than r2?


# Comparing the historical and current communities

treat = as.matrix(c(rep("Historical", 6), rep("Current", 6)))

ordiplot(nmdsBird, type = 'n', xlim = c(-.05, .5), ylim = c(-0.5, 0.5))
orditorp(nmdsBird, display = "sites", col = c(rep("green", 6), rep("blue", 6)), 
         air = 0.01, cex = 1.25)
legend(-0.55, 0.5, c("Historical", "Current"), cex = 0.8, col = c("green", "blue"), 
       pch = 15:15)

ordihull(nmdsBird, treat, display = "si", lty = 1, col = "green", show.groups = "Historical")
ordihull(nmdsBird, treat, display = "si", lty = 1, col = "blue", show.groups = "Current")


#### Correspondence Analysis ####

# uses ca package

?ca

caTree <- ca(tree)

caTree
plot(caTree, xlim = c(-0.5, 1), ylim = c(-0.5, 0.5))


# by hand
# Divide the data matrix caTree by the grand total of the matrix:

p <- as.matrix(tree/sum(tree))

# Cross tabulate row and column sums to be used in calculating expected
# values for the Chi Square values:

rs <- as.vector(apply(p, 1, sum))
cs <- as.vector(apply(p, 2, sum))

# Calculate expected values for the Chi Square calculation:
cp <- rs %*% t(cs)

# Calculate Chi Square values and check them out:
Qbar <- as.matrix((p - cp)/sqrt(cp))

# Conduct singular value decomposition (svd):

Q.svd <- svd(Qbar)


# Scale eigenvectors for rows and columns by the square root of row and
# column sums respectively:
V <- diag(1/sqrt(cs)) %*% Q.svd$v
Vhat <- diag(1/sqrt(rs)) %*% Q.svd$u

# Calculate ordination coordinates for both rows and columns:

F <- diag(1/rs) %*% p %*% V
Fhat <- diag(1/cs) %*% t(p) %*% Vhat


# Plot row and column coordinates in ordination space.

plot(Fhat[, 1:2], xlim = c(-0.5, 1), ylim = c(-0.5, 0.5), type = "n", xlab = "Coordinate 1", 
     ylab = "Coordinate 2", lwd = 2)
text(Fhat[, 1:2], labels = colnames(tree), cex = 0.7)
text(F[, 1:2], labels = rownames(tree), cex = 0.7)
abline(h = 0, lty = 2)
abline(v = 0, lty = 2)


rm(list=ls())

#### Portal data ####

library(portalr)
library(dplyr)
library(tidyr)

allPortal <- load_data(path = 'repo')

rodents <- allPortal[[1]]
trapping <- allPortal[[3]]
plots <- allPortal[[5]]

trapping <- trapping %>%
  select(month, year, period, plot, sampled)

plotsjoined <- left_join(plots, trapping, by = c('plot', 'month', 'year')) %>%
  filter(sampled == 1) %>%
  select(year, month, plot, treatment, period)

rodents <- rodents %>% 
  left_join(plotsjoined, by = c('year', 'month', 'plot', 'period')) %>%
  filter(treatment == 'control', year == 2016) %>%
  select(plot, species)


rodentCounts <- rodents %>%
  group_by(plot, species) %>%
  tally() %>%
  ungroup() %>%
  group_by(plot) %>%
  spread(key=species, value = n, fill = 0)

rodentBin <- rodentCounts
rodentBin[,2:17] <- ceiling(rodentBin[,2:17] / (rodentBin[,2:17] + 1))


# rodentBin is suitable for the same code as above (from lab). 
# rodentCounts is continuous.

