library(raster)
library(vegan)

source('week-10/biostats.r')


# call in data
data(varespec)

data(varechem)

str(varespec)
summary(varespec)

str(varechem)
summary(varechem)

# select species: identify and remove extremely common or rare species
occur <- foa.plots(varespec)
# ^^ this is a pretty fun plotting function!

# remove the rare/common ones
rare <- which(occur[, 2] < 5)

common <- which(occur[, 2] > 95)

reduced <- varespec[, -c(rare, common)]

# are species' abundances normally distributed?
mapply(hist, as.data.frame(varespec[, 1:44]), main = colnames(varespec[, 1:44]), 
       xlab = "abundance")

# ..... n o p e --> log transform + 1

log.full <- log1p(varespec)
log.red <- log1p(reduced)

# check rowsum
rsum <- rowSums(log.full)
csum <- colSums(log.full)
hist(rsum)
hist(csum)

cv(rsum)
## [1] 15.6402
cv(csum)
## [1] 154.4498
# Reduced data set:
rsumRed <- rowSums(log.red)
csumRed <- colSums(log.red)
hist(rsumRed)
hist(csumRed)

cv(rsumRed)
cv(csumRed)

# If either the row or column sums have cv >50, standardize by the total:
cSpec <- sweep(log.full, 2, csum, "/")
cSpecRed <- sweep(log.red, 2, csumRed, "/")


# determine if species abundances show a linear (RDA) or 
# a unimodal (CCA) relationship with the underlying gradient.

# use Detrended Correspondence Analysis (DCA) to determine the length 
# of the canonical axes. You will use the decorana function 
# focus on axis length (1st axis). An axis length > 3 is evidence of a unimodal 
# relationship. An axis length of <3 is evidence of a linear relationship. 

decorana(cSpec) 
decorana(cSpecRed)

# axis > 3 --> use CCA


Vars <- varechem[, c(1, 2, 7)]
env <- as.data.frame(scale(Vars))
# Run CCA:
  sp.CCA <- cca(cSpec ~ ., data = env)

# Function for plotting species abundances vs. CCA Axis 1:
  f2 <- function(x) {
    plot(x ~ sp.CCA$CC$wa[, 1], xlab = "CCA AXIS 1", ylab = "Abundance ")
  }

# Apply the function across all the species:

mapply(f2, varespec)

# EXPLANATORY VARIABLES
# avoid multicolinearity
# a priori knowledge: focus on AL, P, N

# pairwise correlations?
Vars <- varechem[, c(1, 2, 7)]
Vars
round(as.dist(cor(Vars)), 2)

# check for high CV (indicates a need for scale)
cv(colSums(Vars))
## [1] 91.17077
# scale
env <- as.data.frame(scale(Vars))

# Before running the constrained model, 
# run an unconstrained ordination 
# (i.e. a regular Correspondence Analysis) 
# a measure of the amount of variation in the site 
# by species matrix that you will try to explain 
# with the explanatory variables (i.e. constraints). 

# Full Data
ca <- cca(cSpec)
plot(ca)
summary(ca)

# Reduced Data
ca <- cca(cSpecRed)
plot(ca)
summary(ca)


# The permutation allows you to test if your constrained axes 
# explain more variation than would be expected randomly. 
# Yse the anova.cca function in vegan to conduct the permutation.
# It is “anova-like” but not an anova. 
# Global Test (i.e. all variables together):
anova(sp.CCA)

# Axes Tests (i.e. each axis individually):
anova(sp.CCA, by = "axis")

# Variable Tests (i.e. each variable individually):
anova(sp.CCA, by = "terms")

# Observed (F matrix) and predicted (Z matrix) site scores
summary(sp.CCA)

# The matrix labeled “Site scores (weighted averages of species scores)” is
# the F matrix and 
# the matrix labeled “Site constraints (linear combinations of 
# constraining variables)”is the Z matrix. L
# ook at these two sets of site scores projected in ordination space:
  par(mfrow = c(1, 2))
plot(sp.CCA$CC$wa[, 1], sp.CCA$CC$wa[, 2], xlab = "CCA AXIS 1", ylab = "CCA AXIS 2")
plot(sp.CCA$CC$u[, 1], sp.CCA$CC$u[, 2], xlab = "CCA AXIS 1", ylab = "CCA AXIS 2")

# Look at the correlation between these two matrices. 
# These correlations can lend insight as to how well the 
# predicted site locations match the observed ones. 
# However, they are not to be trusted as the only line of evidence.
spenvcor(sp.CCA)

# Correlations between the Z matrix (predicted site scores) and 
# the environmental variables provide information on which 
# variables have the largest influence on the constrained ordination. 
# These also denote the placement of the environmental variables as 
# vectors on the CCA tri-plot.

sp.CCA$CCA$biplot
##The Tri-Plot (using the site scores from the F matrix)
plot(sp.CCA, choices = c(1, 2), display = c("wa", "sp", "bp"), scaling = 2)
# and using the site scores from the Z matrix:
plot(sp.CCA, choices = c(1, 2), display = c("lc", "sp", "bp"), scaling = 2)



# Rerunning using reduced dataset (remove common and rare species)

# Run CCA:
sp.CCA.reduced <- cca(cSpecRed ~ ., data = env)
# Global Test (i.e. all variables together):
anova(sp.CCA.reduced)

# Axes Tests (i.e. each axis individually):
anova(sp.CCA.reduced, by = "axis")

# Variable Tests (i.e. each variable individually):
anova(sp.CCA.reduced, by = "terms")

# Observed (F matrix) and predicted (Z matrix) site scores
summary(sp.CCA.reduced)

# The matrix labeled “Site scores (weighted averages of species scores)” is
# the F matrix and 
# the matrix labeled “Site constraints (linear combinations of 
# constraining variables)”is the Z matrix. L
# ook at these two sets of site scores projected in ordination space:
par(mfrow = c(1, 2))
plot(sp.CCA.reduced$CC$wa[, 1], sp.CCA.reduced$CC$wa[, 2], xlab = "CCA AXIS 1", ylab = "CCA AXIS 2")
plot(sp.CCA.reduced$CC$u[, 1], sp.CCA.reduced$CC$u[, 2], xlab = "CCA AXIS 1", ylab = "CCA AXIS 2")

# Look at the correlation between these two matrices. 
# These correlations can lend insight as to how well the 
# predicted site locations match the observed ones. 
# However, they are not to be trusted as the only line of evidence.
spenvcor(sp.CCA.reduced)

# Correlations between the Z matrix (predicted site scores) and 
# the environmental variables provide information on which 
# variables have the largest influence on the constrained ordination. 
# These also denote the placement of the environmental variables as 
# vectors on the CCA tri-plot.

sp.CCA.reduced$CCA$biplot
##The Tri-Plot (using the site scores from the F matrix)
plot(sp.CCA.reduced, choices = c(1, 2), display = c("wa", "sp", "bp"), scaling = 2)
# and using the site scores from the Z matrix:
plot(sp.CCA.reduced, choices = c(1, 2), display = c("lc", "sp", "bp"), scaling = 2)




