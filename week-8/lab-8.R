library(MASS)
library(candisc)
library(ade4)
library(vegan)

# Testing assumptions of discriminant analysis
# homogeneity of variances - fligner-killeen tests
# significant differences --> reason to transform
fligner.test(iris$Sepal.Length, iris$Species)
fligner.test(iris$Sepal.Width, iris$Species)
fligner.test(iris$Petal.Length, iris$Species)
fligner.test(iris$Petal.Width, iris$Species)

# transform
log <- cbind.data.frame(apply(iris[, 1:4] + 1, 2, log), iris$Species)
names(log)[5] <- "Species"
fligner.test(log$Sepal.Length, log$Species)
fligner.test(log$Sepal.Width, log$Species)
fligner.test(log$Petal.Length, log$Species)
fligner.test(log$Petal.Width, log$Species)

# transformation got rid of the significantly not homogeneous variances

# Check for multicolinearity
# Correlations > 0.7 can be problematic.
cor(iris[, 1:4])
# There is some multicolinearity. Later we will re-run this but remove petal length.


# Check for outliers
# Calculate a withi-group distance matrix:

brayDist <- vegdist(iris[1:50, 1:4], "bray")

# Calculate the average distance of each sample to all other samples (i.e.
# column average) and turn the means in z-scores:

multOut <- scale(colMeans(as.matrix(brayDist)))

# Now look at a histogram of these data to identify samples that are > 3 sd
# from the mean:

hist(multOut)

# and get the number of those samples:

multOut[multOut > 3, ]

# Repeat for the other two groups (i.e.,species):

brayDist <- vegdist(iris[51:100, 1:4], "bray")
multOut <- scale(colMeans(as.matrix(brayDist)))
hist(multOut)
multOut[multOut > 3, ]

##      61 
## 3.07037
brayDist <- vegdist(iris[101:150, 1:4], "bray")
multOut <- scale(colMeans(as.matrix(brayDist)))
hist(multOut)

multOut[multOut > 3, ]
##      107 
## 3.816638
# Finally, make a vector of the outliers to pull out of the data set later:

Outliers <- c(42, 61, 107)


### with petal.length and outliers ####

# check that the variables change linearly along canonical axes
pairs(iris[, 1:4])
# these are, apparently, "linear"

# Discriminant analysis

# Create training and analysis datasets
set.seed(11)
train <- sample(1:150, 75)
prior <- table(iris$Sp[train])
iris.LDA <- lda(Species ~ ., iris, prior = cbind(prior/75), subset = train)

# lda in this context is 
# LINEAR DISCRIMINANT ANALYSIS

iris.LDA

# 99% of trace is in LD1. So 1 meaningful axis 


iris.LDA.p <- predict(iris.LDA, iris[train, ])

corTest <- lm(iris.LDA.p$x ~ iris$Sp[train])
summary(corTest)

plot(iris.LDA, xlim = c(-11, 11), ylim = c(-6, 6))



# Test the classification accuracy
ct <- table(iris[train, ]$Species, predict(iris.LDA, iris[train, ])$class)

# Change to a table of proportions:

pct <- prop.table(ct)
pct

# Calculate classification rate by summing the diagonal:

sum(diag(pct))

# Interpreting the canonical axes
# MASS::lda doesn't have good functionality for this
# So use candisc::candisc 

# Build a linear model using lm:

iris.mod <- lm(cbind(Sepal.Length, Sepal.Width, Petal.Length, Petal.Width) ~ 
                 Species, data = iris[train, ])

# candisc does a generalized canonical discriminant analysis
# for one term in a multivariate linear model

iris.can <- candisc(iris.mod, data = iris[train, ])


iris.can$coeffs.raw # same as coefficients from MASS above
iris.can$coeffs.std # describes variables' contributions to the axis
iris.can$structure # describes how well the axis describes the variable

# Use structure coefficients to interpret axes.

# Validate canonical axes - how well does the da classify 
# new data?

iris.LDA.new <- predict(iris.LDA, iris[-train, ])

ynew.table <- table(iris[-train, ]$Species, iris.LDA.new$class)
ynew.table

sum(diag(prop.table(ynew.table)))


# MANOVA
# DA provides groups, MANOVA tests whether the groups are different.
Y <- as.matrix(iris[, 1:4])
Sp <- factor(iris[, 5])

fit <- manova(Y ~ Sp)
summary(fit, test = "Wilks")
Yset <- as.matrix(iris[1:50, 1:4])
Yversi <- as.matrix(iris[51:100, 1:4])
Yvirg <- as.matrix(iris[101:150, 1:4])
Sp <- factor(iris[, 5])

fit1 <- manova(rbind(Yset, Yversi) ~ Sp[1:100])
summary(fit1, test = "Hotelling-Lawley")

fit2 <- manova(rbind(Yversi, Yvirg) ~ Sp[51:150])
summary(fit2, test = "Hotelling-Lawley")

fit3 <- manova(rbind(Yset, Yvirg) ~ Sp[-c(51:100)])
summary(fit3, test = "Hotelling-Lawley")



# all fits are significant



### pulling out petal.length and outliers ####

iris2 <- iris[-Outliers ,c(1,2,4,5)]

# check that the variables change linearly along canonical axes
pairs(iris2[, 1:3])
# these are, apparently, "linear"

# Discriminant analysis

# Create training and analysis datasets
set.seed(11)
train <- sample(1:nrow(iris2), 75)
prior <- table(iris2$Sp[train])
iris2.LDA <- lda(Species ~ ., iris2, prior = cbind(prior/75), subset = train)

# lda in this context is 
# LINEAR DISCRIMINANT ANALYSIS

iris2.LDA

# 99% of trace is in LD1. So 1 meaningful axis 


iris2.LDA.p <- predict(iris2.LDA, iris2[train, ])

corTest <- lm(iris2.LDA.p$x ~ iris2$Sp[train])
summary(corTest)

plot(iris2.LDA, xlim = c(-11, 11), ylim = c(-6, 6))



# Test the classification accuracy
ct <- table(iris2[train, ]$Species, predict(iris2.LDA, iris2[train, ])$class)

# Change to a table of proportions:

pct <- prop.table(ct)
pct

# Calculate classification rate by summing the diagonal:

sum(diag(pct))

# Interpreting the canonical axes
# MASS::lda doesn't have good functionality for this
# So use candisc::candisc 

# Build a linear model using lm:

iris2.mod <- lm(cbind(Sepal.Length, Sepal.Width,Petal.Width) ~ 
                 Species, data = iris2[train, ])

# candisc does a generalized canonical discriminant analysis
# for one term in a multivariate linear model

iris2.can <- candisc(iris2.mod, data = iris2[train, ])


iris2.can$coeffs.raw # same as coefficients from MASS above
iris2.can$coeffs.std # describes variables' contributions to the axis
iris2.can$structure # describes how well the axis describes the variable

# Use structure coefficients to interpret axes.

# Validate canonical axes - how well does the da classify 
# new data?

iris2.LDA.new <- predict(iris2.LDA, iris2[-train, ])

ynew.table <- table(iris2[-train, ]$Species, iris2.LDA.new$class)
ynew.table

sum(diag(prop.table(ynew.table)))


# MANOVA
# DA provides groups, MANOVA tests whether the groups are different.
Y <- as.matrix(iris2[, 1:3])
Sp <- factor(iris2[, 4])

fit <- manova(Y ~ Sp)
summary(fit, test = "Wilks")
Yset <- as.matrix(iris2[which(iris2$Species == 'setosa'), 1:3])
Yversi <- as.matrix(iris2[which(iris2$Species == 'versicolor'), 1:3])
Yvirg <- as.matrix(iris2[which(iris2$Species == 'virginica'), 1:3])
Sp <- factor(iris2[, 4])

fit1 <- manova(rbind(Yset, Yversi) ~ Sp[c(which(iris2$Species == 'setosa'), which(iris2$Species == 'versicolor'))])
summary(fit1, test = "Hotelling-Lawley")

fit2 <- manova(rbind(Yversi, Yvirg) ~ Sp[c(which(iris2$Species == 'versicolor'), which(iris2$Species == 'virginica'))])
summary(fit2, test = "Hotelling-Lawley")

fit3 <- manova(rbind(Yset, Yvirg) ~ Sp[c(which(iris2$Species == 'setosa'), which(iris2$Species == 'virginica'))])
summary(fit3, test = "Hotelling-Lawley")

# overall and all pairwise differences are significant
