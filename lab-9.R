# Classification and regression trees

# using iris data

iris

# and ozone.csv

oz <- read.csv('week-9/ozone.csv', header = TRUE)

# packages
library(MASS)
library(rpart)
library(ade4)
library(vegan)


# CART
# Subset training and test data
set.seed(51)
train <- sample(1:150, 75)
freq <- table(iris$Sp[train])
freq

# Specify the model
# ~ . defines a model w/ all variables.
model <- Species ~ .

# Use rpart 
iris_rpart <- rpart(model, data = iris[train, ],
                    method = "class", control = rpart.control(minsplit = 10))
post(iris_rpart, file = "", title = "Iris Classification Tree")

summary(iris_rpart)



# Regression tree using ozone

