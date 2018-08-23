library(MVA)

# using usAir dataset

usAir <- USairpollution

usAir_mod <- read.csv("week-1/data/usAir_mod.csv", row = 1, header = TRUE)

dim(usAir)
str(usAir)

usAir[1,1]

usAir[1,]
usAir[,1]

usAir[-1,]
usAir[,-1]

usAir[1:5, ]
usAir[, 1:5]

usAir[-(1:5),]
usAir[,-(1:5)]


usAir[usAir[,2] < 50, ]

t(usAir)

matrix_usAir <- as.matrix(usAir$SO2)
t(matrix_usAir)

temp <- usAir$temp

ranks <- rank(temp)
sorted <- sort(temp)
ordered <- order(temp)

table <- data.frame(temp, ranks, sorted, ordered)

des_sort <- rev(sort(temp))

temp_ordered <- usAir[order(temp),]

mean(usAir[, 3])

median(usAir[, 3])

var(usAir[, 3])

sum(usAir[,3])

mean(t(usAir[3,]))
median(t(usAir[3, ]))
var(t(usAir[3,]))
sum(t(usAir[3,]))

columnSum <- colSums(usAir)

rowSums(usAir)
colMeans(usAir)
rowMeans(usAir)

write.csv(columnSum, "week-1/outputs/usAir_columnSums.csv")

rm(list=ls())

# using pitcher plant dataset

pitcher <- read.csv('week-1/data/Pitcher.csv', header = T, stringsAsFactors = F)


dim(pitcher)
str(pitcher)

pitcher[1,1]

pitcher[1,]
pitcher[,1]

pitcher[-1,]
pitcher[,-1]

pitcher[1:5, ]
pitcher[, 1:5]

pitcher[-(1:5),]
pitcher[,-(1:5)]


pitcher[pitcher[,2] < 50, ]

t(pitcher)

matrix_pitcher <- as.matrix(pitcher$NEMATO)
t(matrix_pitcher)

MOSQUI <- pitcher$MOSQUI

ranks <- rank(MOSQUI)
sorted <- sort(MOSQUI)
ordered <- order(MOSQUI)

table <- data.frame(MOSQUI, ranks, sorted, ordered)

des_sort <- rev(sort(MOSQUI))

MOSQUI_ordered <- pitcher[order(MOSQUI),]

mean(pitcher[, 3])

median(pitcher[, 3])

var(pitcher[, 3])

sum(pitcher[,3])

mean(t(pitcher[3, -1]))
median(t(pitcher[3, -1]))
var(t(pitcher[3, -1]))
sum(t(pitcher[3, -1]))

columnSum <- colSums(pitcher[,-1])

rowSums(pitcher[,-1])
colMeans(pitcher[,-1])
rowMeans(pitcher[,-1])

write.csv(columnSum, "week-1/outputs/pitcher_columnSums.csv")


