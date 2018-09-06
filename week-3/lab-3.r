library(MVA)
library(psych)
library(Hmisc)
library(vegan)
library(StatMatch)
library(MASS)
library(dplyr)

#### Primer of matrix algebra ####
newMatrix <- matrix(c(1,4,5,4,5,6,9,1,9), nrow = 3, ncol = 3)
newMatrix
dim(newMatrix)

oneMatrix <- matrix(c(1), nrow = 3, ncol = 3)
oneMatrix

# matrix addition
newMatrix + oneMatrix

# subtraction
newMatrix - oneMatrix

# scalar multiplication
3*newMatrix

# matrix multiplication
oneMatrix %*% newMatrix
newMatrix %*% oneMatrix

# matrix transposition
transMatrix <- t(newMatrix)
transMatrix

# identity matrix
identityMatrix <- diag(3)
identityMatrix

# matrix inversion
invMatrix <- solve(newMatrix)
invMatrix

invMatrix %*% newMatrix
round(invMatrix %*% newMatrix)

# eigenvalues & eigenvectors
eig <- eigen(newMatrix)
eig

#### PCA #### 
# using US air pollution data
usAir <- USairpollution
head(usAir)

hist(usAir$SO2) # not normal
hist(log(usAir$SO2)) # log is better but not perfect
hist(usAir$temp) # kind of normal
hist(log(usAir$SO2)) # log transformed looks worse
hist(usAir$manu) # not normal
hist(log(usAir$manu)) # log is better
hist(usAir$popul) # not normal
hist(log(usAir$popul)) # log is better
hist(usAir$wind) # kind of normal
hist(log(usAir$wind)) # transformed looks worse
hist(usAir$precip) # kind of normal
hist(log(usAir$precip)) # transformed doesn't look better
hist(usAir$predays) # kind of normal
hist(log(usAir$predays)) # transformed looks worse

# log transform SO2, manu, popul.

colnames(usAir)
usAir[,c(1, 3, 4)] <- log(usAir[,c(1,3,4)])

usAir_z <- scale(usAir)

usAir_pca <- princomp(usAir_z, cor = F)
summary(usAir_pca)

eigenVal <- (usAir_pca$sdev * sqrt(41/40)) ^ 2

propVar <- eigenVal/sum(eigenVal)
cumVar <- cumsum(propVar)
pca_Table <- t(rbind(eigenVal, propVar, cumVar))
pca_Table
loadings(usAir_pca)
scores(usAir_pca)

plot(usAir_pca, type = 'lines')

pca_Table


# factor loadings function
# function from lab doc
sigpca2<-function (x, permutations=1000, ...)
{
  pcnull <- princomp(x, ...)
  res <- pcnull$loadings
  out <- matrix(0, nrow=nrow(res), ncol=ncol(res))
  N <- nrow(x)
  for (i in 1:permutations) {
    pc <- princomp(x[sample(N, replace=TRUE), ], ...)
    pred <- predict(pc, newdata = x)
    r <-  cor(pcnull$scores, pred)
    k <- apply(abs(r), 2, which.max)
    reve <- sign(diag(r[k,]))
    sol <- pc$loadings[ ,k]
    sol <- sweep(sol, 2, reve, "*")
    out <- out + ifelse(res > 0, sol <=  0, sol >= 0)
  }
  out/permutations
}

sigpca2(usAir_z, permutations=1000)

pcnull<-princomp(usAir_z)
res <- pcnull$loadings
out <- matrix(0, nrow=nrow(res), ncol=ncol(res))
N <- nrow(usAir_z)
pc<-princomp(usAir_z[sample(N, replace=TRUE), ])
pred <- predict(pc, newdata = usAir_z)        
r <-  cor(pcnull$scores, pred)
k <- apply(abs(r), 2, which.max)
reve <- sign(diag(r[k,]))
sol <- pc$loadings[ ,k]
sol <- sweep(sol, 2, reve, "*")
out <- out + ifelse(res > 0, sol <=  0, sol >= 0)

# pca plot
summary(usAir_pca, loadings = T, cutoff = 0.3)

plot(usAir_pca$loadings[,1:2],type="n",xlab="PC 1, 34%", ylab="PC 2, 24%",ylim=c(-.8,.8), xlim=c(-.6,.6))
text(usAir_pca$loadings, labels=as.character(colnames(usAir_z)), pos=1, cex=1)

usAir_pca$loadings[,1:3]

# PC 1 corresponds to cold, densely populated, windy cities with 
# many precipitation days and lots of manufacturing and sulfur emissions.

# PC 2 corresponds to warm cities with high sulfur and precip days
# but low population density and manufacturing.

# pca plot with cities
plot(usAir_pca$scores,type="n",xlab="PC 1, 34%", ylab="PC 2, 24%",ylim=c(-4,4), xlim=c(-4,8))
text(usAir_pca$scores, labels=as.character(rownames(usAir_z)), pos=1, cex=1)

biplot(usAir_pca$scores,usAir_pca$loading,xlab="PC 1, 34%", ylab="PC 2, 24%",ylim=c(-2,6), xlim=c(-4,7))

biplot(usAir_pca$scores,usAir_pca$loading,xlabs= rep("*",41),xlab="PC 1, 34%", ylab="PC 2, 24%",ylim=c(-2,6), xlim=c(-4,7))

# eigen analysis
eig<-eigen(cov(usAir_z))
eig
eigVec<-as.matrix(eig$vectors[,1:2])
rownames(eigVec) <- rownames(cov(usAir_z))
eigVec


scores<-t(rbind(eigVec[,1]%*%t(usAir_z),eigVec[,2]%*%t(usAir_z)))#####hand calculated scores
biplot(scores,eigVec,xlab="PC 1, 34%", ylab="PC 2, 24%",ylim=c(-2,6), xlim=c(-4,7))


