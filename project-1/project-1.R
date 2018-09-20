library(dplyr)

darling <- read.csv("project-1/Darlingtonia.csv", header = T,
                    stringsAsFactors = F, row.names =1 )

# scale

darling <- scale(darling)

darling_pca <- princomp(darling, cor = F)
summary(darling_pca)

# latent root criterion
eigenVal <- (darling_pca$sdev * sqrt(nrow(darling)/(nrow(darling) - 1))) ^ 2

propVar <- eigenVal/sum(eigenVal)
cumVar <- cumsum(propVar)
pca_Table <- t(rbind(eigenVal, propVar, cumVar))
pca_Table

# latent root criterion (eliminate factors with eigenvalue < 1) 
# would indicate to keep 3 axes
# this is plausible given the scree plot:
# by axis 3 you have diminishing returns

plot(darling_pca, type = 'lines')

# look at loadings
loading_summary <- loadings(darling_pca)

loading_summary[which(loading_summary < .3)] <- NA

loading_summary[,1:3]

plot(darling_pca$loadings[,1:2],type="n",xlab="PC 1, 48%", ylab="PC 2, 21%",ylim=c(-1,1), xlim=c(-1,1))
text(darling_pca$loadings, labels=as.character(colnames(darling)), pos=1, cex=1)

# axis 1 seems to describe keel diam vs mouth diam, hood area

plot(darling_pca$scores,type="n",xlab="PC 1, 48%", ylab="PC 2, 21%",ylim=c(-8,8), xlim=c(-8,8))
text(darling_pca$scores, labels=as.character(rownames(darling)), pos=1, cex=1)

biplot(darling_pca$scores[c(63, 75),],darling_pca$loading,xlab="PC 1, 48%", ylab="PC 2, 21%",ylim=c(-8,8), xlim=c(-8,8))

biplot(darling_pca$scores,darling_pca$loading,xlabs= rep("*",87),xlab="PC 1, 34%", ylab="PC 2, 24%",ylim=c(-2,6), xlim=c(-4,7))
scores(darling_pca)
