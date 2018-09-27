library(dplyr)

#### Darlingtonia PCA #### 

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

biplot(darling_pca$scores,darling_pca$loading,xlabs= rep("*",87),xlab="PC 1, 48%", ylab="PC 2, 21%",ylim=c(-8,8), xlim=c(-8,8))
scores(darling_pca)

#### Birds NMDS ####

library(vegan)

birds <- read.csv('project-1/Atlantic_Caribbean.csv', row.names = 1)

birds_for_nmds <- birds[,4:384]
row.names(birds_for_nmds) <- birds[,3]

jbirds <- vegdist(birds_for_nmds, 'bray')


# Use metaMDS to get the NMDS 

nmdsBird <- metaMDS(birds_for_nmds, distance = "jaccard", k = 2, trace = T, trymax = 100000000000)
stressplot(nmdsBird)

nmdsBird

# High r2 (.99, .96) for the stress plot tells us that 
# the axes accurately preserve the rank-dissimilarlities
# in the original data.
# Min stress of 0.09; this is 'fair'. 
# Stress is more informative than r2?


# Comparing the historical and current communities


ordiplot(nmdsBird, type = 'n', xlim = c(-5, 5), ylim = c(-5, 5))
orditorp(nmdsBird, display = "sites", col = c(rep("green", 35), rep("blue", 19)), 
         air = 1, cex = .8)
legend(-4, -3.5, c("Atlantic", "Caribbean"), cex = 0.8, col = c("green", "blue"), 
       pch = 15:15)
title("Atlantic & Caribbean birds NMDS - colored by ocean")


arch_names <- birds %>%
  select(Archipelago) %>%
  distinct() %>%
  mutate(plot_col = rainbow(7))

plot_colors <- birds %>%
  select(Archipelago) %>%
  left_join(arch_names, by = 'Archipelago')


ordiplot(nmdsBird, type = 'n', xlim = c(-7, 7), ylim = c(-7, 7))
orditorp(nmdsBird, display = "sites", col = plot_colors$plot_col, 
         air = 1, cex = .8)
legend(-8, -1.5, arch_names$Archipelago, cex = 0.7, col = arch_names$plot_col, 
       pch = 15:15)
title("Atlantic & Caribbean birds NMDS - colored by archipelago")

