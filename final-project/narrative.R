## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(dplyr)
library(vegan)
library(ca)

setwd('final-project')

## ----summer pcoa, echo = T-----------------------------------------------
summer_plants_c <- read.csv('data/summer-control-plants-adjusted.csv',
                          stringsAsFactors = F)

summer_plants_c_wis <- vegan::wisconsin(summer_plants_c[,2:ncol(summer_plants_c)])

summer_dist_mat <- vegdist(summer_plants_c_wis, 'bray')
  
summer_pcoa <- cmdscale(summer_dist_mat, k = nrow(summer_plants_c_wis) - 1, eig = T)


  
# Proportion of variance table 
  eigenvalues <- summer_pcoa$eig[1:nrow(summer_plants_c_wis)-1]
  propVar <- eigenvalues/sum(eigenvalues)
  cumVar <- cumsum(propVar)
  Summer_PCoA_Table <- cbind(eigenvalues, propVar, cumVar)
  
  Summer_PCoA_Table[1:15,]
  
  
# Scree plot:
  plot(eigenvalues)
  lines(lowess(eigenvalues))
  
  

ordiplot(scores(summer_pcoa)[, c(1, 2)], type = "n", cex = 1, main = "Summer plant PCoA")
## species scores not available
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)

# Add species

species_pc <- wascores(summer_pcoa$points[, 1:2], summer_plants_c_wis)
text(species_pc, rownames(species_pc), cex = 0.7, col = "red")
  

## ----save summer axes, echo = F------------------------------------------
summer_vals <- cbind(summer_plants_c$year, summer_pcoa$points[,1:3])
  
  season_names <- rep("SummerPCoAxis_", 3) %>%
    paste0(1:3)
  season_names <- c('year', season_names)
  colnames(summer_vals) <- season_names
  
  
  write.csv(summer_vals, 'models/summer_pcoa_vals_3.csv', row.names = F)
 

## ----winter pcoa, echo = T-----------------------------------------------
winter_plants_c <- read.csv('data/winter-control-plants-adjusted.csv',
                          stringsAsFactors = F)

winter_plants_c_wis <- vegan::wisconsin(winter_plants_c[,2:ncol(winter_plants_c)])


winter_dist_mat <- vegdist(winter_plants_c_wis, 'bray')
  
winter_pcoa <- cmdscale(winter_dist_mat, k = nrow(winter_plants_c_wis) - 1, eig = T)
  
# Proportion of variance table 
  eigenvalues <- winter_pcoa$eig[1:nrow(winter_plants_c_wis)-1]
  propVar <- eigenvalues/sum(eigenvalues)
  cumVar <- cumsum(propVar)
  winter_PCoA_Table <- cbind(eigenvalues, propVar, cumVar)
  
  winter_PCoA_Table[1:15,]
  
  
# Scree plot:
  plot(eigenvalues)
  lines(lowess(eigenvalues))
  

ordiplot(scores(winter_pcoa)[, c(1, 2)], type = "n", cex = 1, main = "Winter plant PCoA")
## species scores not available
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)

# Add species

species_pc <- wascores(winter_pcoa$points[, 1:2], winter_plants_c_wis)
text(species_pc, rownames(species_pc), cex = 0.7, col = "red")

## ----save winter axes, echo = F------------------------------------------
winter_vals <- cbind(winter_plants_c$year, winter_pcoa$points[,1:3])
  
  season_names <- rep("WinterPCoAxis_", 3) %>%
    paste0(1:3)
  season_names <- c('year', season_names)
  colnames(winter_vals) <- season_names
  
  
  write.csv(winter_vals, 'models/winter_pcoa_vals_3.csv', row.names = F)
 

#rm(list=ls())


## ----RDA, echo = T-------------------------------------------------------

rodents <- read.csv('data/rodents-adjusted.csv', 
                    stringsAsFactors = F)

summer_axes <- read.csv('models/Summer_pcoa_vals_3.csv', 
                        stringsAsFactors = F)

winter_axes <- read.csv('models/Winter_pcoa_vals_3.csv',
                        stringsAsFactors = F)


pred_vals <- inner_join(winter_axes, summer_axes, by = 'year')

rodents <- filter(rodents, year %in% pred_vals$year) %>%
  select(-year)

rodents_hel <- decostand(rodents, 'hellinger')

rodents_rda <- rda(rodents_hel ~ ., pred_vals)

R2 <- RsquareAdj(rodents_rda)$r.squared
R2adj <- RsquareAdj(rodents_rda)$adj.r.squared

R2
R2adj


anova(rodents_rda, step = 1000)
anova(rodents_rda, by = "axis", step = 1000)



## ----reduce vars, echo = T-----------------------------------------------

set.seed(11)
step.forward <- ordiR2step(rda(rodents_hel ~ 1, data = pred_vals), scope = formula(rodents_rda), 
                           R2scope = F, direction = "forward", pstep = 1000)

# Most parsimonious is Call: rodents_hel ~ WinterPCoAxis_1 + WinterPCoAxis_3 + WinterPCoAxis_2 + SummerPCoAxis_2

rod_rda_pars <- rda(rodents_hel ~ WinterPCoAxis_1 + WinterPCoAxis_3 + WinterPCoAxis_2 + SummerPCoAxis_2, pred_vals)

R2p <- RsquareAdj(rod_rda_pars)$r.squared
R2adjp <- RsquareAdj(rod_rda_pars)$adj.r.squared

R2p
R2adjp


anova(rod_rda_pars, step = 1000)
anova(rod_rda_pars, by = "axis", step = 1000)



## ----partial, echo = T---------------------------------------------------

rod_part <- varpart(rodents_hel, ~ WinterPCoAxis_1, ~ WinterPCoAxis_3, ~  WinterPCoAxis_2, ~ SummerPCoAxis_2, data = pred_vals)
rod_part

plot(rod_part, digits = 2)

## ----plot winter pcoa1 v year--------------------------------------------
 

plot(pred_vals$year, pred_vals$WinterPCoAxis_1)
abline(v = 1990, col = 'red')



colnames(rodents)

rodents_props <- rodents/rowSums(rodents)


plot(pred_vals$year, rodents_props$DS, col = 'blue', ylim = c(0, 1.1))
points(pred_vals$year, rodents_props$PE, col = 'green')
points(pred_vals$year, rodents_props$RM, col = 'purple')
points(pred_vals$year, rodents_props$DM, col = 'pink')
abline(v = 1990, col = 'red')



## ----winter plants exclosures, echo = T----------------------------------
winter_plants_e <- read.csv('data/winter-exclosure-plants-adjusted.csv',
                          stringsAsFactors = F)


winter_plants_e_wis <- vegan::wisconsin(winter_plants_e[,2:ncol(winter_plants_e)])

# which plants dominate winter axis 1?

impt_species <- as.data.frame(species_pc) %>%
  select(V1) %>%
  mutate(abs_score = abs(V1), species = row.names(species_pc)) %>%
  arrange(desc(abs_score))

valyears_c <- which(winter_plants_c$year %in% pred_vals$year)
valyears_e <- which(winter_plants_e$year %in% pred_vals$year)

winter_plants_c_plot <- winter_plants_c_wis[valyears_c, ]
winter_plants_e_plot <- winter_plants_e_wis[valyears_e, ]

par(mfrow=c(3,1))
plot(pred_vals$year, pred_vals$WinterPCoAxis_1)
plot(pred_vals$year, winter_plants_c_plot$eria.diff, col = 'red', ylim = c(0, .5))
points(pred_vals$year, winter_plants_c_plot$hapl.grac, col = 'red')
points(pred_vals$year, winter_plants_c_plot$esch.mexi, col = 'red')
points(pred_vals$year, winter_plants_c_plot$erig.drive, col = 'red')
points(pred_vals$year, winter_plants_c_plot$step.exig, col = 'red')
points(pred_vals$year, winter_plants_c_plot$pect.recu, col = 'blue')
points(pred_vals$year, winter_plants_c_plot$amsi.tess, col = 'blue')
points(pred_vals$year, winter_plants_c_plot$ambr.arte, col = 'blue')
points(pred_vals$year, winter_plants_c_plot$desc.obtu, col = 'blue')
points(pred_vals$year, winter_plants_c_plot$sisy.irio, col = 'blue')
points(pred_vals$year, winter_plants_c_plot$laen.coul, col = 'blue')
abline(v = 1990, col = 'red')
plot(pred_vals$year, winter_plants_e_plot$eria.diff, col = 'red', ylim = c(0, .5))
points(pred_vals$year, winter_plants_e_plot$hapl.grac, col = 'red')
points(pred_vals$year, winter_plants_e_plot$esch.mexi, col = 'red')
points(pred_vals$year, winter_plants_e_plot$erig.drive, col = 'red')
points(pred_vals$year, winter_plants_e_plot$step.exig, col = 'red')
points(pred_vals$year, winter_plants_e_plot$pect.recu, col = 'blue')
points(pred_vals$year, winter_plants_e_plot$amsi.tess, col = 'blue')
points(pred_vals$year, winter_plants_e_plot$ambr.arte, col = 'blue')
points(pred_vals$year, winter_plants_e_plot$desc.obtu, col = 'blue')
points(pred_vals$year, winter_plants_e_plot$sisy.irio, col = 'blue')
points(pred_vals$year, winter_plants_e_plot$laen.coul, col = 'blue')
abline(v = 1990, col = 'red')

# it's not exactly a slam dunk