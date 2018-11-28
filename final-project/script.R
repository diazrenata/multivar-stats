#### Data processing ####
source('final-project/R/store_plant_data.R')
source('final-project/R/store_rodent_data.R')
source('final-project/R/store_adjusted_data.R')

library(dplyr)

store_plant_data()
store_rodent_data()

rodents <- read.csv('final-project/data/summer-rodents-raw.csv', 
                    stringsAsFactors = F)
summer_plants <- read.csv('final-project/data/summer-plants-raw.csv', 
                          stringsAsFactors = F)
winter_plants <- read.csv('final-project/data/winter-plants-raw.csv', 
                          stringsAsFactors = F)
all_plants <- rbind(summer_plants, winter_plants)

store_adjusted_data(plant_data = all_plants, rodent_data = rodents,
                    focal_type = 'rodent')
store_adjusted_data(plant_data = summer_plants, focal_type = 'plant',
                    season = 'summer')
store_adjusted_data(plant_data = winter_plants, focal_type = 'plant',
                    season = 'winter')

rm(list=ls())

#### Data exploration ####
summer_plants <- read.csv('final-project/data/summer-plants-adjusted.csv', 
                          stringsAsFactors = F)

 # logColSums <- log(colSums(summer_plants))
 # hist(logColSums)
# 
# veryrare <- which(logColSums <= 2)
# rare <- which(logColSums <= 4)
# common <- which(logColSums >= 8)
# verycommon <- which(logColSums >= 10)

summer_plants <- select(summer_plants, -year)
species <- colnames(summer_plants)
species <- as.data.frame(species)
species$occurrences <- 0

for(i in 1:nrow(species)) {
  not0 <- as.integer(length(which(summer_plants[,i] >= 1)))
  species$occurrences[i] <- not0
}

hist(species$occurrences)


length(which(species$occurrences >= 5))
filter(species, occurrences >= 30)

species$maximum <- 0

for(i in 1:nrow(species)){
  species$maximum[i] <- max(summer_plants[,i], na.rm = T)
}

hist(species$maximum[ which(species$maximum < 50)])

filter(species, maximum < 25)

# OK so we're going to go with removing species whose maximum abundance
# in any year is less than 50. 

rm(list=ls())

#### PCoA ####
source('final-project/R/run_pcoa.R')
summer_plants <- read.csv('final-project/data/summer-plants-adjusted.csv', 
                          stringsAsFactors = F)

run_pcoa(summer_plants, season = 'Summer')

winter_plants <- read.csv('final-project/data/winter-plants-adjusted.csv', 
                          stringsAsFactors = F)

run_pcoa(winter_plants, season = 'Winter')

# looks like first 3 axes are the most useful ones?

rm(list=ls())

#### RDA #### 

rodents <- read.csv('final-project/data/rodents-adjusted.csv', 
                    stringsAsFactors = F)

summer_axes <- read.csv('final-project/models/Summer_pcoa_vals.csv', 
                        stringsAsFactors = F)
summer_axes <- summer_axes[,1:4]

winter_axes <- read.csv('final-project/models/Winter_pcoa_vals.csv',
                        stringsAsFactors = F)
winter_axes <- winter_axes[,1:4]

pred_vals <- inner_join(winter_axes, summer_axes, by = 'year')

rodents <- filter(rodents, year %in% pred_vals$year) %>%
  select(-year)
rodents_hel <- decostand(rodents, 'hellinger')

rodents.rda <- rda(rodents_hel ~ ., pred_vals)

# unadjusted and adjusted r2
R2 <- RsquareAdj(rodents.rda)$r.squared
R2adj <- RsquareAdj(rodents.rda)$adj.r.squared

# plot w f scores
plot(rodents.rda, scaling = 1, main = "Triplot RDA rodents.rda ~ pred_vals - scaling 1 - w a scores")
spe.sc <- scores(spe.rda, choices = 1:2, scaling = 1, display = "sp")
arrows(0, 0, spe.sc[,1], spe.sc[,2], length = 0, lty = 1, col = 'red')

# plot w z scores
plot(rodents.rda, scaling = 1, display = c("sp", "lc", "cn"), main = "Triplot RDA spe.hel ~ env2 - scaling 1 - lc scores")
arrows(0, 0, spe.sc[, 1], spe.sc[, 2], length = 0, lty = 1, col = "red")

anova(rodents.rda, step = 1000)
anova(rodents.rda, by = "axis", step = 1000)


# Reduce the number of variables for the most parsimonious model.

set.seed(11)
step.forward <- ordiR2step(rda(rodents_hel ~ 1, data = pred_vals), scope = formula(rodents.rda), 
                           R2scope = F, direction = "forward", pstep = 1000)

# Most parsimonious is Call: rodents_hel ~ Winter 1 + year + Summer 2 + Winter 2 

rod.rda.pars <- rda(rodents_hel ~ WinterPCoAxis_1 + year + SummerPCoAxis_2 + WinterPCoAxis_2 + SummerPCoAxis_1 , pred_vals)

R2p <- RsquareAdj(rod.rda.pars)$r.squared
R2adjp <- RsquareAdj(rod.rda.pars)$adj.r.squared

# plot w f scores
plot(rod.rda.pars, scaling = 1, main = "Triplot")
spe.sc.p <- scores(rod.rda.pars, choices = 1:2, scaling = 1, display = "sp")
arrows(0, 0, spe.sc.p[,1], spe.sc.p[,2], length = 0, lty = 1, col = 'red')

# plot w z scores
plot(rod.rda.pars, scaling = 1, display = c("sp", "lc", "cn"), main = "Triplot RDA spe.hel ~ alt + oxy + dbo, env2 - scaling 1 - lc scores")
arrows(0, 0, spe.sc.p[, 1], spe.sc.p[, 2], length = 0, lty = 1, col = "red")

anova(rod.rda.pars, step = 1000)
anova(rod.rda.pars, by = "axis", step = 1000)


# Partial RDA

partial.alt <- rda(rodents_hel ~ year + Condition(WinterPCoAxis_1  + SummerPCoAxis_2 + WinterPCoAxis_2 + SummerPCoAxis_1), data = pred_vals)
anova(partial.alt, step = 1000)

# Variance partitioning
??varpart
rod_part <- varpart(rodents_hel, ~WinterPCoAxis_1, ~SummerPCoAxis_2, ~WinterPCoAxis_2, ~SummerPCoAxis_1, data = pred_vals)
rod_part

plot(rod_part, digits = 2)
# WinterPCoAxis_1 has the largetst chunk (.4)
