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

#### PCA or PCoA? ####
