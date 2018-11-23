source('final-project/R/store_plant_data.R')
source('final-project/R/store_rodent_data.R')
source('final-project/R/store_adjusted_data.R')

library(dplyr)

store_plant_data()
store_rodent_data()

summer_plants <- read.csv('final-project/data/summer-plants-raw.csv', 
                   stringsAsFactors = F)
summer_rodents <- read.csv('final-project/data/summer-rodents-raw.csv',
                           stringsAsFactors = F)

store_adjusted_data(summer_rodents, summer_plants, 'summer')

summer_plants <- read.csv('final-project/data/summer-plants-adjusted.csv', 
                   stringsAsFactors = F)
summer_rodents <- read.csv('final-project/data/summer-rodents-adjusted.csv', 
                   stringsAsFactors = F)

library(LDATS)

source('final-project/R/run_LDA.R')

plant_comm <- summer_plants %>%
  select(-year)

plant_LDA <- run_LDA(plant_comm)
