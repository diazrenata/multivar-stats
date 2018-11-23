#' @importFrom magrittr "%>%"
#' @title store portal plant dataframes
#'
#' Import Portal rodent data using portalr functions.
#'
#' @export

library(dplyr)
store_rodent_data <- function(){
  primary_tables <- portalr::load_data(path = 'repo')
  
  rodents <- primary_tables[[1]]
  species <- primary_tables[[2]]
  trapping <- primary_tables[[3]]
  newmoons <- primary_tables[[4]]
  plots <- primary_tables[[5]]
  
  plots_history <- left_join(plots, trapping, 
                     by = c('year', 'month', 'plot')) 
  plots_history <- plots_history %>%
    select(year, month, plot, treatment, period, sampled)
  
  focal_spp <- species %>%
    filter(granivore == 1, unidentified == 0) %>%
    select(species)
  
  rodents_control <- rodents %>%
    filter(species %in% focal_spp$species) %>%
    select(month, year, period, plot, species) %>%
    group_by(month, year, period, plot, species) %>%
    summarize(n = n()) %>%
    ungroup() %>%
    left_join(plots_history, by = c('year', 'month', 'period', 'plot')) %>%
    filter(treatment == 'control')
  
  rodents_summer <- rodents_control %>%
    filter(month %in% c(4:9))
  
  rodents_winter <- rodents_control %>%
    filter(month %in% c(1:3, 10:12))
  
  write.csv(rodents_summer, 'final-project/data/summer-rodents-raw.csv', row.names = F)
  write.csv(rodents_winter, 'final-project/data/winter-rodents-raw.csv', row.names = F)
  
}