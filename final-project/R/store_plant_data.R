#' @importFrom magrittr "%>%"
#' @title store portal plant dataframes
#'
#' Import Portal rodent data using portalr functions.
#'
#' @export

library(dplyr)
store_plant_data <- function(){
primary_tables <- portalr::load_data(path = 'repo')

plot_treatments_summer <- primary_tables[[5]] %>%
  dplyr::filter(month %in% 4:9) %>%
  dplyr::select(year, plot, treatment) %>%
  dplyr::distinct()

summer <- portalr::plant_abundance(path = 'repo', level = 'Plot', type = 'Summer Annuals', plots = 'Longterm', unknowns = F, 
                          correct_sp = T, shape = 'flat', na_drop = T, effort = T)

summer <- summer %>% 
  dplyr::filter(season == 'summer') %>%
  dplyr::left_join(plot_treatments_summer, by = c('year', 'plot'))

winter <- portalr::plant_abundance(path = 'repo', level = 'Plot', type = 'Winter Annuals', plots = 'Longterm', unknowns = F, 
                             correct_sp = T, shape = 'flat', na_drop = T, effort = T)

plot_treatments_winter <- primary_tables[[5]] %>%
  dplyr::filter(month %in% c(1:3, 10:12)) %>%
  dplyr::select(year, plot, treatment) %>%
  dplyr::distinct()

winter <- winter %>% 
  dplyr::filter(season == 'winter') %>%
  dplyr::left_join(plot_treatments_winter, by = c('year', 'plot'))

summer_c <- summer %>%
  dplyr::filter(treatment == 'control')

winter_c <- winter %>%
  dplyr::filter(treatment == 'control')

summer_e <- summer %>%
  dplyr::filter(treatment == 'exclosure')

winter_e <- winter %>%
  dplyr::filter(treatment == 'exclosure')

write.csv(summer_c, 'final-project/data/summer-plants-raw-c.csv', row.names = F)
write.csv(winter_c, 'final-project/data/winter-plants-raw-c.csv', row.names = F)


write.csv(summer_e, 'final-project/data/summer-plants-raw-e.csv', row.names = F)
write.csv(winter_e, 'final-project/data/winter-plants-raw-e.csv', row.names = F)


}