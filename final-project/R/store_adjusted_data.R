
store_adjusted_data <- function(rodent_data, plant_data, season){ 
  rodent_data <- rodent_data %>%
    select(-month, -treatment) %>%
    filter(year %in% plant_data$year)
  
  rodent_data_effort <- rodent_data %>%
    select(year, period, plot, sampled) %>%
    distinct() %>%
    group_by(year, period) %>%
    summarize(effort = sum(sampled)) %>%
    ungroup()
  
  rodent_data_adjusted <- rodent_data %>%
    select(-sampled) %>%
    group_by(year, period, species) %>%
    summarize(n = sum(n)) %>% 
    ungroup() %>%
    left_join(rodent_data_effort, 
              by = c("year", "period")) %>%
    mutate(adjusted_n =as.integer(ceiling((n / effort) * 8))) %>%
    select(year, species, adjusted_n) %>%
    group_by(year, species) %>%
    summarize(adjusted_total = sum(adjusted_n)) %>%
    ungroup() %>%
    tidyr::spread('species', 'adjusted_total', fill = 0)
  
  plant_data_effort <- plant_data %>%
    select(year, plot, quads) %>%
    distinct() %>%
    group_by(year) %>%
    summarize(total_quads = sum(quads)) %>%
    ungroup()
  
  plant_data_adjusted <- plant_data %>%
    select(-treatment, -season) %>%
    group_by(year, species) %>%
    summarize(total_n = sum(abundance)) %>%
    ungroup() %>%
    left_join(plant_data_effort, by = 'year') %>%
    mutate(adjusted_n = as.integer(ceiling((total_n / total_quads) * 64))) %>%
    select(year, species, adjusted_n) %>%
    tidyr::spread(species, adjusted_n, fill = 0)
  
  
  write.csv(rodent_data_adjusted, paste0('final-project/data/', season,
                                         '-rodents-adjusted.csv'),
            row.names = F)
  write.csv(plant_data_adjusted, paste0('final-project/data/', season,
                                         '-plants-adjusted.csv'),
            row.names = F)
  
}