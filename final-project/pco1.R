library(dplyr)
controls_scores <-read.csv('final-project/models/winter_controls_species_scores.csv',
                           stringsAsFactors = F)
exclosures_scores <- read.csv('final-project/models/winter_exclosures_species_scores.csv',
                              stringsAsFactors = F)

colnames(controls_scores) <- c('species', 'cpco1', 'pc02')
colnames(exclosures_scores) <- c('species', 'epco1', 'pc02')

exclosures_scores <- select(exclosures_scores, species, epco1)

pco1 <- controls_scores %>%
  select(species, cpco1) %>%
  inner_join(exclosures_scores, by = 'species')

plot(pco1$epco1, pco1$epco1, type = 'l')
points(pco1$epco1,pco1$cpco1)


controls_impt <- controls_scores %>%
  select(species, cpco1) %>%
  mutate(c_abs_score = abs(cpco1)) %>%
  arrange(desc(c_abs_score))


exclosures_impt <- exclosures_scores %>%
  select(species, epco1) %>%
  mutate(e_abs_score = abs(epco1)) %>%
  arrange(desc(e_abs_score))

head(controls_impt)
head(exclosures_impt)

impt <- full_join(exclosures_impt, controls_impt, by = 'species')
