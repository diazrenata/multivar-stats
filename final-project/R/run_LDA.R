run_LDA <- function(community_data, topics_vector = 2:12, nseeds = 200, ncores = 4){
  
  #### Run LDAs ####
  LDA_models = LDATS::parLDA(data =community_data, ntopics =  topics_vector,
                             nseeds = nseeds, ncores = ncores)
  
  #### Select the best LDA (AICc) ####
  selected = LDATS:::LDA_select(lda_models = LDA_models, LDA_eval = quote(AIC), correction = TRUE,
                                LDA_selector = quote(min))
  # plot best lda topics (skip for now)
  # plot(selected)
  return(selected)
}