
## FUNCTION TO RANDOMIZE SPECIES' ABUNDANCES ACROSS SITES

make.random.community <- function(community_data) {
  
  # number of individual records to be randomized
  number_of_observations <- nrow(community_data)
  
  # randomizing species occurrences across bromeliads
  ## each bromeliad has x slots, that will be repopulated by invertebrate species according to their 
  ## total abundance across all bromeliads
  random_community <- data.frame(site = sample(x = community_data$site, size = number_of_observations, replace = FALSE),
                                 species = sample(x = community_data$species, size = number_of_observations, replace = FALSE)) %>% 
    group_by(site, species) %>% 
    summarise(abundance = n()) %>% 
    ungroup %>% 
    spread(key = species, value = abundance, fill = 0) %>% 
    data.frame %>% 
    `rownames<-`(.$site) %>% 
    select(-site)
  
  return(random_community)
  
}