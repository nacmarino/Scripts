#' @name randomize.community.abundance
#' 
#' @title Shuffles individuals from each species across sites in proportion to their total abundance and site marginal sum
#'
#' @description This function is used to randomize a community matrix containing species abundance data, preserving the total number of 
#' individuals found in each sample as well as species' total abundance across all samples.
#' 
#' 
## FUNCTION BODY
randomize.community.abundance <- function(community_data) {
  
  # loading packages, if needed
  if(!"tibble" %in% loadedNamespaces()) {
    suppressPackageStartupMessages(suppressMessages(library(tidyverse)))
  }
  
  # transform community data on a large format into long format, with each record of a individual from a given species as a separate row
  community_long <- suppressWarnings(
    community_data %>% 
      rownames_to_column(var = "site") %>% 
      gather(key = "species", value = "abundance", 2:ncol(.)) %>% 
      filter(abundance > 0) %>% 
      split(.$site) %>% 
      map(~ select(., -site)) %>% 
      map(~ data.frame(species = rep(.$species, .$abundance))) %>% 
      bind_rows(.id = "site")
  )
  
  # number of individual records to be randomized
  number_of_observations <- nrow(community_long)
  
  # randomizing species occurrences across sites
  ## each bromeliad has x slots, that will be repopulated by species according to their total abundance across all sites
  random_community <- data.frame(site = sample(x = community_long$site, size = number_of_observations, replace = FALSE),
                                 species = sample(x = community_long$species, size = number_of_observations, replace = FALSE)) %>% 
    group_by(site, species) %>% 
    summarise(abundance = n()) %>% 
    ungroup %>% 
    spread(key = species, value = abundance, fill = 0) %>% 
    data.frame %>% 
    `rownames<-`(.$site) %>% 
    select(-site)
  
  return(random_community)
  
}
#' 
#' @usage 
#' randomize.community.abundance(community_data)
#' 
#' @param community_data a species by site community data.frame or matrix containing species' abundances, with sites as rows and species as columns. Species abundances are expected to be integers or rounded.
#' 
#' @details 
#' This function starts by taking a species by site community data.frame or matrix in a large format and coverting it into a long format. Here, the data in the long format is also modified, in order that 
#' each individual from a given species at a site occupies on row in the new data.frame/matrix. In other words, the first step is to create a data.frame or matrix containing the records of each individual
#' organism in each site it occurs.
#' 
#' In the next step, the function randomized the sequence of sites that will be repopulated and then shuffles the individuals of each species in proportion to their total abundance across all sites. That is,
#' it fixes the number of individuals each site can support, rearrange the order through which individuals will be allocated to each site, and then shuffles the individuals across these sites.
#' 
#' Finally, the function recalculates species' abundances at each site and recreates the community data.frame/matrix in the large format.
#' 
#' @return a data.frame with the same number of sites (rows) and species (columns) as the original input, with the same marginal sums for the rows (total abundance per site)
#' and columns (species' total abundance), with species' abundances  distributed differently across samples.
#' 
#' @example 
#'
#' # creating a random community with sites as rows and species as columns
#' set.seed(33)
#' foo <- data.frame(sp1 = rnbinom(n = 9, size = 10, prob = 0.8),
#'                  sp2 = rbinom(n = 9, size = 100, prob = 0.2),
#'                  sp3 = rnbinom(n = 9, size = 20, prob = 0.5),
#'                  sp4 = rnbinom(n = 9, size = 2, prob = 0.5),
#'                  sp5 = rbinom(n = 9, size = 4, prob = 0.75),
#'                  sp6 = rnbinom(n = 9, size = 5, prob = 0.4),
#'                  sp7 = rbinom(n = 9, size = 1, prob = 0.5))
#' rownames(foo) <- paste0("site", rownames(foo))
#' foo
#'
#' # randomizing original community
#' set.seed(99)
#' foo1 <- randomize.community.abundance(foo)
#' foo1
#' 
#' # total abundance per species
#' colSums(foo1);colSums(foo)
#' 
#' # total abundance per site 
#' rowSums(foo); rowSums(foo1)