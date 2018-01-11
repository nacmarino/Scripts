#' @title es_factorial
#' 
#' @description This function is used to determine the nature of a trophic interaction: positive, negative or none. It is on the standardized mean 
#' difference: Hedge's d. These calculations are based on the formulae presented in the book "Handbook of Meta-analysis in Ecology and Evolution", 
#' by Koricheva, Gurevitch & Mengersen; formulae may be found in Chapter 6 - "Effect sizes: conventional choices and calculations".

# internal functions ------------------------------------------------------------------------------------------------------------------
# this function is used to calculate the correction factor for small sample bias from simple treatment and control differences
small_bias_correction <- function(n1i, n0i) {
  1 - (3/((4*(n1i + n0i - 2)) - 1))
}


# this function is used to calculate the correction factor for small sample bias from the interaction of the four treatments
bias_correction_interaction <- function(n1i1, n1i0, n0i1, n0i0) {
  1 - (3/((4*(n1i1 + n1i0 + n0i1 + n0i0 - 4)) - 1))
}

# function call -----------------------------------------------------------

es_fatorial <- function(m1i1, m1i0, m0i1, m0i0, sd1i1, sd1i0, sd0i1, sd0i0, n1i1, n1i0, n0i1, n0i0, data, int_direction) {
  
  # pooled standard deviation -----------------------------------------------
  
  pooled_numerator <- ((n1i1 - 1) * (sd1i1^2)) + ((n1i0 - 1) * (sd1i0^2)) + ((n0i1 - 1) * (sd0i1^2)) + ((n0i0 - 1) * (sd0i0^2))
  pooled_denominator <- n1i1 + n1i0 + n0i1 + n0i0 - 4
  pooled_sd <- sqrt(pooled_numerator/pooled_denominator)
  
  # effect sizes ------------------------------------------------------------
  # effect sizes for each scenario ------------------------------------------------------------------------------------------------------
  
  ES.m0i <- ((m0i1 - m0i0)/pooled_sd) * small_bias_correction(n1i = n0i1, n0i = n0i0) # consumer on ambient scenario
  ES.m1i <- ((m1i1 - m1i0)/pooled_sd) * small_bias_correction(n1i = n1i1, n0i = n1i0) # consumer on climate change scenario
  
  # variance of the effect sizes for each scenario
  
  VES.m0i <- (1/n0i1) + (1/n0i1) + ((ES.m0i^2)/(2*(n0i1 + n0i0)))
  VES.m1i <- (1/n1i1) + (1/n1i1) + ((ES.m1i^2)/(2*(n1i1 + n1i0)))
  
  
  # effect sizes for main effects -------------------------------------------------------------------------------------------------------
  # of climate change on responses
  ES.m1 <- (((m1i1 + m1i0) - (m0i1 + m0i0))/(2*pooled_sd)) * bias_correction_interaction(n1i1 = n1i1, n1i0 = n1i0, n0i1 = n0i1, n0i0 = n0i0)
  
  # of consumers on responses
  ES.i1 <- (((m0i1 + m1i1) - (m0i0 + m1i0))/(2*pooled_sd)) * bias_correction_interaction(n1i1 = n1i1, n1i0 = n1i0, n0i1 = n0i1, n0i0 = n0i0)
  
  # variance of effect size of climate change on responses
  VES.m1 <- ((1/n1i1) + (1/n1i0) + (1/n0i1) + (1/n0i0) + ((ES.m1^2)/(2*(n1i1 + n1i0 + n0i1 + n0i0))))*(1/4)
  
  # variance of effect size of consumers on responses
  VES.i1 <- ((1/n1i1) + (1/n1i0) + (1/n0i1) + (1/n0i0) + ((ES.i1^2)/(2*(n1i1 + n1i0 + n0i1 + n0i0))))*(1/4)
  
  bci <- bias_correction_interaction(n1i1 = n1i1, n1i0 = n1i0, n0i1 = n0i1, n0i0 = n0i0)
  
  ES.interaction <- ifelse(int_direction != "inverteu", 
                           (((abs(m1i1 - m1i0) - abs(m0i1 - m0i0))/pooled_sd) * bci),
                           (((m1i1 - m1i0) - (m0i1 - m0i0))/pooled_sd) * bci)
  
  
  VES.interaction <-  (1/n1i1) + (1/n1i0) + (1/n0i1) + (1/n0i0) + ((ES.interaction^2)/(2*(n1i1 + n1i0 + n0i1 + n0i0)))
  
  final_data <- cbind(data, ES.amb = round(ES.m0i, digits = 4), VES.amb = round(VES.m0i, digits = 4), 
                      ES.mc = round(ES.m1i, digits = 4), VES.mc = round(VES.m1i, digits = 4),
                      ES.mc.only = round(ES.m1, digits = 4), VES.mc.only = round(VES.m1, digits = 4), 
                      ES.cons.only = round(ES.i1, digits = 4), VES.cons.only = round(VES.i1, digits = 4),
                      ES.interaction = round(ES.interaction, digits = 4), 
                      VES.interaction = round(VES.interaction, digits = 4))
  
  return(final_data)
}

#' @example 
#' library(dplyr)
#' source("R functions/find_interaction.R")
#' teste_data <- data.frame(amb_pres = c(5, 10, 3, 11, 3, 4, 15, 5, 2.2, 3, 10),
#'                         amb_aus = c(10, 5, 7, 4, 3.1, 12, 5, 5.2, 2.3, 10, 3),
#'                         erro_amb_pres = rep(1, 11),
#'                         erro_amb_aus = rep(1, 11),
#'                         n_amb_pres = rep(10, 11),
#'                         n_amb_aus = rep(10, 11),
#'                         mc_pres = c(1, 9, 11, 3, 3.1, 2, 4.2, 6, 20, 7, 8),
#'                         mc_aus = c(9, 1, 4, 7, 3, 2.1, 4, 13, 7, 10, 4),
#'                         erro_mc_pres = rep(1, 11),
#'                         erro_mc_aus = rep(1, 11),
#'                         n_mc_pres = rep(10, 11),
#'                         n_mc_aus = rep(10, 11))
#'
#'teste_data <- find_interaction(x_presence = teste_data$amb_pres, sd_presence = teste_data$erro_amb_pres, n_presence = teste_data$n_amb_pres,
#'                               x_absence = teste_data$amb_aus, sd_absence = teste_data$erro_amb_aus, n_absence = teste_data$n_amb_aus,
#'                               cutoff = 0, data = teste_data) %>%
#'  rename(int_ambient = int_type) %>%
#'  find_interaction(x_presence = teste_data$mc_pres, sd_presence = teste_data$erro_mc_pres, n_presence = teste_data$n_mc_pres,
#'                   x_absence = teste_data$mc_aus, sd_absence = teste_data$erro_mc_aus, n_absence = teste_data$n_mc_aus,
#'                   cutoff = 0, data = .) %>%
#'  rename(int_climate = int_type) %>%
#'  mutate(int_direction = ifelse(int_ambient == int_climate, "manteve",
#'                                ifelse((int_ambient == "negative" | int_ambient == "positive") & int_climate == "none", "eliminou",
#'                                       ifelse(int_ambient == "none" & (int_climate == "negative" | int_climate == "positive"), "criou", "inverteu"))))
#'
#'
#'es_fatorial(m1i1 = teste_data$mc_pres, m1i0 = teste_data$mc_aus, m0i1 = teste_data$amb_pres, m0i0 = teste_data$amb_aus,
#'            sd1i1 = teste_data$erro_mc_pres, sd1i0 = teste_data$erro_mc_aus, sd0i1 = teste_data$erro_amb_pres, sd0i0 = teste_data$erro_amb_aus,
#'            n1i1 = teste_data$n_mc_pres, n1i0 = teste_data$n_mc_aus, n0i1 = teste_data$n_amb_pres, n0i0 = teste_data$n_amb_aus,
#'            data = teste_data, int_direction = teste_data$int_direction)
