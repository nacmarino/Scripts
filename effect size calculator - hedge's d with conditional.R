# function to calculate effect size and sampling variance for factorial meta-analysis
# this function is based on the paper by Gurevitch et al, 2000, Am Nat, The interaction between competition and predation: a meta-analysis of field experiments
# all the formulae are in the Appendix of the paper


# this function is used to calculate the correction factor for small sample bias from
# simple treatment and control differences
small_bias <- function(n_treatment, n_control) {
  1 - (3/((4*(n_control + n_treatment - 2)) - 1))
}


# this function is used to calculate the correction factor for small sample bias from
# the interaction of the four treatments
small_bias_interaction <- function(n_treatment.natural, n_control.natural, n_treatment.altered, n_control.altered) {
  1 - (3/((4*(n_treatment.natural + n_control.natural + n_treatment.altered + n_control.altered - 4)) - 1))
}


# this function calculates the weighted variance for each treatment
wgh_variance <- function(sample_size, std_deviation) {
  (sample_size - 1) * (std_deviation^2)
}


factorial_hedge <- function(x_natural = dados$mean_nopred_Controle, sd_natural = dados$error_nopred_Controle, n_natural = dados$n_nopred_Controle, 
                            x_natural_treatment = dados$mean_pred_Controle, sd_natural_treatment = dados$error_pred_Controle, n_natural_treatment = dados$n_pred_Controle,
                            x_altered = dados$mean_nopred_Tratamento, sd_altered = dados$error_nopred_Tratamento, n_altered = dados$n_nopred_Tratamento, 
                            x_altered_treatment = dados$mean_pred_Tratamento, sd_altered_treatment = dados$error_pred_Tratamento, n_altered_treatment = dados$n_pred_Tratamento, 
                            ID_column = dados$par_id, factorial = TRUE) { 
  
  
  
  if (factorial == FALSE) {
    sample_variance_natural <- sqrt((wgh_variance(n_natural, sd_natural) + wgh_variance(n_natural_treatment, sd_natural_treatment))/(n_natural + n_natural_treatment -2))
    
    sample_variance_altered <- sqrt((wgh_variance(n_altered, sd_altered) + wgh_variance(n_altered_treatment, sd_altered_treatment))/(n_altered + n_altered_treatment -2))
    
    
  }
  
  else {
    # pooled sample variance of the hedge's d effect size metric for all factors
    sample_variance_natural <- sqrt((wgh_variance(n_natural, sd_natural) + wgh_variance(n_natural_treatment, sd_natural_treatment) +
                                       wgh_variance(n_altered, sd_altered) + wgh_variance(n_altered_treatment, sd_altered_treatment))/
                                      (n_natural + n_natural_treatment + n_altered + n_altered_treatment - 4))
    
    sample_variance_altered <- sample_variance_natural
  }
  
  
  
  
  # calculations of the hedge's d effect size
  # hedge's d for the treatment under natural conditions
  d_treat.natural <- ((x_natural_treatment - x_natural)/sample_variance_natural) * small_bias(n_natural, n_natural_treatment)
  
  # hedge's d for the treatment under altered conditions
  d_treat.altered <- ((x_altered_treatment - x_altered)/sample_variance_altered) * small_bias(n_altered, n_altered_treatment)
  
  # this object is the sampling variance of adding a predator under the ambient condtion
  d_var_treat.natural <- (1/n_natural) + (1/n_natural_treatment) + ((d_treat.natural^2)/(2*(n_natural_treatment + n_natural)))
  
  # this object is the sampling variance of adding a predator under the climate change condition
  d_var_treat.altered <- (1/n_altered) + (1/n_altered_treatment) + ((d_treat.altered^2)/(2*(n_altered_treatment + n_altered)))
  
  # sample size treatment, natural
  size_n_treat.natural <- mean(c(n_natural, n_natural_treatment))
  
  # sample size treatment, altered
  size_n_treat.altered <- mean(c(n_altered, n_altered_treatment))
  
  
  
  # calculation for the overall effect of each factor in the experiment
  # effect size of the pure effect of the treatment
  d_treatment <- (((x_natural_treatment + x_altered_treatment) - (x_natural + x_altered))/(2*sample_variance_natural)) * small_bias_interaction(n_natural, n_natural_treatment, n_altered, n_altered_treatment)
  
  # effect size of the pure effect of the altered condition
  d_altered <- (((x_altered_treatment + x_altered) - (x_natural_treatment + x_natural))/(2*sample_variance_altered)) * small_bias_interaction(n_natural, n_natural_treatment, n_altered, n_altered_treatment)
  
  # variance of the pure effect of the treatment
  d_var_treatment <- ((1/n_altered) + (1/n_altered_treatment) + (1/n_natural) + (1/n_natural_treatment) + ((d_treatment^2)/(2*(n_altered + n_altered_treatment + n_natural + n_natural_treatment))))*(1/4)
  
  # variance of the pure effect of the altered condition
  d_var_altered <- ((1/n_altered) + (1/n_altered_treatment) + (1/n_natural) + (1/n_natural_treatment) + ((d_altered^2)/(2*(n_altered + n_altered_treatment + n_natural + n_natural_treatment))))*(1/4)
  
  
  
  
  
  
  
  # calculation for the effect size of the interaction according to hedge's d
  # this object is the result of adding a predator under climate change
  d_interaction <- ((d_treat.altered - d_treat.natural)/sample_variance_natural) * small_bias_interaction(n_natural, n_natural_treatment, n_altered, n_altered_treatment)
  
  # this object is the sampling variance of the interaction
  d_var_interaction<- (1/n_altered) + (1/n_altered_treatment) + (1/n_natural) + (1/n_natural_treatment) + ((d_interaction^2)/(2*(n_altered + n_altered_treatment + n_natural + n_natural_treatment)))
  
  # sample size for interaction
  sample_size <- mean(c(n_natural, n_natural_treatment, n_altered, n_altered_treatment))
  
  
  
  
  
  # output
  data.frame(par_id = ID_column, sample_variance_natural, sample_variance_altered, sample_size, 
             d_treat.natural, d_var_treat.natural, size_n_treat.natural,
             d_treat.altered, d_var_treat.altered, size_n_treat.altered, 
             d_treatment, d_var_treatment, 
             d_altered, d_var_altered, 
             d_interaction, d_var_interaction)
  
}