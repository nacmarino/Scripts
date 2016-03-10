# a are the rows, b are the columns
# similarly, 0 is absence, 1 is presence

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
small_bias_interaction <- function(n_a0b1, n_a0b0, n_a1b1, n_a1b0) {
  1 - (3/((4*(n_a0b1 + n_a0b0 + n_a1b1 + n_a1b0 - 4)) - 1))
}


# this function calculates the weighted variance for each treatment
wgh_variance <- function(sample_size, std_deviation) {
  (sample_size - 1) * (std_deviation^2)
}


factorial_hedge <- function(x_a0 = mean_a0b0, sd_a0 = error_a0b0, n_a0 = n_a0b0, 
                            x_a0b1 = mean_a0b1, sd_a0b1 = error_a0b1, n_a0b1 = n_a0b1,
                            x_a1 = mean_a1b0, sd_a1 = error_a1b0, n_a1 = n_a1b0, 
                            x_a1b1 = mean_a1b1, sd_a1b1 = error_a1b1, n_a1b1 = n_a1b1, 
                            ID_column = par_id) { 
  
  
  
  # pooled sample variance of the hedge's d effect size metric for all factors
  sample_variance <- sqrt((wgh_variance(n_a0, sd_a0) + wgh_variance(n_a0_treatment, sd_a0_treatment) +
                             wgh_variance(n_a1, sd_a1) + wgh_variance(n_a1b1, sd_a1b1))/
                            (n_a0 + n_a0_treatment + n_a1 + n_a1b1 - 4))
  
  
  
  # calculations of the hedge's d effect size
  # hedge's d for the treatment under a0 conditions
  d_a0b1 <- ((x_a0b1 - x_a0)/sample_variance) * small_bias(n_a0, n_a0b1)
  
  # hedge's d for the treatment under a1 conditions
  d_a1b1 <- ((x_a1b1 - x_a1)/sample_variance) * small_bias(n_a1, n_a1b1)
  
  # this object is the sampling variance of adding a b1 under the a0 condition
  d_var_a0b1 <- (1/n_a0) + (1/n_a0b1) + ((d_a0b1^2)/(2*(n_a0b1 + n_a0)))
  
  # this object is the sampling variance of adding a b1 under the a1 condition
  d_var_a1b1 <- (1/n_a1) + (1/n_a1b1) + ((d_a1b1^2)/(2*(n_a1b1 + n_a1)))
  
  # sample size treatment, a0
  size_n_a0b1 <- mean(c(n_a0, n_a0b1))
  
  # sample size treatment, a1
  size_n_a1b1 <- mean(c(n_a1, n_a1b1))
  
  
  
  # calculation for the overall effect of each factor in the experiment
  # effect size of the pure effect of the treatment
  d_treatment <- (((x_a0b1 + x_a1b1) - (x_a0 + x_a1))/(2*sample_variance)) * small_bias_interaction(n_a0, n_a0b1, n_a1, n_a1b1)
  
  # effect size of the pure effect of the a1 condition
  d_a1 <- (((x_a1b1 + x_a1) - (x_a0b1 + x_a0))/(2*sample_variance)) * small_bias_interaction(n_a0, n_a0b1, n_a1, n_a1b1)
  
  # variance of the pure effect of the treatment
  d_var_treatment <- ((1/n_a1) + (1/n_a1b1) + (1/n_a0) + (1/n_a0b1) + ((d_treatment^2)/(2*(n_a1 + n_a1b1 + n_a0 + n_a0b1))))*(1/4)
  
  # variance of the pure effect of the a1 condition
  d_var_a1 <- ((1/n_a1) + (1/n_a1b1) + (1/n_a0) + (1/n_a0b1) + ((d_a1^2)/(2*(n_a1 + n_a1b1 + n_a0 + n_a0b1))))*(1/4)
  
  
  
  
  
  
  
  # calculation for the effect size of the interaction according to hedge's d
  # this object is the result of adding a b1 under a1
  d_interaction <- ((d_a1b1 - d_a0b1)/sample_variance) * small_bias_interaction(n_a0, n_a0b1, n_a1, n_a1b1)
  
  # this object is the sampling variance of the interaction
  d_var_interaction<- (1/n_a1) + (1/n_a1b1) + (1/n_a0) + (1/n_a0b1) + ((d_interaction^2)/(2*(n_a1 + n_a1b1 + n_a0 + n_a0b1)))
  
  # sample size for interaction
  sample_size <- mean(c(n_a0, n_a0b1, n_a1, n_a1b1))
  
  
  
  
  
  # output
  data.frame(par_id = ID_column, sample_variance, sample_size, 
             d_a0b1, d_var_a0b1, size_n_a0b1,
             d_a1b1, d_var_a1b1, size_n_a1b1, 
             d_treatment, d_var_treatment, 
             d_a1, d_var_a1, 
             d_interaction, d_var_interaction)
  
}