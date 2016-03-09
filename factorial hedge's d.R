# function to calculate effect size and sampling variance for factorial meta-analysis
# this function is based on the paper by Gurevitch et al, 2000, Am Nat, The interaction between competition and predation: a meta-analysis of field experiments
# all the formulae are in the Appendix of the paper


# this function is used to calculate the correction factor for small sample bias from
# simple Aplus and control differences
small_bias <- function(n_Aplus, n_control) {
  1 - (3/((4*(n_control + n_Aplus - 2)) - 1))
}


# this function is used to calculate the correction factor for small sample bias from
# the interaction of the four Apluss
small_bias_interaction <- function(n_Aplus.Bminus, n_control.Bminus, n_Aplus.Bplus, n_control.Bplus) {
  1 - (3/((4*(n_Aplus.Bminus + n_control.Bminus + n_Aplus.Bplus + n_control.Bplus - 4)) - 1))
}


# this function calculates the weighted variance for each Aplus
wgh_variance <- function(sample_size, std_deviation) {
  (sample_size - 1) * (std_deviation^2)
}


factorial_hedge <- function(x_Bminus = mean_Bminus_Aminus, sd_Bminus = error_Bminus_Aminus, n_Bminus = n_Bminus_Aminus, 
                    x_Bminus_Aplus = mean_Bminus_Aplus, sd_Bminus_Aplus = error_Bminus_Aplus, n_Bminus_Aplus = n_Bminus_Aplus,
                    x_Bplus = mean_Bplus_Aminus, sd_Bplus = error_Bplus_Aminus, n_Bplus = n_Bplus_Aminus, 
                    x_Bplus_Aplus = mean_Bplus_Aplus, sd_Bplus_Aplus = error_Bplus_Aplus, n_Bplus_Aplus = n_Bplus_Aplus, 
                    ID_column = par_id) { 
  

  
  # pooled sample variance of the hedge's d effect size metric for all factors
  sample_variance <- sqrt((wgh_variance(n_Bminus, sd_Bminus) + wgh_variance(n_Bminus_Aplus, sd_Bminus_Aplus) +
      wgh_variance(n_Bplus, sd_Bplus) + wgh_variance(n_Bplus_Aplus, sd_Bplus_Aplus))/
    (n_Bminus + n_Bminus_Aplus + n_Bplus + n_Bplus_Aplus - 4))
  
  
  
  # calculations of the hedge's d effect size
  # hedge's d for the Aplus under Bminus conditions
  d_Aplus.Bminus <- ((x_Bminus_Aplus - x_Bminus)/sample_variance) * small_bias(n_Bminus, n_Bminus_Aplus)
  
  # hedge's d for the Aplus under Bplus conditions
  d_Aplus.Bplus <- ((x_Bplus_Aplus - x_Bplus)/sample_variance) * small_bias(n_Bplus, n_Bplus_Aplus)
  
  # this object is the sampling variance of adding a predator under the ambient condtion
  d_var_Aplus.Bminus <- (1/n_Bminus) + (1/n_Bminus_Aplus) + ((d_Aplus.Bminus^2)/(2*(n_Bminus_Aplus + n_Bminus)))
  
  # this object is the sampling variance of adding a predator under the climate change condition
  d_var_Aplus.Bplus <- (1/n_Bplus) + (1/n_Bplus_Aplus) + ((d_Aplus.Bplus^2)/(2*(n_Bplus_Aplus + n_Bplus)))
  
  # sample size Aplus, Bminus
  size_n_Aplus.Bminus <- mean(c(n_Bminus, n_Bminus_Aplus))
  
  # sample size Aplus, Bplus
  size_n_Aplus.Bplus <- mean(c(n_Bplus, n_Bplus_Aplus))
  
  
  
  # calculation for the overall effect of each factor in the experiment
  # effect size of the pure effect of the Aplus
  d_Aplus <- (((x_Bminus_Aplus + x_Bplus_Aplus) - (x_Bminus + x_Bplus))/(2*sample_variance)) * small_bias_interaction(n_Bminus, n_Bminus_Aplus, n_Bplus, n_Bplus_Aplus)
  
  # effect size of the pure effect of the Bplus condition
  d_Bplus <- (((x_Bplus_Aplus + x_Bplus) - (x_Bminus_Aplus + x_Bminus))/(2*sample_variance)) * small_bias_interaction(n_Bminus, n_Bminus_Aplus, n_Bplus, n_Bplus_Aplus)
  
  # variance of the pure effect of the Aplus
  d_var_Aplus <- ((1/n_Bplus) + (1/n_Bplus_Aplus) + (1/n_Bminus) + (1/n_Bminus_Aplus) + ((d_Aplus^2)/(2*(n_Bplus + n_Bplus_Aplus + n_Bminus + n_Bminus_Aplus))))*(1/4)
  
  # variance of the pure effect of the Bplus condition
  d_var_Bplus <- ((1/n_Bplus) + (1/n_Bplus_Aplus) + (1/n_Bminus) + (1/n_Bminus_Aplus) + ((d_Bplus^2)/(2*(n_Bplus + n_Bplus_Aplus + n_Bminus + n_Bminus_Aplus))))*(1/4)
  
  
  
  
  
  
  
  # calculation for the effect size of the interaction according to hedge's d
  # this object is the result of adding a predator under climate change
  d_interaction <- ((d_Aplus.Bplus - d_Aplus.Bminus)/sample_variance) * small_bias_interaction(n_Bminus, n_Bminus_Aplus, n_Bplus, n_Bplus_Aplus)
  
  # this object is the sampling variance of the interaction
  d_var_interaction<- (1/n_Bplus) + (1/n_Bplus_Aplus) + (1/n_Bminus) + (1/n_Bminus_Aplus) + ((d_interaction^2)/(2*(n_Bplus + n_Bplus_Aplus + n_Bminus + n_Bminus_Aplus)))
  
  # sample size for interaction
  sample_size <- mean(c(n_Bminus, n_Bminus_Aplus, n_Bplus, n_Bplus_Aplus))
  
  
  
  
  
  # output
  data.frame(par_id = ID_column, sample_variance, sample_size, 
             d_Aplus.Bminus, d_var_Aplus.Bminus, size_n_Aplus.Bminus,
             d_Aplus.Bplus, d_var_Aplus.Bplus, size_n_Aplus.Bplus, 
             d_Aplus, d_var_Aplus, 
             d_Bplus, d_var_Bplus, 
             d_interaction, d_var_interaction)

}