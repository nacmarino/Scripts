# functions to calculate the effect size as the log response ratio
# the calculations performed here are based on the formulations presented on:
# Lajeunesse, 2011, Ecology, On the meta-analysis of response ratios for studies with correlated and multi-group designs
# Lajeunesse, 2015, Ecology, Bias and correction for the log response ratio in ecological meta-analysis


# this function takes two numeric vectors as the input, one of the vectors are the means
# of the treatment group and the other are the mean of a control group
# the function then takes the log ratio of the two vectors and return it to the user
# the definition of treatment and control depends on the research question:
# treatment is the factor you are interested on
make_RR <- function(x_treatment, x_control) {
  log(x_treatment/x_control)
}


# this function calculates the weighted variance for a given treatment
# it is used to calculate the individual weighted variance that will be used in the
# formation for the variance of the RR
make_weight_var <- function(x, std_dev, sample_size){
  (std_dev^2)/(sample_size*(x^2))
}


# this function calculates the bias correction for the weighted variance for a given treatment
# it is used to calculate the bias corrected individual weighted variance that will be used in the
# formation for the bias corrected variance of the RR
make_corrweight_var <- function(x, std_dev, sample_size){
  (std_dev^4)/((sample_size^2)*(x^4))
}




factorial_RR <- function(x_Bminus = dados$mean_B_minus.A_minus, sd_Bminus = dados$error_B_minus.A_minus, n_Bminus = dados$n_B_minus.A_minus, 
                            x_Bminus_Aplus = dados$mean_B_minus.A_plus, sd_Bminus_Aplus = dados$error_B_minus.A_plus, n_Bminus_Aplus = dados$n_B_minus.A_plus,
                            x_Bplus = dados$mean_B_plus.A_minus, sd_Bplus = dados$error_B_plus.A_minus, n_Bplus = dados$n_B_plus.A_minus, 
                            x_Bplus_Aplus = dados$mean_B_plus.A_plus, sd_Bplus_Aplus = dados$error_B_plus.A_plus, n_Bplus_Aplus = dados$n_B_plus.A_plus, 
                            ID_column = dados$par_id) { 
  
  
  
  # objects that will not go into the output of the function, but are middle steps to calculate other things
  # starts here
  # weighted variance for each factor of the factorial meta-analysis
  wghvar_Bminus.Aplus <- make_weight_var(x_Bminus_Aplus, sd_Bminus_Aplus, n_Bminus_Aplus)
  wghvar_Bminus.Aminus <- make_weight_var(x_Bminus, sd_Bminus, n_Bminus)
  wghvar_Bplus.Aplus <- make_weight_var(x_Bplus_Aplus, sd_Bplus_Aplus, n_Bplus_Aplus)
  wghvar_Bplus.Aminus <- make_weight_var(x_Bplus, sd_Bplus, n_Bplus)

  # bias correction for the individual weighted variance of the response ratio
  corrwghvar_Bminus.Aplus <- make_corrweight_var(x_Bminus_Aplus, sd_Bminus_Aplus, n_Bminus_Aplus)
  corrwghvar_Bminus.Aminus <- make_corrweight_var(x_Bminus, sd_Bminus, n_Bminus)
  corrwghvar_Bplus.Aplus <- make_corrweight_var(x_Bplus_Aplus, sd_Bplus_Aplus, n_Bplus_Aplus)
  corrwghvar_Bplus.Aminus <- make_corrweight_var(x_Bplus, sd_Bplus, n_Bplus)
  
  # ends here
  
  
  
  
  # calculation of the log response ratio and its variance
  # response ratio for the treatment under natural condition
  RR_Bminus.Aplus <- make_RR(x_Bminus_Aplus, x_Bminus)
  
  # response ratio for the treatment under altered condition
  RR_Bplus.Aplus <- make_RR(x_Bplus_Aplus, x_Bplus)
  
  # variance of the response ratio for the treatment under natural condition  
  varRR_Bminus.Aplus <-  wghvar_Bminus.Aplus + wghvar_Bminus.Aminus
  
  # variance of the response ratio for the treatment under altered condition
  varRR_Bplus.Aplus <- wghvar_Bplus.Aplus + wghvar_Bplus.Aminus  
  
  # sample size under natural condition
  size_n_Bminus <- (n_Bminus + n_Bminus_Aplus)/2
  
  # sample size under altered condition
  size_n_Bplus <- (n_Bplus + n_Bplus_Aplus)/2
  
  
  
  
  
  # calculation of the bias corrected log response ratio and its variance
  # bias corrected response ratio for the treatment under natural condition
  corrected_RR_Bminus.Aplus <- RR_Bminus.Aplus + 1/2*(wghvar_Bminus.Aplus - wghvar_Bminus.Aminus)
  
  # bias corrected response ratio for the treatment under altered condition
  corrected_RR_Bplus.Aplus <- RR_Bplus.Aplus + 1/2*(wghvar_Bplus.Aplus - wghvar_Bplus.Aminus)
  
  # bias corrected variance for the response ratio of the treatment under natural condition
  corrected_varRR_Bminus.Aplus <- varRR_Bminus.Aplus + 1/2*(corrwghvar_Bminus.Aplus + corrwghvar_Bminus.Aminus)
  
  # bias corrected variance for the response ratio of the treatment under altered condition
  corrected_varRR_Bplus.Aplus <- varRR_Bplus.Aplus + 1/2*(corrwghvar_Bplus.Aplus + corrwghvar_Bplus.Aminus)
  
  
  
  
  
  # calculation of the log response ratio for the interaction term between treatment and condition
  # log response ratio for the interaction
  RR_interaction <- RR_Bplus.Aplus - RR_Bminus.Aplus
  
  # bias correction of the log response ratio for the interaction
  corrected_RR_interaction <- corrected_RR_Bplus.Aplus - corrected_RR_Bminus.Aplus
  
  # variance of the log response ratio for the interaction
  varRR_interaction <- varRR_Bminus.Aplus + varRR_Bplus.Aplus
  
  # bias corrected variance of the log response ratio for the interaction
  corrected_varRR_interaction <-  corrected_varRR_Bminus.Aplus +  corrected_varRR_Bplus.Aplus
  
  # mean sample size of the experiment
  mean_n <- (n_Bminus + n_Bminus_Aplus + n_Bplus + n_Bplus_Aplus)/4

  
  
  
  
  # output
  data.frame(par_id = ID_column, RR_Bminus.Aplus, varRR_Bminus.Aplus, size_n_Bminus, RR_Bplus.Aplus, varRR_Bplus.Aplus,
             size_n_Bplus, RR_interaction, varRR_interaction, mean_n, corrected_RR_Bminus.Aplus, corrected_varRR_Bminus.Aplus,
             corrected_RR_Bplus.Aplus, corrected_varRR_Bplus.Aplus, corrected_RR_interaction, corrected_varRR_interaction)
  
}