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




factorial_RR <- function(x_amb = dados$mean_nopred_Controle, sd_amb = dados$error_nopred_Controle, n_amb = dados$n_nopred_Controle, 
                            x_amb_pred = dados$mean_pred_Controle, sd_amb_pred = dados$error_pred_Controle, n_amb_pred = dados$n_pred_Controle,
                            x_clim = dados$mean_nopred_Tratamento, sd_clim = dados$error_nopred_Tratamento, n_clim = dados$n_nopred_Tratamento, 
                            x_clim_pred = dados$mean_pred_Tratamento, sd_clim_pred = dados$error_pred_Tratamento, n_clim_pred = dados$n_pred_Tratamento, 
                            ID_column = dados$par_id) { 
  
  
  
  # objects that will not go into the output of the function, but are middle steps to calculate other things
  # starts here
  # weighted variance for each factor of the factorial meta-analysis
  wghvar_treat.natural <- make_weight_var(x_amb_pred, sd_amb_pred, n_amb_pred)
  wghvar_cont.natural <- make_weight_var(x_amb, sd_amb, n_amb)
  wghvar_treat.altered <- make_weight_var(x_clim_pred, sd_clim_pred, n_clim_pred)
  wghvar_cont.altered <- make_weight_var(x_clim, sd_clim, n_clim)

  # bias correction for the individual weighted variance of the response ratio
  corrwghvar_treat.natural <- make_corrweight_var(x_amb_pred, sd_amb_pred, n_amb_pred)
  corrwghvar_cont.natural <- make_corrweight_var(x_amb, sd_amb, n_amb)
  corrwghvar_treat.altered <- make_corrweight_var(x_clim_pred, sd_clim_pred, n_clim_pred)
  corrwghvar_cont.altered <- make_corrweight_var(x_clim, sd_clim, n_clim)
  
  # ends here
  
  
  
  
  # calculation of the log response ratio and its variance
  # response ratio for the treatment under natural condition
  RR_treat.natural <- make_RR(x_amb_pred, x_amb)
  
  # response ratio for the treatment under altered condition
  RR_treat.altered <- make_RR(x_clim_pred, x_clim)
  
  # variance of the response ratio for the treatment under natural condition  
  varRR_treat.natural <-  wghvar_treat.natural + wghvar_cont.natural
  
  # variance of the response ratio for the treatment under altered condition
  varRR_treat.altered <- wghvar_treat.altered + wghvar_cont.altered  
  
  # sample size under natural condition
  size_n_natural <- (n_amb + n_amb_pred)/2
  
  # sample size under altered condition
  size_n_altered <- (n_clim + n_clim_pred)/2
  
  
  
  
  
  # calculation of the bias corrected log response ratio and its variance
  # bias corrected response ratio for the treatment under natural condition
  corrected_RR_treat.natural <- RR_treat.natural + 1/2*(wghvar_treat.natural - wghvar_cont.natural)
  
  # bias corrected response ratio for the treatment under altered condition
  corrected_RR_treat.altered <- RR_treat.altered + 1/2*(wghvar_treat.altered - wghvar_cont.altered)
  
  # bias corrected variance for the response ratio of the treatment under natural condition
  corrected_varRR_treat.natural <- varRR_treat.natural + 1/2*(corrwghvar_treat.natural + corrwghvar_cont.natural)
  
  # bias corrected variance for the response ratio of the treatment under altered condition
  corrected_varRR_treat.altered <- varRR_treat.altered + 1/2*(corrwghvar_treat.altered + corrwghvar_cont.altered)
  
  
  
  
  
  # calculation of the log response ratio for the interaction term between treatment and condition
  # log response ratio for the interaction
  RR_interaction <- RR_treat.altered - RR_treat.natural
  
  # bias correction of the log response ratio for the interaction
  corrected_RR_interaction <- corrected_RR_treat.altered - corrected_RR_treat.natural
  
  # variance of the log response ratio for the interaction
  varRR_interaction <- varRR_treat.natural + varRR_treat.altered
  
  # bias corrected variance of the log response ratio for the interaction
  corrected_varRR_interaction <-  corrected_varRR_treat.natural +  corrected_varRR_treat.altered
  
  # mean sample size of the experiment
  mean_n <- (n_amb + n_amb_pred + n_clim + n_clim_pred)/4

  
  
  
  
  # output
  data.frame(par_id = ID_column, RR_treat.natural, varRR_treat.natural, size_n_natural, RR_treat.altered, varRR_treat.altered,
             size_n_altered, RR_interaction, varRR_interaction, mean_n, corrected_RR_treat.natural, corrected_varRR_treat.natural,
             corrected_RR_treat.altered, corrected_varRR_treat.altered, corrected_RR_interaction, corrected_varRR_interaction)
  
}