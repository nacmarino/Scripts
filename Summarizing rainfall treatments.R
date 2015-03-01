# what does the rainfall distribution mean? -------------------------------

#rainfall is the data frame I am working on
#it describes the rainfall pattern imposed in the experiment
#the columns in the data frame are
#tratamento = the rainfall treatment (categorical)
#dia = the the day of the experiment (numerical)
#precipitacao = how much water was added in each day (numerical)
head(rainfall)


#function to calculate the length of a sequence of zeroes in a vector
testvec <- rnbinom(30,size = 1, prob = 0.8)

#a function to count the largest number of zeros in a seq
#by Andrew MacDonald
n_max_zero <- function(vec){
  testvec_list <- rle(vec)
  where_zero <- which(testvec_list$values == 0)
  testvec_list$lengths[where_zero]
}
testvec
n_max_zero(testvec)

#pipe the function to get the total number and duration of chunks of zero rain
zerodays<-rainfall %>%
  filter(dia > 12) %>% #exclude the initial days, since all treatments share the initial period of no rain (we did it that way)
  group_by(tratamento) %>%
  do(nzero = n_max_zero(.$precipitacao)) %>%
  mutate(max_seq = max(nzero), zero_events = length(nzero), 
         sdzero = sd(nzero, na.rm=TRUE), meanzero = mean(nzero))

#given the control treatment, what is the distribution of rainfall?
rain_control <- rainfall %>%
  filter(precipitacao > 0, tratamento =="Controle") %>%
  summarise(quantile10 = quantile(precipitacao, 0.1),
            quantile25 = quantile(precipitacao, 0.25),
            quantile50 = quantile(precipitacao, 0.5),
            quantile75 = quantile(precipitacao, 0.75),
            quantile90 = quantile(precipitacao, 0.9))

#create a data frame only with the days that rained
rains <- rainfall %>%
  filter(precipitacao > 0) %>%
  group_by(tratamento) %>%
  summarise(small.event = sum(precipitacao <= rain_control$quantile10), 
            event.25 = sum(precipitacao > rain_control$quantile10 & precipitacao <= rain_control$quantile25),
            event.50 = sum(precipitacao > rain_control$quantile25 & precipitacao <= rain_control$quantile50),
            event.75 = sum(precipitacao > rain_control$quantile50 & precipitacao <= rain_control$quantile75),
            event.90 = sum(precipitacao > rain_control$quantile75 & precipitacao <= rain_control$quantile90),
            big.event = sum(precipitacao > rain_control$quantile90),
            total.rainfall = sum(precipitacao), rain_event = mean(precipitacao),
            max_event = max(precipitacao), min_event = min(precipitacao),
            q10 = mean(quantile(precipitacao, 0.1)), 
            q25 = mean(quantile(precipitacao, 0.25)),
            q50 = mean(quantile(precipitacao, 0.5)),
            q75 = mean(quantile(precipitacao, 0.75)), 
            q90 = mean(quantile(precipitacao, 0.9))) %>%
  left_join(zerodays)
head(rains)

#now, with your data in hand, summarise the characteristics of your rainfall
#distribution
raindist <- rainfall %>%
  group_by(tratamento) %>%
  summarise(no_rain = sum(precipitacao == 0), yes_rain = sum(precipitacao > 0)) %>%
  left_join(rains) %>%
  mutate(prop.small = round(small.event/yes_rain, digits = 2), 
         prop.25 = round(event.25/yes_rain, digits = 2), 
         prop.50 = round(event.50/yes_rain, digits = 2), 
         prop.75 = round(event.75/yes_rain, digits = 2), 
         prop.90 = round(event.90/yes_rain, digits = 2), 
         prop.big = round(big.event/yes_rain, digits = 2))
raindist

#no_rain = number of days with no rainfall
#yes_rain = number of days with rainfall
#small.event = number of small precipitation events (departing from control)
#event.25 = number of precipitation events between 10% and 25% quantile (departing from control)
#event.50 = number of precipitation events between 25% and 50% quantile (departing from control)
#event.75 = number of precipitation events between 50% and 75% quantile (departing from control)
#event.90 = number of precipitation events between 75% and 90% quantile (departing from control)
#big.event = number of precipitation events greater than 90% quantile (departing from control)
#total.rainfall = total volume of water added to the plant
#rain_event = mean volume of water added to the plant per event
#max_event = maximum volume of water added to the plant in a single event
#min_event = minimum volume of water added to the plant in a single event
#q10 = size of a small precipitation event for that treatment
#q25 = size of a precipitation even smaller than the 25% quantile for that treatment
#q50 = size of a precipitation even smaller than the 50% quantile for that treatment
#q75 = size of a precipitation even smaller than the 75% quantile for that treatment
#q90 = size of a precipitation even smaller than the 90% quantile for that treatment
#max_seq = maximum length of days where the plant received no water
#zero_events = number of times where the plant received no water
#sdzero = standard deviation of the number of days where the plant got no water
#meanzero = mean number of days where the plant got no water
#prop.x = proportion of events on the x size quantile