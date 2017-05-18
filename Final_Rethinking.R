library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(forcats)
library(lsmeans)
library(car)
library(broom)
library(lme4)
library(RLRsim)
library(sjPlot)
library(nlme)
library(knitr)
library(rethinking)

lipids <- read.csv("lipid_mod_bayes.csv") %>%
  na.omit(lipids) %>%
  mutate(Date_Collected = fct_inorder(Date_Collected)) %>%
  mutate(colony_idx = as.numeric(Colony)) %>%
  mutate(time_idx = as.numeric(Date_Collected)) %>%
  mutate(area_idx = as.numeric(Area)) %>%
  mutate(bleaching_idx = as.numeric(bleaching_level)) %>%
  mutate(upwelling_idx = as.numeric(Upwelling)-1)


lipids$log_lipids <- log(lipids$Total_Energetic_Lipid)

str(lipids)

base_plot <- ggplot(data =lipids, mapping = aes(x = Date_Collected, 
                                                y = log_lipids)) +
  geom_point() +
  theme_bw(base_size=17) +
  labs( x = "Date Collected", 
        y = "log(Total lipid (μg/mg))")

base_plot_2 <- ggplot(data =lipids, mapping = aes(x = Upwelling, 
                                                  y = log_lipids)) +
  geom_point() +
  theme_bw(base_size=17) +
  labs( x = "Date Collected", 
        y = "log(Total lipid (μg/mg))")


lipid_mod_13<- alist(
  #likelihood
  log_lipids ~ dnorm(mu, sigma),
  
  #DGP
  mu <- a_bar + a_area[area_idx] + a_time[time_idx] + a_colony[colony_idx] + (b_bar_upwelling + b_upwelling[upwelling_idx] + b_time[time_idx]) + (b_bar_colony + b_colony[colony_idx] + b_time[time_idx])*bleaching_level,
  
  #Random effects
  a_area[area_idx] ~ dnorm(0, sigma_area),
  c(a_colony,b_colony)[colony_idx] ~ dmvnormNC(sigma_colony, Rho_colony),
  
  #priors
  a_bar ~ dnorm(0,10),
  b_bar_upwelling ~ dnorm(0,10),
  b_bar_colony ~ dnorm(0,10),
  a_time[time_idx] ~ dnorm(0, 10),
  b_upwelling[upwelling_idx] ~ dnorm(0,10),
  b_time[time_idx] ~ dnorm(0,10),
  
  sigma ~ dcauchy(0,2),
  sigma_area ~ dcauchy(0,2),
  sigma_time ~ dcauchy(0,2),
  sigma_colony ~ dcauchy(0,2),
  Rho_colony ~ dlkjcorr(2)
)


fit_mod_13 <- map2stan(lipid_mod_13, data=lipids, iter=4000, chains=3, warmup=2000)

lipid_mod_14<- alist(
  #likelihood
  log_lipids ~ dnorm(mu, sigma),
  
  #DGP
  mu <- a_bar + a_area[area_idx] + a_colony[colony_idx] + (b_bar_upwelling + b_upwelling[upwelling_idx] + b_time[time_idx]) + (b_bar_colony + b_colony[colony_idx] + b_time[time_idx])*bleaching_level,
  
  #Random effects
  a_area[area_idx] ~ dnorm(0, sigma_area),
  c(a_colony,b_colony)[colony_idx] ~ dmvnormNC(sigma_colony, Rho_colony),
  
  #priors
  a_bar ~ dnorm(0,10),
  b_bar_upwelling ~ dnorm(0,10),
  b_bar_colony ~ dnorm(0,10),
  b_upwelling[upwelling_idx] ~ dnorm(0,10),
  b_time[time_idx] ~ dnorm(0,10),
  
  sigma ~ dcauchy(0,2),
  sigma_area ~ dcauchy(0,2),
  sigma_time ~ dcauchy(0,2),
  sigma_colony ~ dcauchy(0,2),
  Rho_colony ~ dlkjcorr(2)
)

fit_mod_14 <- map2stan(lipid_mod_14, data=lipids, iter=4000, chains=3, warmup=2000)

lipid_mod_15<- alist(
  #likelihood
  log_lipids ~ dnorm(mu, sigma),
  
  #DGP
  mu <- a_bar + a_colony[colony_idx] + (b_bar_upwelling + b_upwelling[upwelling_idx] + b_time[time_idx]) + (b_bar_colony + b_colony[colony_idx] + b_time[time_idx])*bleaching_level,
  
  #Random effects
  c(a_colony,b_colony)[colony_idx] ~ dmvnormNC(sigma_colony, Rho_colony),
  
  #priors
  a_bar ~ dnorm(0,10),
  b_bar_upwelling ~ dnorm(0,10),
  b_bar_colony ~ dnorm(0,10),
  b_upwelling[upwelling_idx] ~ dnorm(0,10),
  b_time[time_idx] ~ dnorm(0,10),
  
  sigma ~ dcauchy(0,2),
  sigma_area ~ dcauchy(0,2),
  sigma_time ~ dcauchy(0,2),
  sigma_colony ~ dcauchy(0,2),
  Rho_colony ~ dlkjcorr(2)
)

fit_mod_15 <- map2stan(lipid_mod_15, data=lipids, iter=4000, chains=3, warmup=2000)

lipid_mod_16<- alist(
  #likelihood
  log_lipids ~ dnorm(mu, sigma),
  
  #DGP
  mu <- a_bar + a_time[time_idx] + a_colony[colony_idx] + (b_bar_upwelling + b_upwelling[upwelling_idx] + b_time[time_idx]) + (b_bar_colony + b_colony[colony_idx] + b_time[time_idx])*bleaching_level,
  
  #Random effects
  c(a_colony,b_colony)[colony_idx] ~ dmvnormNC(sigma_colony, Rho_colony),
  
  #priors
  a_bar ~ dnorm(0,10),
  b_bar_upwelling ~ dnorm(0,10),
  b_bar_colony ~ dnorm(0,10),
  a_time[time_idx] ~ dnorm(0, 10),
  b_upwelling[upwelling_idx] ~ dnorm(0,10),
  b_time[time_idx] ~ dnorm(0,10),
  
  sigma ~ dcauchy(0,2),
  sigma_area ~ dcauchy(0,2),
  sigma_time ~ dcauchy(0,2),
  sigma_colony ~ dcauchy(0,2),
  Rho_colony ~ dlkjcorr(2)
)


fit_mod_16 <- map2stan(lipid_mod_16, data=lipids, iter=4000, chains=3, warmup=2000)

compare(fit_mod_13, fit_mod_14, fit_mod_15, fit_mod_16)
plot(compare(fit_mod_13, fit_mod_14, fit_mod_15, fit_mod_16))

post <- extract.samples(fit_mod_14)
str(post)

time_df <- as.data.frame(post$b_time)
names(time_df) <- levels(lipids$Date_Collected)

time_df_est <- time_df %>%
  gather(Date_Collected, log_lipids) %>%
  group_by(Date_Collected) %>%
  mutate(Pot="z")

base_plot +
  geom_point() +
  geom_point(data = time_df_est)


upwelling_df <- as.data.frame(post$b_upwelling)
names(upwelling_df) <- levels(lipids$Upwelling)

upwelling_df_est <- upwelling_df %>%
  gather(Upwelling, log_lipids) %>%
  group_by(Upwelling) %>%
  mutate(Pot="z")

base_plot_2 +
  geom_point() +
  geom_point(data = upwelling_df_est)

#conditional differnce

###Is May different than July?
May_v_July <- post$b_time[,1] - post$b_time[,2] 
plot(density(May_v_July))
#conditional distribution k
1-sum(May_v_July<0)/length(May_v_July)

ggplot(data=time_df) +
  geom_density(aes(x=time_df$`15-May`), alpha=.3) +
  geom_density(aes(x=time_df$`15-Jul`), alpha=.3)

###Is July different than 

July_v_March <- post$b_time[,2] - post$b_time[,3] 
plot(density(July_v_March))
#conditional distribution k
1-sum(July_v_March<0)/length(July_v_March)

#how to do this temporally vary by area? Upwelling vs No Upwelling
Upwelling_v_None <- (post$b_upwelling[,2] - post$b_upwelling[,1])
plot(density(Upwelling_v_None))

1-sum(Upwelling_v_None<0)/length(Upwelling_v_None)


#bayesian p value differences between time points by exaiming conditional differences 


