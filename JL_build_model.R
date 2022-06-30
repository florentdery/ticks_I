# script to build a random walk model based only on week



# libraries ---------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(rjags)
#remotes::install_github("njtierney/mmcc")
library(mmcc) # this package helps to tidy-verse-ify coda outputs


# download data -----------------------------------------------------------
ticks <- readr::read_csv("https://data.ecoforecast.org/targets/ticks/ticks-targets.csv.gz", guess_max = 1e6)
site_data <- readr::read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-ticks/master/Ticks_NEON_Field_Site_Metadata_20210928.csv")



# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# 1. make  model with no zeros -----------------------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

ticks_no_zero <- ticks %>%
  filter(amblyomma_americanum > 0)


RandomWalk= "
model{
  
  #### Data Model
  for(t in 1:n){
    y[t] ~ dnorm(x[t],tau_obs)
  }
  
  #### Process Model
  for(t in 2:n){
    x[t]~dnorm(x[t-1],tau_add)
  }
  
  #### Priors
  x[1] ~ dnorm(x_ic,tau_ic)
  tau_obs ~ dgamma(a_obs,r_obs)
  tau_add ~ dgamma(a_add,r_add) #  
}
"



## a. format data -------------------------------------------------------------

data <- list(y=ticks_no_zero$amblyomma_americanum,
             n=nrow(ticks_no_zero),     ## data
             x_ic=100,tau_ic=10,        ## initial condition prior
             a_obs=1,r_obs=1,           ## obs error prior
             a_add=1,r_add=1            ## process error prior
)



# X is predictied value!
# y is real value that we observed

## b. build JAGS framework ----------------------------------------------------

nchain = 3
init <- vector(mode ="list", length=3)




## c. run model ---------------------------------------------------------------
j.model   <- jags.model (file = textConnection(RandomWalk),
                         data = data,
                         n.chains = 3)



jags.out   <- coda.samples (model = j.model,
                            variable.names = c("tau_add","tau_obs"),
                            n.iter = 10000)
plot(jags.out)
# ok these plots don't look so good, but let's leave it for now
gelman.diag(jags.out)


# sample again including x variables
jags.out   <- coda.samples (model = j.model,
                            variable.names = c("tau_add","tau_obs","x"),
                            n.iter = 10000)



## d. tidy and plot data ------------------------------------------------------
tidy_out <- tidy(jags.out)

tidy_out_edited <- tidy_out %>%
  dplyr::filter(stringr::str_detect(parameter,"x")) %>% 
  mutate(timestep = stringr::str_replace(parameter,"x\\[",""),
         timestep = as.double(stringr::str_replace(timestep,"\\]",""))) %>%
  as_tibble() %>%
  mutate(time_from_dataset = ticks_no_zero  %>% pull(time),
         count_from_dataset = ticks_no_zero  %>% pull(amblyomma_americanum),
         site_from_dataset = ticks_no_zero  %>% pull(siteID))


tidy_out_edited %>%
  ggplot(aes(x=timestep, y = median)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = "lightblue") +
  geom_line() +
  theme_minimal() +
  geom_point(aes(x=timestep, 
                 y=count_from_dataset, 
                 color = site_from_dataset), size=.5)





# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# 2. repeat in a forecasting framework ---------------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---


## a. format data -------------------------------------------------------------


# let's take out some ticks from the dataset
ticks_no_zero_forecast <- ticks_no_zero %>%
  group_by(siteID) %>%
  mutate(amblyomma_americanum =
           case_when(
             row_number() >= (nrow(.)-5) ~ as.double(NA),
             TRUE ~ amblyomma_americanum
           )) %>%
  ungroup()

data_forecast <- list(y=ticks_no_zero_forecast$amblyomma_americanum,
             n=nrow(ticks_no_zero_forecast),     ## data
             x_ic=100,tau_ic=10,        ## initial condition prior
             a_obs=1,r_obs=1,           ## obs error prior
             a_add=1,r_add=1            ## process error prior
)


## b. build JAGS framework ----------------------------------------------------

nchain = 3
init <- vector(mode= "list", length = 3)



## c. run model ---------------------------------------------------------------


j.model_forecast   <- jags.model (file = textConnection(RandomWalk),
                         data = data_forecast,
                         n.chains = 3)



jags.out_forecast   <- coda.samples (model = j.model_forecast,
                            variable.names = c("tau_add","tau_obs"),
                            n.iter = 10000)
plot(jags.out_forecast)
# ok, these look good
gelman.diag(jags.out_forecast)


# resample model with x
jags.out_forecast   <- coda.samples (model = j.model_forecast,
                                     variable.names = c("tau_add","tau_obs","x"),
                                     n.iter = 10000)

## d. tidy and visualize ------------------------------------------------------

tidy_out_forecast <- tidy(jags.out_forecast)

tidy_out_edited_forecast <- tidy_out_forecast %>%
  dplyr::filter(stringr::str_detect(parameter,"x")) %>% 
  mutate(timestep = stringr::str_replace(parameter,"x\\[",""),
         timestep = as.double(stringr::str_replace(timestep,"\\]",""))) %>%
  as_tibble() %>%
  mutate(time_from_dataset = ticks_no_zero_forecast %>% pull(time),
         count_from_dataset = ticks_no_zero_forecast %>% pull(amblyomma_americanum),
         site_from_dataset = ticks_no_zero_forecast %>% pull(siteID))


tidy_out_edited_forecast %>%
  ggplot(aes(x=timestep, y = median)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = "lightblue") +
  geom_line() +
  theme_minimal() +
  geom_point(aes(x=timestep, y=count_from_dataset, color = site_from_dataset), size=.5)








# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# 3. try again grouped -------------------------------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
# here, I'm going to try to separate each site (of 9)
# 



## a. filter data -------------------------------------------------------------


ticks_no_zero <- ticks %>%
  filter(amblyomma_americanum > 0)


## b. build JAGS model ------------------------------------------------------

RandomWalk_grouped= "
model{
  
  #### Data Model
  for(t in 1:n){
    for (site in 1:n_site){
        y[t,site] ~ dnorm(x[t,site],tau_obs)

    }
  }
  
  #### Process Model
  for(t in 2:n){
    for(site in 1:n_site){
      x[t,site]~dnorm(x[t-1,site],tau_add)
    }
  }
  
  #### Priors
  for(site in 1:n_site){
    x[1,site] ~ dnorm(x_ic,tau_ic)
  }
  tau_obs ~ dgamma(a_obs,r_obs)
  tau_add ~ dgamma(a_add,r_add) #  
}
"


## c. make dataset ------------------------------------------------------------
ticks_no_zero <- ticks_no_zero %>%
  mutate(year = lubridate::year(time)) %>%
  mutate(year_week = paste0(year,"_",mmwrWeek)) %>%
  arrange(year_week)
 

# need to reformat to a matrix, so that JAGS can loop [1:time_n, 1:site_n]
# but, easiest way to do for now, is take all sites, arrange by date, 
# and assign timestep 1:n, meaning all will artificially start on timestep 1
# even if they are on different dates, and timesteps will still be assumed to 
# be evenly spaced. What will be better is if we can go by month/year, or something
# so we can *actually* regularize timesteps, and just keep NAs in where they belong
ticks_matrix <- ticks_no_zero %>% 
  select(siteID, year_week,amblyomma_americanum) %>% 
  group_by(siteID) %>%
  arrange(year_week) %>%
  mutate(timestep = 1:n()) %>%
  ungroup() %>% arrange(timestep) %>%
  select(-year_week) %>%
  tidyr::pivot_wider(names_from = siteID,
                     values_from = amblyomma_americanum) %>%
  select(-timestep) 



data_groups <- list(y=as.matrix(ticks_matrix),                ## data in a matrix 
             n_site=ncol(ticks_matrix),     ## number of sites
             n = nrow(ticks_matrix),        ## number of times (58)
             x_ic=100,tau_ic=10,        ## initial condition prior
             a_obs=1,r_obs=1,           ## obs error prior
             a_add=1,r_add=1            ## process error prior
)


## d. build JAGS framework ----------------------------------------------------

nchain = 3
init <- vector(mode= "list", length = 3)



## e. run model ---------------------------------------------------------------


j.model_groups   <- jags.model (file = textConnection(RandomWalk_grouped),
                                  data = data_groups,
                                  n.chains = 3)



jags.out_groups   <- coda.samples (model = j.model_groups,
                                     variable.names = c("tau_add","tau_obs"),
                                     n.iter = 10000)
plot(jags.out_groups)
# ok, these don't look so good
gelman.diag(jags.out_groups)

# run again, extracting x vars
jags.out_groups   <- coda.samples (model = j.model_groups,
                                   variable.names = c("tau_add","tau_obs","x"),
                                   n.iter = 10000)


## f. tidy and visualize ------------------------------------------------------

tidy_out_groups <- tidy(jags.out_groups)

tidy_out_groups
# ok, this came out really inconvenient for tidyverse, where the 
# variables came out in one column as x[1,1] meaning time1, site1,
# so I'm going to use a bunch of string operations with stringr
# to take the value before the comma as time, and the value 
# after comma as site


tidy_out_groups_edited <- tidy_out_groups %>% 
  filter(stringr::str_detect(parameter,"x")) %>%
  mutate(timestep = stringr::str_replace(parameter, "x\\[",""),
         timestep = stringr::str_replace(timestep,",[0-9]\\]",""),
         site_num = stringr::str_replace(parameter, "x\\[[0-9],",""),
         site_num = stringr::str_replace(site_num, "x\\[[0-9][0-9],",""),
         site_num = stringr::str_replace(site_num, "\\]","")) %>%
  mutate(timestep = as.numeric(timestep)) %>%
  as_tibble() 

# add back group names
tidy_out_groups_edited <- tidy_out_groups_edited %>%
  group_by(site_num) %>%
  tidyr::nest() %>%
  ungroup() %>%
  mutate(site_name = colnames(ticks_matrix)) %>%
  tidyr::unnest(data)

tidy_out_groups_edited %>%
  ggplot(aes(x=timestep, y = mean, group = site_num, fill = site_num)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha=.5)+
  geom_line(size=.15) +
  facet_grid(site_num~.) +
  ggthemes::theme_few()




# add in real data
tidy_out_groups_edited_merged <- 
  tidy_out_groups_edited  %>%
  left_join(
    # join back to original dataset with real values
    ticks_matrix %>% 
      mutate(timestep = row_number())%>%
      tidyr::pivot_longer(cols = -timestep,
                          names_to = "site_name",
                          values_to = "tick_density")
  )


# now plot with real values
tidy_out_groups_edited_merged %>%
  ggplot(aes(x=timestep, y = mean, group = site_num, fill = site_num)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha=.5)+
  geom_line(size=.15) +
  facet_grid(site_name~.) +
  ggthemes::theme_few()+
  geom_point(aes(x=timestep, y=tick_density), size=.1)+
  labs()



# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# 4. Fix matrix to even times ------------------------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 

ticks <- readr::read_csv("https://data.ecoforecast.org/targets/ticks/ticks-targets.csv.gz", guess_max = 1e6)
site_data <- readr::read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-ticks/master/Ticks_NEON_Field_Site_Metadata_20210928.csv")



## a. format data -------------------------------------------------------------


ticks_standardized <- ticks %>% 
  
  select(siteID, time ,amblyomma_americanum ) %>%
  
  mutate(year = lubridate::year(time),
         month = lubridate::month(time)) %>%
  
  # add leading zero to month
  mutate(month = case_when(nchar(month) == 1 ~ paste0("0",month),
                           TRUE ~ as.character(month))) %>%
  select(-time) %>%
  
  # since we are always missing data for jan/feb, add in those artificiatlly as NAs
  tidyr::complete(month = c("01","02","03","04",
                            "05","06","07","08",
                            "09","10","11","12"), year, siteID) %>%
  group_by(year, month,siteID) %>% 
  mutate(year_month = paste0(year,"-",month)) %>%
  group_by(siteID,year_month) %>%
  add_count(name = "n_samples") %>% #arrange(desc(n_samples)) %>%
  
  # when there are more than one sample in a site in a year/month, 
  # find the average density in the samples
  mutate(month  = as.double(month)) %>%
  summarise_all( .funs = "mean") %>%
  ungroup()  %>%
  
  select(-n_samples) %>%
  # pivot to a matrix of sites by year_month
  tidyr::pivot_wider(names_from = siteID,
                     values_from = amblyomma_americanum)  %>%
  arrange(year_month) %>%
  
  # add a numeric timestep based on year_month column
  group_by(year_month) %>%
  tidyr::nest() %>%
  ungroup() %>%
  mutate(timestep = c(1:n())) %>%
  tidyr::unnest(data) %>%
  ungroup() %>%
  relocate(timestep)


# plot the standardized dataset
ticks_standardized %>%
  mutate(date = as.Date(paste0(year_month,"-01"))) %>%
  tidyr::pivot_longer(cols = c("BLAN":"UKFS"),
                      names_to = "siteID",
                      values_to = "tick_density") %>%
  tidyr::drop_na(tick_density) %>%
  ggplot(aes(x=date,
             y = tick_density,
             group = siteID,
             color = siteID)) +
  geom_point() +
  geom_path()
  

ticks_standardized_matrix <- ticks_standardized %>%
  select(-year_month, 
         -month,
         -year,
         -timestep) %>%
  # replace all zeros with NAs
  mutate_all(~replace(., . == 0, NA))



## b. build JAGS model ------------------------------------------------------

RandomWalk_grouped_standardized= "
model{
  
  #### Data Model
  for(t in 1:n){
    for (site in 1:n_site){
        y[t,site] ~ dnorm(x[t,site],tau_obs)

    }
  }
  
  #### Process Model
  for(t in 2:n){
    for(site in 1:n_site){
      x[t,site]~dnorm(x[t-1,site],tau_add)
    }
  }
  
  #### Priors
  for(site in 1:n_site){
    x[1,site] ~ dnorm(x_ic,tau_ic)
  }
  tau_obs ~ dgamma(a_obs,r_obs)
  tau_add ~ dgamma(a_add,r_add) #  
}
"


## c. make dataset ------------------------------------------------------------
data_groups_standardized <- list(y=as.matrix(ticks_standardized_matrix),     ## data in a matrix 
                                 #n_site=ncol(ticks_standardized_matrix),     ## number of sites
                                 n_site = 1,
                                 n = nrow(ticks_standardized_matrix),        ## number of times (58)
                                 x_ic=100,tau_ic=.001,        ## initial condition prior
                                 a_obs=.005,r_obs=.005,           ## obs error prior
                                 a_add=.005,r_add=.005            ## process error prior
)



## d. build JAGS framework ----------------------------------------------------

nchain = 3



## e. run model ---------------------------------------------------------------


j.model_groups_standardized   <- jags.model (file = textConnection(RandomWalk_grouped_standardized),
                                data = data_groups_standardized,
                                n.chains = 3)



jags.out_groups_standardized   <- coda.samples (model = j.model_groups_standardized,
                                   variable.names = c("tau_add","tau_obs"),
                                   n.iter = 10000)
plot(jags.out_groups_standardized)
# ok, these don't look so good
gelman.diag(jags.out_groups_standardized)



# run again extracting x variables
jags.out_groups_standardized   <- coda.samples (model = j.model_groups_standardized,
                                                variable.names = c("tau_add","tau_obs","x"),
                                                n.iter = 10000)



## f. tidy and visualize ------------------------------------------------------

tidy_out_groups_standardized <- tidy(jags.out_groups_standardized)

tidy_out_groups_standardized
# ok, this came out really inconvenient for tidyverse, where the 
# variables came out in one column as x[1,1] meaning time1, site1,
# so I'm going to use a bunch of string operations with stringr
# to take the value before the comma as time, and the value 
# after comma as site


tidy_out_groups_standardized_edited <- tidy_out_groups_standardized %>% 
  filter(stringr::str_detect(parameter,"x")) %>%
  mutate(timestep = stringr::str_replace(parameter, "x\\[",""),
         timestep = stringr::str_replace(timestep,",[0-9]\\]",""),
         site_num = stringr::str_replace(parameter, "x\\[[0-9],",""),
         site_num = stringr::str_replace(site_num, "x\\[[0-9][0-9],",""),
         site_num = stringr::str_replace(site_num, "\\]","")) %>%
  mutate(timestep = as.numeric(timestep)) %>%
  as_tibble() 

# add back group names
tidy_out_groups_standardized_edited <- tidy_out_groups_standardized_edited %>%
  group_by(site_num) %>%
  tidyr::nest() %>%
  ungroup() %>%
  mutate(site_name = colnames(ticks_standardized_matrix[1])) %>%
  tidyr::unnest(data)

# plot group ribbons just to see
#tidy_out_groups_standardized_edited %>%
#  ggplot(aes(x=timestep, y = mean, group = site_num, fill = site_num)) +
#  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha=.5)+
#  geom_line(size=.15) +
#  facet_grid(site_num~.) +
#  ggthemes::theme_few()




# add in real data
tidy_out_groups_standardized_edited_merged <- 
  tidy_out_groups_standardized_edited  %>%
  left_join(
    # join back to original dataset with real values
    ticks_standardized_matrix %>% 
      mutate(timestep = row_number())%>%
      tidyr::pivot_longer(cols = -timestep,
                          names_to = "site_name",
                          values_to = "tick_density")
  )


# add in original dates
tidy_out_groups_standardized_edited_merged <- tidy_out_groups_standardized_edited_merged %>%
  group_by(timestep) %>%
  tidyr::nest() %>%
  ungroup() %>%
  mutate(date = as.Date(paste0(ticks_standardized$year_month,"-01")))%>%
  relocate(date) %>%
  tidyr::unnest(data)

# now plot with real values
tidy_out_groups_standardized_edited_merged %>%
  ggplot(aes(x=date, y = mean, group = site_num, fill = site_num)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = .15) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha=.5)+
  geom_line(size=.15) +
  facet_grid(site_name~.,
             scales = "free_y") +
  ggthemes::theme_few()+
  geom_point(aes(x=date, y=tick_density), 
             shape = 21,
             size=1.5, alpha=.75,
             stroke = .4)+
  labs()+
  theme(legend.position = "none")

