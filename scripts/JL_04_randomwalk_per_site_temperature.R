
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# Forecast with a fixed effect  ----------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 

# in this attempt, I'll use the edit the per-site forecast to add a fixed effect
# of temperature to predict in each site. 

ticks <- readr::read_csv("https://data.ecoforecast.org/targets/ticks/ticks-targets.csv.gz", guess_max = 1e6)
site_data <- readr::read_csv("https://raw.githubusercontent.com/eco4cast/neon4cast-ticks/master/Ticks_NEON_Field_Site_Metadata_20210928.csv")

site_temps <- read.csv(here::here("data","daymetChallengeSites.csv")) %>%
  janitor::clean_names() %>%
  mutate(date = as.Date(date))
range(site_temps$date) 
# so this data stops at december 2020,
# so we'll have to use forecasted temps to project past that? 


# get monthly script for site temps
site_temps_monthly <- site_temps %>%
  dplyr::select(site_id,
                date,
                max_temperature,
                precipitation)%>%
  mutate(year = lubridate::year(date),
         month = lubridate::month(date)) %>%
  
  # add leading zero to month
  mutate(month = case_when(nchar(month) == 1 ~ paste0("0",month),
                           TRUE ~ as.character(month))) %>%
  
  # find monthly values
  group_by(site_id, month, year) %>%
  summarize(
    monthly_max_temp = max(max_temperature),
    mean_max_temp = mean(max_temperature),
    mean_precip = mean(precipitation),
    total_precip = sum(precipitation)
  ) %>%
  ungroup() %>%
  
  # format back to year-month column
  mutate(year_month = paste0(year,"-",month)) %>%
  select(-month, -year) %>%
  relocate(year_month) %>%
  
  # rename site_id to match other dataset
  rename(siteID = site_id)






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
  
  select(-n_samples) 



ticks_standardized_w_temps <- 
  ticks_standardized %>%
  left_join(site_temps_monthly)

# plot tick densities against counts, just to see

ticks_standardized_w_temps %>%
  mutate(amblyomma_americanum = case_when(amblyomma_americanum == 0 ~ as.numeric(NA),
                                          TRUE ~ amblyomma_americanum)) %>%
  ggplot(aes(x=mean_max_temp, y = amblyomma_americanum, color = siteID)) +
  geom_point()+
  ggthemes::theme_few()

ticks_standardized_w_temps %>%
  ggplot(aes(x=monthly_max_temp, y = amblyomma_americanum, color = siteID)) +
  geom_point()+
  ggthemes::theme_few()



ticks_standardized_w_temps %>%
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






# remove the last few years of points 
# (change to NA in all sites after timestep 85, 
# which is Jan 2021

ticks_standardized_trimmed <- ticks_standardized %>%
  mutate(across(BLAN:UKFS, 
                ~  case_when(timestep > 85 ~ as.double(NA),
                             TRUE ~ .)))



# plot the standardized dataset
ticks_standardized_trimmed %>%
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


ticks_standardized_matrix_trimmed <- ticks_standardized_trimmed %>%
  select(-year_month, 
         -month,
         -year,
         -timestep) %>%
  # replace all zeros with NAs
  mutate_all(~replace(., . == 0, NA))



## b. build JAGS model ------------------------------------------------------

RandomWalk_grouped_month_effect= "
model{
  
  #### Data Model
  for(t in 1:n){
    for (site in 1:n_site){
        y[t,site] ~ dnorm(x[t,site],tau_obs[site])
    }
  }
  
  #### Process Model
  for(t in 2:n){
    for(site in 1:n_site){
      x[t,site]~dnorm(x[t-1,site],tau_add[site])
    }
  }
  
  #### Priors
  for(site in 1:n_site){
    x[1,site] ~ dnorm(x_ic,tau_ic)
    
    tau_obs[site] ~ dgamma(a_obs,r_obs)
    tau_add[site] ~ dgamma(a_add,r_add)

  }
}
"


## c. make dataset ------------------------------------------------------------
data_groups_standardized_trimmed <- list(y=as.matrix(ticks_standardized_matrix_trimmed),     ## data in a matrix 
                                         n_site=ncol(ticks_standardized_matrix_trimmed),     ## number of sites
                                         #n_site = 2,
                                         n = nrow(ticks_standardized_matrix_trimmed),        ## number of times (58)
                                         x_ic=100,tau_ic=.001,        ## initial condition prior
                                         a_obs=.005,r_obs=.005,           ## obs error prior
                                         a_add=.005,r_add=.05            ## process error prior
)



## d. build JAGS framework ----------------------------------------------------

nchain = 3



## e. run model ---------------------------------------------------------------


j.model_groups_standardized_trimmed   <- jags.model (file = textConnection(RandomWalk_grouped_standardized),
                                                     data = data_groups_standardized_trimmed,
                                                     n.chains = 3)



jags.out_groups_standardized_trimmed   <- coda.samples (model = j.model_groups_standardized_trimmed,
                                                        variable.names = c("tau_add","tau_obs"),
                                                        n.iter = 10000)
plot(jags.out_groups_standardized_trimmed)
# ok, these don't look so good
gelman.diag(jags.out_groups_standardized_trimmed)



# run again extracting x variables
jags.out_groups_standardized_trimmed   <- coda.samples (model = j.model_groups_standardized_trimmed,
                                                        variable.names = c("tau_add","tau_obs","x"),
                                                        n.iter = 30000)



## f. tidy and visualize ------------------------------------------------------

tidy_out_groups_standardized_trimmed <- tidy(jags.out_groups_standardized_trimmed)

head(tidy_out_groups_standardized_trimmed)
# ok, this came out really inconvenient for tidyverse, where the 
# variables came out in one column as x[1,1] meaning time1, site1,
# so I'm going to use a bunch of string operations with stringr
# to take the value before the comma as time, and the value 
# after comma as site


tidy_out_groups_standardized_edited_trimmed <- tidy_out_groups_standardized_trimmed %>% 
  filter(stringr::str_detect(parameter,"x")) %>%
  mutate(timestep = stringr::str_replace(parameter, "x\\[",""),
         timestep = stringr::str_replace(timestep,",[0-9]\\]",""),
         site_num = stringr::str_replace(parameter, "x\\[[0-9],",""),
         site_num = stringr::str_replace(site_num, "x\\[[0-9][0-9],",""),
         site_num = stringr::str_replace(site_num, "\\]","")) %>%
  mutate(timestep = as.numeric(timestep)) %>%
  as_tibble() 

# add back group names
tidy_out_groups_standardized_edited_trimmed <- tidy_out_groups_standardized_edited_trimmed %>%
  group_by(site_num) %>%
  tidyr::nest() %>%
  ungroup() %>%
  mutate(site_name = colnames(ticks_standardized_matrix_trimmed)) %>%
  tidyr::unnest(data)

# plot group ribbons just to see
tidy_out_groups_standardized_edited_trimmed %>%
  ggplot(aes(x=timestep, y = mean, group = site_num, fill = site_num)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha=.5)+
  geom_line(size=.15) +
  facet_grid(site_num~.) +
  ggthemes::theme_few()



# add in real data from the untrimmed dataset
ticks_standardized_matrix_full <- ticks_standardized %>%
  select(-year_month, 
         -month,
         -year,
         -timestep) %>%
  # replace all zeros with NAs
  mutate_all(~replace(., . == 0, NA))


tidy_out_groups_standardized_edited_merged_trimmed <- 
  tidy_out_groups_standardized_edited_trimmed  %>%
  left_join(
    # join back to original dataset with real values
    ticks_standardized_matrix_full %>% 
      mutate(timestep = row_number())%>%
      tidyr::pivot_longer(cols = -timestep,
                          names_to = "site_name",
                          values_to = "tick_density") %>%
      mutate(forecast = case_when(timestep >= 85 ~ "yes",
                                  TRUE ~ "no"))
  )


# add in original dates
tidy_out_groups_standardized_edited_merged_trimmed <- tidy_out_groups_standardized_edited_merged_trimmed %>%
  group_by(timestep) %>%
  tidyr::nest() %>%
  ungroup() %>%
  mutate(date = as.Date(paste0(ticks_standardized_trimmed$year_month,"-01")))%>%
  relocate(date) %>%
  tidyr::unnest(data)



## g. plot and save ----------------------------------------------------

# now plot with real values
forecast_vis <- tidy_out_groups_standardized_edited_merged_trimmed %>% 
  ggplot(aes(x=date, y = mean, group = site_num, fill = site_num)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = .15) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha=.5)+
  geom_line(size=.15) +
  facet_grid(site_name~.,
             scales = "free_y") +
  ggthemes::theme_few()+
  # add in non-forecasted points (which the model was fit on)
  geom_point(data = . %>% filter(forecast == "no"),
             aes(x=date, y=tick_density), 
             shape = 21,
             size=1.5, alpha=.75,
             stroke = .4)+
  geom_vline(xintercept = as.Date("2021-01-01"),
             linetype = "dashed") +
  # add in forecasted points (after Jan 2021)
  geom_point(data = . %>% filter(forecast == "yes"),
             aes(x=date, y=tick_density), 
             shape = 8,
             size=1.5, alpha=.5,
             stroke = .4)+
  # add a shaded box to show forecasts
  geom_rect(inherit.aes=F,
            aes(xmin = as.Date("2021-01-01"),
                xmax = as.Date(Inf),
                ymin = -Inf,
                ymax = Inf),
            stat = "unique",
            fill = "grey10",
            alpha = .1) +
  labs(y = "Model Prediction",
       x = "Date",
       title ="Tick density predictions across 9 sites",
       subtitle = "Random walk models fitting tick densities across time, forecasting 2021")+
  theme(legend.position = "none",
        plot.title.position = "plot")


ggsave(
  here::here("plots","grouped_forecast_plot.png"),
  forecast_vis,
  dpi = 300,
  height = 8,
  width = 6,
  units = "in"
)



