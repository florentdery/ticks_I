
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
# project from 2020 to 2022  ----------------------------------
# --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 

# in this script, we'll use real tick data up to december 2020 to project
# tick densities up to december 2022 using a random walk model. 
# our assumption is that uncertainties will be very large in 2022,
# which later we will shrink with 2021 real data


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



# remove the last few years of points 
# (change to NA in all sites after timestep 85, 
# which is Jan 2021
ticks_standardized %>% glimpse()
ticks_standardized_trimmed <- ticks_standardized %>%
  # change every value past december 2020 to NA
  mutate(across(BLAN:UKFS, 
                ~  case_when(year_month > "2020-12" ~ as.double(NA),
                             TRUE ~ .))) 


# make a blank df for 2022 that we will bind to the dataframe above
blank_2022 <- 
  data.frame(
    year = rep(2022, 12),
    month = c(1:12)
  ) %>%
  mutate(year_month = paste0(year, "-",month)) %>%
  mutate("BLAN" = as.numeric(NA),
         "KONZ" = as.numeric(NA),
         "LENO" = as.numeric(NA),
         "ORNL" = as.numeric(NA),
         "OSBS" = as.numeric(NA),
         "SCBI" = as.numeric(NA),
         "SERC" = as.numeric(NA),
         "TALL" = as.numeric(NA),
         "UKFS" = as.numeric(NA),
         timestep = NA
         )



ticks_standardized_extended <- 
  ticks_standardized_trimmed %>%
  rbind(blank_2022) %>%
  mutate(timestep = 1:n())

rm(ticks_standardized_trimmed)

ticks_standardized_extended %>% glimpse()


# plot the standardized dataset
ticks_standardized_extended %>%
  mutate(date = as.Date(paste0(year_month,"-01"))) %>%
  tidyr::pivot_longer(cols = c("BLAN":"UKFS"),
                      names_to = "siteID",
                      values_to = "tick_density") %>%
  ggplot(aes(x=date,
             y = tick_density,
             group = siteID,
             color = siteID)) +
  geom_point() +
  geom_path()


ticks_standardized_matrix_extended <- ticks_standardized_extended %>%
  select(-year_month, 
         -month,
         -year,
         -timestep)# %>%
  # replace all zeros with NAs
 # mutate_all(~replace(., . == 0, NA))



## b. build JAGS model ------------------------------------------------------

RandomWalk_grouped_standardized= "
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
data_groups_standardized_extended <- list(y=as.matrix(ticks_standardized_matrix_extended),     ## data in a matrix 
                                         n_site=ncol(ticks_standardized_matrix_extended),     ## number of sites
                                         #n_site = 2,
                                         n = nrow(ticks_standardized_matrix_extended),        ## number of times (58)
                                         x_ic=100,tau_ic=.001,        ## initial condition prior
                                         a_obs=.005,r_obs=.005,           ## obs error prior
                                         a_add=.005,r_add=.05            ## process error prior
)



## d. build JAGS framework ----------------------------------------------------

nchain = 3



## e. run model ---------------------------------------------------------------


j.model_groups_standardized_extended   <- jags.model (file = textConnection(RandomWalk_grouped_standardized),
                                                     data = data_groups_standardized_extended,
                                                     n.chains = 3)



jags.out_groups_standardized_extended   <- coda.samples (model = j.model_groups_standardized_extended,
                                                        variable.names = c("tau_add","tau_obs"),
                                                        n.iter = 10000)
plot(jags.out_groups_standardized_extended)
# ok, these don't look so good
gelman.diag(jags.out_groups_standardized_extended)



# run again extracting x variables
jags.out_groups_standardized_extended   <- coda.samples (model = j.model_groups_standardized_extended,
                                                        variable.names = c("tau_add","tau_obs","x"),
                                                        n.iter = 30000)



## f. tidy and visualize ------------------------------------------------------

tidy_out_groups_standardized_extended <- tidy(jags.out_groups_standardized_extended)

head(tidy_out_groups_standardized_extended)
# ok, this came out really inconvenient for tidyverse, where the 
# variables came out in one column as x[1,1] meaning time1, site1,
# so I'm going to use a bunch of string operations with stringr
# to take the value before the comma as time, and the value 
# after comma as site


tidy_out_groups_standardized_edited_extended <- tidy_out_groups_standardized_extended %>% 
  filter(stringr::str_detect(parameter,"x")) %>%
  mutate(timestep = stringr::str_replace(parameter, "x\\[",""),
         timestep = stringr::str_replace(timestep,",[0-9]\\]",""),
         site_num = stringr::str_replace(parameter, "x\\[[0-9],",""),
         site_num = stringr::str_replace(site_num, "x\\[[0-9][0-9],",""),
         site_num = stringr::str_replace(site_num, "x\\[[0-9][0-9][0-9],",""),
         site_num = stringr::str_replace(site_num, "\\]","")) %>%
  mutate(timestep = as.numeric(timestep)) %>%
  as_tibble() 

# add back group names
tidy_out_groups_standardized_edited_extended <- tidy_out_groups_standardized_edited_extended %>%
  group_by(site_num) %>%
  tidyr::nest() %>%
  ungroup() %>%
  mutate(site_name = colnames(ticks_standardized_matrix_extended)) %>%
  tidyr::unnest(data)

# plot group ribbons just to see
tidy_out_groups_standardized_edited_extended %>%
  ggplot(aes(x=timestep, y = mean, group = site_num, fill = site_num)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha=.5)+
  geom_line(size=.15) +
  facet_grid(site_num~.) +
  ggthemes::theme_few()



# add in real data from the full dataset
ticks_standardized_matrix_full <- ticks_standardized %>%
  select(-year_month, 
         -month,
         -year,
         -timestep) %>%
  # replace all zeros with NAs
  mutate_all(~replace(., . == 0, NA))


tidy_out_groups_standardized_edited_merged_extended <- 
  tidy_out_groups_standardized_edited_extended  %>%
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
tidy_out_groups_standardized_edited_merged_extended <- 
  tidy_out_groups_standardized_edited_merged_extended %>%
  group_by(timestep) %>%
  tidyr::nest() %>%
  ungroup() %>%
  mutate(date = as.Date(paste0(ticks_standardized_extended$year_month,"-01")))%>%
  relocate(date) %>%
  tidyr::unnest(data)



## g. plot and save ----------------------------------------------------

# now plot with real values
forecast_vis <- tidy_out_groups_standardized_edited_merged_extended %>% 
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
       title ="Tick density 2-year forecast",
       subtitle = "Random walk models using values through Dec 2020 projecting to Dec 2022")+
  theme(legend.position = "none",
        plot.title.position = "plot",
        plot.title = element_text(face = "bold"))


forecast_vis

ggsave(
  here::here("plots","ticks_two_year_forecast.png"),
  forecast_vis,
  dpi = 300,
  height = 8,
  width = 6,
  units = "in"
)




