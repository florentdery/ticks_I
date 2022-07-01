
require(tidyverse)
ticks <- readr::read_csv("https://data.ecoforecast.org/targets/ticks/ticks-targets.csv.gz", guess_max = 1e6)
require(lubridate)

ticks=ticks %>% 
  rename(aa=amblyomma_americanum,
         ww=mmwrWeek) %>% 
  mutate(yy=year(time),
         mm=month(time))

head(ticks)
str(ticks)

ggplot(ticks)+
  geom_line(aes(ww, aa, color=factor(year(time))))+
  facet_wrap(~siteID)+
  theme_bw()

summary(ticks)
length(unique(ticks$siteID))
table(ticks$siteID)
table(ticks$yy)


#Ready to modify model :

ticks_sub=ticks %>% 
  subset(yy %in% c(min(yy):2019))

ggplot(ticks_sub)+
  geom_point(aes(y=aa, x=ww))

time_sub = as.Date(ticks_sub$time)
ww_sub=ticks_sub$ww
y_sub = ticks_sub$aa
# time_sub = as.Date(ticks_sub[ticks_sub$aa>0,]$time)
# ww_sub=ticks_sub[ticks_sub$aa>0,]$ww
# y_sub = ticks_sub[ticks_sub$aa>0,]$aa
# 
y_sub
data <- list(y=y_sub,n=length(y_sub),      ## data
             x_ic=log(1000),tau_ic=100, ## initial condition prior
             a_obs=1,r_obs=1,           ## obs error prior
             a_add=1,r_add=1            ## process error prior
)


require(R2jags)

RandomWalk = "
model{

  #### Data Model
  for(t in 1:n){
    y[t]~ dnorm(x[t],tau_obs)
  }
  
  #### Process Model
  for(t in 2:n){
    z[t]~ dnorm(x[t-1],tau_add)
    x[t]<-max(0,z[t])
  }
  
  #### Priors
  x[1] ~ dnorm(x_ic,tau_ic)
  tau_obs ~ dgamma(a_obs,r_obs)
  tau_add ~ dgamma(a_add,r_add)
}
"

nchain = 3
init <- list()
for(i in 1:nchain){
  y.samp = sample(y_sub,length(y_sub),replace=TRUE)
  init[[i]] <- list(tau_add=1/var(diff(log(y.samp))),  ## initial guess on process precision
                    tau_obs=5/var(log(y.samp)))        ## initial guess on obs precision
}

j.modelsub   <- jags.model (file = textConnection(RandomWalk),
                            data = data,
                            inits = init,
                            n.chains = 3)

jags.outsub   <- coda.samples (model = j.modelsub,
                               variable.names = c("tau_add","tau_obs"),
                               n.iter = 5000)


#zero_inflated model: 


zirw = "
model{

  #### Data Model
  for(t in 1:n){
    y[t]~ dlnorm(x[t],tau_obs)
  }
  
  #### Process Model
  for(t in 2:n){
    w[t]~ dbern(psi)
    z[t] <- w[t]*x[t-1]
    x[t] ~ dlnorm(z[t], tau_add)
  }
  
  #### Priors
  x[1] ~ dpois(x_ic)
  tau_obs ~ dgamma(a_obs,r_obs)
  tau_add ~ dgamma(a_add,r_add)
  psi ~ dunif(0, 1)
}
"


data <- list(y=y_sub,n=length(y_sub),      ## data
             x_ic=log(1000),tau_ic=100, ## initial condition prior, tau not useful since poisson, but kept it
             a_obs=1,r_obs=1,           ## obs error prior
             a_add=1,r_add=1            ## process error prior
)

nchain = 3
init <- list()
for(i in 1:nchain){
  y.samp = sample(y_sub,length(y_sub),replace=TRUE)
  init[[i]] <- list(tau_add=1/var(diff(log(y.samp+0.00001))),  ## initial guess on process precision
                    tau_obs=5/var(log(y.samp+0.00001)))        ## initial guess on obs precision
}

j.modelsub   <- jags.model (file = textConnection(zirw),
                            data = data,
                            inits = init,
                            n.chains = 3)

jags.outsub   <- coda.samples (model = j.modelsub,
                               variable.names = c("psi", "tau_add","tau_obs"),
                               n.iter = 5000)





time.rng = c(1,length(time_sub))       ## adjust to zoom in and out
out <- as.matrix(jags.outsub)         ## convert from coda to matrix  
x.cols <- c(grep("^tau",colnames(out)), grep("psi",colnames(out)) )## grab all columns that start with the letter x
ci <- apply(exp(out[,x.cols]),2,quantile,c(0.025,0.5,0.975)) ## model was fit on log scale

plot(ww_sub,
     ci[2,],
     type='n',
     ylim=range(y_sub,na.rm=TRUE),
     ylab="Flu Index",
     log='y_sub',xlim=ww_sub[time.rng])
## adjust x-axis label to be monthly if zoomed
if(diff(time.rng) < 100){ 
  axis.Date(1, at=seq(time[time.rng[1]],time[time.rng[2]],by='month'), format = "%Y-%m")
}
ecoforecastR::ciEnvelope(time,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(time,y,pch="+",cex=0.5)
plot(jags.outsub)




#Flo' edit : ####

ticks_standardized2 <- ticks %>% 
  
  select(siteID, time ,amblyomma_americanum ) %>%
  subset(siteID=="TALL") %>% 
  mutate(year = lubridate::year(time),
         month = lubridate::month(time)) %>%
  
  # add leading zero to month
  mutate(month = case_when(nchar(month) == 1 ~ paste0("0",month),
                           TRUE ~ as.character(month)),
         aa=case_when(year==2020~NA_real_,
                      TRUE~amblyomma_americanum)) %>%
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
                     values_from = aa)  %>%
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
ticks_standardized2 %>%
  mutate(date = as.Date(paste0(year_month,"-01"))) %>%
  tidyr::pivot_longer(cols = c("TALL"),
                      names_to = "siteID",
                      values_to = "tick_density") %>%
  tidyr::drop_na(tick_density) %>%
  ggplot(aes(x=date,
             y = tick_density,
             group = siteID,
             color = siteID)) +
  geom_point() +
  geom_path()


ticks_standardized_matrix2 <- ticks_standardized2 %>%
  select(-year_month, 
         -month,
         -year,
         -timestep) %>%
  # replace all zeros with NAs
  mutate_all(~replace(., . == 0, NA))


# ticks_standardized_matrix2 <- ticks_standardized2 %>%
#   select(-year_month, 
#          -month,
#          -year,
#          -timestep) %>%
#   # replace all zeros with NAs
#   mutate_all(~replace(., . == 0, NA))

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
data_groups_standardized2 <- list(y=as.matrix(ticks_standardized_matrix2),     ## data in a matrix 
                                  n_site=ncol(ticks_standardized_matrix2),     ## number of sites
                                  #n_site = 8,
                                  n = nrow(ticks_standardized_matrix2),        ## number of times (58)
                                  x_ic=100,tau_ic=.001,        ## initial condition prior
                                  a_obs=.005,r_obs=.005,           ## obs error prior
                                  a_add=.005,r_add=.005            ## process error prior
)



## d. build JAGS framework ----------------------------------------------------

nchain = 3



## e. run model ---------------------------------------------------------------


j.model_groups_standardized   <- jags.model (file = textConnection(RandomWalk_grouped_standardized),
                                             data = data_groups_standardized2,
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
  mutate(site_name = colnames(ticks_standardized_matrix)) %>%
  tidyr::unnest(data)

# plot group ribbons just to see
# tidy_out_groups_standardized_edited %>%
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
    ticks_standardized_matrix2 %>% 
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
  mutate(date = as.Date(paste0(ticks_standardized2$year_month,"-01")))%>%
  relocate(date) %>%
  tidyr::unnest(data)

# now plot with real values
tidy_out_groups_standardized_edited_merged %>%
  ggplot(aes(x=timestep, y = mean, group = site_num, fill = site_num)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = .15) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha=.5)+
  geom_line(size=.15) +
  # facet_grid(site_name~.,
  #            scales = "free_y") +
  ggthemes::theme_few()+
  geom_point(aes(x=timestep, y=tick_density), 
             shape = 21,
             size=1.5, alpha=.75,
             stroke = .4)+
  labs()



