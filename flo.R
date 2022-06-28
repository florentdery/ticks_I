
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

time_sub = as.Date(gflu_sub$Date)
y_sub = gflu_sub$Massachusetts

require(R2jags)

RandomWalk = "
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
  tau_add ~ dgamma(a_add,r_add)
}
"
data <- list(y=log(y_sub),n=length(y_sub),      ## data
             x_ic=log(1000),tau_ic=100, ## initial condition prior
             a_obs=1,r_obs=1,           ## obs error prior
             a_add=1,r_add=1            ## process error prior
)


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
                               n.iter = 1000)
plot(jags.outsub)

jags.outsub   <- coda.samples (model = j.modelsub,
                               variable.names = c("x","tau_add","tau_obs"),
                               n.iter = 10000)




