
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






