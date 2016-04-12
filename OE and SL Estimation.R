library(plyr)
library (dplyr)
library(nlme)
library(beanplot)
library(MuMIn)
library(influence.ME)
library(arm)
library(msm)
library(magrittr)

#Import SL and OE values from literature
oe<- read.csv(file="~/Dropbox (Instream)/Projects/TWN Indian River/Data/IFR Cleaned Data/Observer_Efficiency.csv",header=T, stringsAsFactor=F)
sl<- read.csv(file="~/Dropbox (Instream)/Projects/TWN Indian River/Data/IFR Cleaned Data/Survey_Life.csv",header=T, stringsAsFactor=F)

#detach(package:plyr)

#str(sl)

#sl1<- sl %>%
#  group_by(spp) %>%
#  dplyr::mutate(mean_sl=mean(x),sd_sl=sd(x),n_sl=count(x))

#Select SL values for each species
chsl<- filter(sl,spp=="chum")

pisl<- filter(sl,spp=="pink")

cosl<- filter(sl,spp=="coho")

#Function to generate data from specified mean and sd from OE data
oe_dist <- ddply(oe, c("stream"), function(x){ 
  dist <- x$mean+x$sd*scale(rnorm(5))
  data.frame(dist)
})

#Calculate mean, sd, N for OE & spp-specific SL below
mean_oe<- mean(oe_dist$dist)
sd_oe<- sd(oe_dist$dist)
n_oe<- nrow(oe)
se_oe<- sd_oe/sqrt(n_oe)

mean_chsl<- mean(chsl$x)
sd_chsl<- sd(chsl$x)
n_chsl<- nrow(chsl)
se_chsl<- sd_chsl/sqrt(n_chsl)

mean_cosl<- mean(cosl$x)
sd_cosl<- sd(cosl$x)
n_cosl<- nrow(cosl)
se_cosl<- sd_cosl/sqrt(n_cosl)

mean_pisl<- mean(pisl$x)
sd_pisl<- sd(pisl$x)
n_pisl<- nrow(pisl)
se_pisl<- sd_pisl/sqrt(n_pisl)
