#Read in salmon count data
d1<-read.csv(file="~/Dropbox (Instream)/Projects/TWN Indian River/Data/IFR Cleaned Data/ir_salmon_esc_coupled_2013_16.csv",header=T, stringsAsFactor=F)

#calculate julian days, must add 365 to any days that overlap Dec. 31 
#(this is because some species are counted over the new year but are still techincally part of the previous years run)
jul<-strptime(d1$date, "%Y-%m-%d")
d1$jday<- ifelse(jul$yday>150,jul$yday,jul$yday+365)

#Order data set by year, species and day (so matches output in next step)
d2<-d1[order(d1$year, d1$species, d1$jday),]

#Create object with first survey date equalling 0 for each year/ species combination
day<-ddply(d2, c("year","species"), function(x){
  day<-x$jday-min(x$jday)
  data.frame(day)
})

#Cbind julian day into salmon count dataset
fish_dat<- cbind(d2,day$day)
names(fish_dat)[[9]]<- c("day")

#Format date
#fish_dat$date <- as.Date(fish_dat$date, format = "%Y-%m-%d")

#########################################################################################################################
#Apply the modeling from Millar et al 2012 CJFAS - Simple estimators of salmonid escapement and its variance using a new area-under-the-curve method.
#The code below applies a ML estimator using AUC from a Gaussian spawner model (GAUC).
#The GAUC model outperformed the Hilborn et al. 1999 and Trapizoidal models when used to calculate AUC escapement estimates from simulated and real data.
#########################################################################################################################
#Here, count and day are data vectors containing the values of ct and t, respectively. 
#sl and oe are the stream life (l) and observer efficiency (v), and sl.se and oe.se are their standard errors, respectively.

#Fit equation (8), and extract coefficients
#beta0, beta1 and beta2 
esc_dat<-ddply(fish_dat, c("year", "species"), .progress = progress_text(char = "."), function(y) {
  g=glm(count~day+I(day^2), data=y, family=quasipoisson); g
  x=coef(g)
  
  #Can use equation 8 to get predicted counts from model by applying the following equation y1=c[1]+(c[2]*days)+c[3]*(days^2); y<-exp(y1)
  #from equation 8, Millar et al. 2012
  #################################################################################################
  #These need to be agreed upon how to incorporate mean and uncertainty for sl and oe by IFR staff. 
  #By November 2014 we will have estimates of sl and oe for 3 years.
  #################################################################################################
  
  #Apply equation (9) to obtain estimated fish-days 
  f = sqrt(-pi / x[3]) * exp(x[1] - x[2] ^ 2 / (4 * x[3]))
  
  #Apply general OE and spp-specific SL estimates to each model
  oe<- mean_oe
  
  if(y$species[1]=="chum"){
    sl <- mean_chsl
    se_sl <- se_chsl
    e = f / (sl * oe)
    se_f=deltamethod(~sqrt(-pi/x3)*exp(x1-x2^2/(4*x3)), mean=x, cov=vcov(g))
    se_e=deltamethod(~x1/(x2*x3), mean=c(F,sl,oe), cov=diag(c(se_f,se_sl,se_oe))^2)
  }
  
  if(y$species[1]=="coho"){
    sl <- mean_cosl
    se_sl <- se_cosl
    e = f / (sl * oe)
    se_f=deltamethod(~sqrt(-pi/x3)*exp(x1-x2^2/(4*x3)), mean=x, cov=vcov(g))
    se_e=deltamethod(~x1/(x2*x3), mean=c(F,sl,oe), cov=diag(c(se_f,se_sl,se_oe))^2)
  }
  
  if(y$species[1]=="pink"){
    sl <- mean_pisl
    se_sl <- se_pisl
    e = f / (sl * oe)
    se_f=deltamethod(~sqrt(-pi/x3)*exp(x1-x2^2/(4*x3)), mean=x, cov=vcov(g))
    se_e=deltamethod(~x1/(x2*x3), mean=c(F,sl,oe), cov=diag(c(se_f,se_sl,se_oe))^2)
  }
  
  data.frame("oe" = round(oe, 3),
             "se_oe" = round(se_oe, 3),
             "sl" = round(sl, 1),
             "se_sl" = round(se_sl, 2),
             "escapement"=round(e, 0), 
             "escapement.se"=round(se_e, 0),
             "lower95CI"=round(e-(1.96*se_e), 0), 
             "upper95CI"=round(e+(1.96*se_e), 0))
})

esc_dat

write.csv(x=esc_dat, file="~/Dropbox (Instream)/Projects/TWN Indian River/4 - Figures & Tables/escapement_summary.csv")
