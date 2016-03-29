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
oe<- read.csv(file="~/Dropbox (Instream)/Projects/TWN Indian River/Data/IFR Cleaned Data/Observer_Efficiency.csv",header=T, stringsAsFactor=T)
sl<- read.csv(file="~/Dropbox (Instream)/Projects/TWN Indian River/Data/IFR Cleaned Data/Survey_Life.csv",header=T, stringsAsFactor=T)

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

#Read in salmon count data
d1<-read.csv(file="~/Dropbox (Instream)/Projects/TWN Indian River/Data/IFR Cleaned Data/ir_salmon_esc_coupled_2013_15.csv",header=T, stringsAsFactor=T)

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
fish_dat$date <- as.Date(fish_dat$date, format = "%Y-%m-%d")

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
g=glm(count~day+I(day^2), data=y, family=quasipoisson); g #subset(z, year==2005)
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

#################
#GOOD TO HERE!!!!
#################

#dat <- dplyr::filter(fish_dat, species == "chum")
#sp<- "chum"

plot_fit <- function(dat, sp){

  xx<- dplyr::filter(dat, species == sp)
  
  fig_name<-paste("~/Dropbox (Instream)/Projects/TWN Indian River/4 - Figures & Tables/ir_auc_model_fit_",sp, ".png")
  png(fig_name, height=1200, width=1200)
  par(mfrow=c(2,1), mar=c(1.5,1.5,1,1.5), oma=c(4,4,0,0), cex=1.5)
  
  vv<- dplyr::arrange(xx,year,day)
  d_ply(vv, c("year"), function(x){
    g <- glm(count~day+I(day^2), data=x, family=quasipoisson); g; c<-coef(g)
	  days<-seq(from=range(x$day)[1], to=range(x$day)[2])
	  y1=c[1]+(c[2]*days)+c[3]*(days^2)
	  y<-exp(y1) #from equation 8, Millar et al. 2012
	  y_limit<-as.numeric(ifelse(max(x$count)>=max(y), max(x$count), max(y)))
	  
	  plot(x$count~x$jday, axes=FALSE, typ="b", xlab="", ylab="", pch=1, cex=1.5, lwd=2, col="blue",
	     ylim=c(0, y_limit*1.2), xlim=c(min(x$jday), max(x$jday)))

	  lines(y=y, x=days+min(x$jday), col="#00000050", lwd=2.5)

	  r <- range(x$date)
	  #r <- as.Date(range(auc$date))
	
	  #axis.POSIXct(1, at = seq(r[1], r[2], by = "weeks"), format = "%b %d", cex.axis=1.25)
	  #axis.Date(1, at = seq(r[1], r[2], by = "weeks"), cex.axis=1.25)
	
	  #axis(1, at = seq(min(x$jday), max(x$jday), by = 7), labels = seq(r[1], r[2], by = "weeks", format = "%b %d"), cex.axis=1.25)
	  axis(1, x$date, format(x$date, "%b %d"), cex.axis = 0.7)
	
	  mtext(x$year[1], side=3, adj=0.02, line=-1, cex=1.5)
	  axis(side=2, las=1, cex.axis=1.25)
	  box()
	  })

	mtext(paste(sp, "spawner count"), side=2, line=2, outer=TRUE, cex=3)
	mtext("Day of year", side=1, line=2, outer=TRUE, cex=3)#this provides the y-axis title.
	
  dev.off()#this closes the "device" and 
}

plot_fit(fish_dat, "chum")
plot_fit(fish_dat, "coho")
plot_fit(fish_dat, "pink")




fig.name<-"/Users/doug/Documents/Instream/Projects/BRGMON-3 - Bridge River/Figures/Bridge Chinook - Abundance by Year - Zeros Added.png"
png(fig.name, width=1200, height=1200)
par(mfrow=c(1,1), mar=c(1.5,1.5,1.5,1.5), oma=c(4,4,0,0), cex=1.5)
y_limit<-max(esc.data$escapement+esc.data$escapement.se, na.rm=TRUE)+1000
plot(escapement~year, data=subset(esc.data, year<1997), type="l", ylim=c(0, y_limit), xlim=c(min(esc.data$year), max(esc.data$year)), lwd=2, col=makeTransparent("black", 100), axes=FALSE); par(new=TRUE)
plot(escapement~year, data=subset(esc.data, year>1995), type="l", ylim=c(0, y_limit), xlim=c(min(esc.data$year), max(esc.data$year)), lwd=2, col=makeTransparent("black", 250), axes=FALSE); par(new=TRUE)
plot(escapement~year, data=subset(esc.data, year<1997), type="p", ylim=c(0, y_limit), xlim=c(min(esc.data$year), max(esc.data$year)), pch=16, cex=2, lwd=2, col=makeTransparent("black", 100), axes=FALSE); par(new=TRUE)
plot(escapement~year, data=subset(esc.data, year>1996), type="p", ylim=c(0, y_limit), xlim=c(min(esc.data$year), max(esc.data$year)), pch=16, cex=2, lwd=2, col=makeTransparent("black", 250), axes=FALSE)

#Add in se lines. Tried 95CI lines but they were to big for one of the estimates and didn't look good. Can change if need be.
for(i in 1:dim(esc.data)[1]){
lines(x=c(esc.data[i,1],esc.data[i,1]), y=c(esc.data[i,6]+esc.data[i,7], esc.data[i,6]-esc.data[i,7]), col=makeTransparent("black", 250), pch=19, lwd=2, lty=1)}
#Add the axes and box.
axis(1, at=seq(1994, 2014, 1), cex.axis=1.25)
axis(2, las=1, cex.axis=1.25)
box()
#Add the axes labels.
mtext("Chinook spawner abundance (Yalakom to Dam)", side=2, cex=3, outer=TRUE, line=2)
mtext("Year", side=1, cex=3, outer=TRUE, line=2)
#Add the legend.
mtext("Fence count", side=3, adj=0.02, line=-2, cex=2.75, col=makeTransparent("black", 100))#; points(x=1993,y=490, col=makeTransparent("blue", 100), pch=19, cex=1.25)
mtext("Visual survey", side=3, adj=0.02, line=-3.75, cex=2.75, col=makeTransparent("black", 250))#; points(x=340,y=460, col=makeTransparent("grey60", 100), lwd=2, pch=1, cex=1.25)

dev.off()









