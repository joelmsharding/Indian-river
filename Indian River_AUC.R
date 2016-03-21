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

chsl<- filter(sl,spp=="chum")

pisl<- filter(sl,spp=="pink")

cosl<- filter(sl,spp=="coho")

#Function to generate data from specified mean and sd from OE data
oe_dist <- ddply(oe, c("stream"), function(x){ 
  dist <- x$mean+x$sd*scale(rnorm(200))
  data.frame(dist)
})

#Create mean and sd for OE spp-specific SL below
mean_oe<- mean(oe_dist$dist)
sd_oe<- sd(oe_dist$dist)

mean_chsl<- mean(chsl$x)
sd_chsl<- sd(chsl$x)

mean_cosl<- mean(cosl$x)
sd_cosl<- sd(cosl$x)

mean_pisl<- mean(pisl$x)
sd_pisl<- sd(pisl$x)


makeTransparent<-function(someColor, alpha=75)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

d1<-read.csv(file="~/Dropbox (Instream)/Projects/TWN Indian River/Data/IFR Cleaned Data/ir_salmon_esc_coupled_2013_15.csv",header=T, stringsAsFactor=T)

#Select species/year

ch13<- dplyr::filter(d1, species == "chum" & year == 2013)
co13<- dplyr::filter(d1, species == "coho" & year == 2013)
pi13<- dplyr::filter(d1, species == "pink" & year == 2013)
ch14<- dplyr::filter(d1, species == "chum" & year == 2014)
co14<- dplyr::filter(d1, species == "coho" & year == 2014)

#START HERE: NEED TO CREATE FUNCITON THAT AUTOMATES STRPTIME FOR EACH SPP/YEAR COMBO AND CALLS OE AND SPP-SPECIFIC SL ESTIMATES

#calculate julian days, must add 365 to any days in 2014 (otherwise will be assigned values prior to 2013 days)
strptime(sp$date, "%Y-%m-%d")
sp$jday<- ifelse(jday$yday<150,jday$yday,jday$yday+365) %>%
sp<-sp[order(sp$year, sp$jday),]

#Create object with first survey date equalling 0
day<-ddply(sp, c("year"), function(x){
	day<-x$jday-min(x$jday)
  data.frame(day)
})

#Cbind julian day into spp dataset
auc<- cbind(sp,day$day)
names(auc)[[9]]<- c("day")




#########################################################################################################################
#Apply the modeling from Millar et al 2012 CJFAS - Simple estimators of salmonid escapement and its variance using a new area-under-the-curve method.
#The code below applies a ML estimator using AUC from a Gaussian spawner model (GAUC).
#The GAUC model outperformed the Hilborn et al. 1999 and Trapizoidal models when used to calculate AUC escapement estimates from simulated and real data.
#########################################################################################################################
#Here, count and day are data vectors containing the values of ct and t, respectively. 
#sl and oe are the stream life (l) and observer efficiency (v), and sl.se and oe.se are their standard errors, respectively.

#Fit equation (8), and extract coefficients
#beta0, beta1 and beta2 
abund.est<-ddply(auc, c("year"), .progress = progress_text(char = "."), function(y) {
g=glm(count~day+I(day^2), data=y, family=quasipoisson); g #subset(z, year==2005)
x=coef(g)
#Can use equation 8 to get predicted counts from model by applying the following equation y1=c[1]+(c[2]*days)+c[3]*(days^2); y<-exp(y1)
#from equation 8, Millar et al. 2012
#################################################################################################
#These need to be agreed upon how to incorporate mean and uncertainty for sl and oe by IFR staff. 
#By November 2014 we will have estimates of sl and oe for 3 years.
#################################################################################################
# I used the mean observer efficiency for all years so that they are comparable. This might not be the best way to do it.
# Alternatively we could use the year specific observer efficiency. 
# Either method should be fine since it won't change the estimates all that much, they are still going to be low.
uncertainty <- data.frame("year" = c(2012, 2013, 2014), 
	"survey.life" = c(10, 11, 16), 
	"observer.efficiency" = c(0.58, 0.28, 0.28))
sl    <- mean(uncertainty$survey.life)
oe    <- mean(uncertainty$observer.efficiency)
sl.se <- sd(uncertainty$survey.life) / (sqrt(dim(uncertainty)[1]))
oe.se <- sd(uncertainty$observer.efficiency) / (sqrt(dim(uncertainty)[1]))

#Apply equation (9) to obtain estimated fish-days 
f = sqrt(-pi / x[3]) * exp(x[1] - x[2] ^ 2 / (4 * x[3]))

# I used year specific information on oe and sl when available. 
# The lowest observed oe was applied to 2014 because of poor water clarity.

if(y$year[1]==2012){
	oe <- 0.58
	sl <- 10
	e = f / (sl * oe)
		f.se=deltamethod(~sqrt(-pi/x3)*exp(x1-x2^2/(4*x3)), mean=x, cov=vcov(g))
	e.se=deltamethod(~x1/(x2*x3), mean=c(F,10,0.58), cov=diag(c(f.se,sl.se,oe.se))^2)
}

if(y$year[1]==2013){
	oe <- 0.28
	sl <- 11
	e = f / (sl * oe)
	f.se=deltamethod(~sqrt(-pi/x3)*exp(x1-x2^2/(4*x3)), mean=x, cov=vcov(g))
	e.se=deltamethod(~x1/(x2*x3), mean=c(F,11,0.28), cov=diag(c(f.se,sl.se,oe.se))^2)
}

if(y$year[1]==2014){
	oe <- 0.28
	sl <- 12
	e = f / (sl * oe)
	f.se=deltamethod(~sqrt(-pi/x3)*exp(x1-x2^2/(4*x3)), mean=x, cov=vcov(g))
	e.se=deltamethod(~x1/(x2*x3), mean=c(F,12,0.28), cov=diag(c(f.se,sl.se,oe.se))^2)
}

else
	e = f / (sl * oe)
	f.se=deltamethod(~sqrt(-pi/x3)*exp(x1-x2^2/(4*x3)), mean=x, cov=vcov(g))
	e.se=deltamethod(~x1/(x2*x3), mean=c(F,sl,oe), cov=diag(c(f.se,sl.se,oe.se))^2)

data.frame("oe" = round(oe, 3),
	"oe.se" = round(oe.se, 3),
	"sl" = round(sl, 1),
	"sl.se" = round(sl.se, 2),
	"escapement"=round(e, 0), 
	"escapement.se"=round(e.se, 0), 
	method="visual survey", 
	"lower95CI"=round(e-(1.96*e.se), 0), 
	"upper95CI"=round(e+(1.96*e.se), 0))
})

abund.est
#Add in the fence data. We assume no error with the fence numbers.
Zfc<-subset(Z2, fence.count==1)

fc<-data.frame("year"=Zfc$year,
	"oe" = NA,
	"oe.se" = NA,
	"sl" = NA,
	"sl.se" = NA,
	"escapement"=Zfc$adult.live, 
	"escapement.se"=0, 
	method="fence count", 
	"lower95CI"=round(Zfc$adult.live, 0), 
	"upper95CI"=round(Zfc$adult.live, 0))

esc.data<-rbind(fc, abund.est)
esc.data[is.na(esc.data)] <- NA
esc.data
write.csv(x=esc.data, file="/Users/doug/Documents/Instream/Projects/BRGMON-3 - Bridge River/Tables/Bridge River Chinook AUC Estimates - Variable oe and sl.csv")

min.jday<-min(d$jday)
max.jday<-max(d$jday) 

fig.name<-"/Users/doug/Documents/Instream/Projects/BRGMON-3 - Bridge River/Figures/Bridge Chinook - Count by Day with Model Fit - Zeros Added.png"
png(fig.name, height=1200, width=1200)
par(mfrow=c(5,3), mar=c(1.5,1.5,1,1.5), oma=c(4,4,0,0), cex=1.5)
d_ply(d, c("year"), function(x){
	g <- glm(count~day+I(day^2), data=x, family=quasipoisson); g; c<-coef(g)
	days<-seq(from=range(x$day)[1], to=range(x$day)[2])
	y1=c[1]+(c[2]*days)+c[3]*(days^2); y<-exp(y1) #from equation 8, Millar et al. 2012

	y_limit<-as.numeric(ifelse(max(x$count)>=max(y), max(x$count), max(y)))
	plot(x$count~x$jday, axes=FALSE, typ="b", xlab="", ylab="", pch=1, cex=2, lwd=2, col="blue", 
	ylim=c(0, y_limit*1.2), xlim=c(min.jday, max.jday))

	lines(y=y, x=days+min(x$jday), col=makeTransparent("black", 175), lwd=2.5)

	axis(1, cex.axis=1.25)
	mtext(x$year[1], side=3, adj=0.02, line=-1, cex=1.5)
	axis(side=2, las=1, cex.axis=1.25)
	box()
	})
	mtext("Chinook spawner count (Yalakom to Dam)", side=2, line=2, outer=TRUE, cex=3)
	mtext("Day of year", side=1, line=2, outer=TRUE, cex=3)#this provides the y-axis title. 
	#Note it is outside of the ddply() so that it is only plotted once.
dev.off()#this closes the "device" and 

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









