##########################
#Escapement estimate plots
##########################

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









