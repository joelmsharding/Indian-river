##########################
#Escapement estimate plots
##########################

#Plot function
spp_plot<- function(dat, spp, est1, est2, est3){
  
  #If estimates are "NA" then make value off plot space
  if(is.null(est1)) {
    est1 <- -50000
  }
  
  if(is.null(est2)) {
    est2 <- -50000
  }
  
  if(is.null(est3)) {
    est3 <- -50000
  }
  xx<- filter(dat, species == spp)
  #Calculate number of plots based on number of years
  yrs<-length(unique(xx$year))
  #Make height of plot based on number of years
  wd<- ifelse(yrs==3,750,ifelse(yrs==2,500,1000))
  fig_name<-paste("~/Dropbox (Instream)/Projects/TWN Indian River/4 - Figures & Tables/ir_auc_pop_estimate_",spp,".png")
  png(fig_name, width=wd, height=800)
  par(mfrow=c(1,1), mar=c(1.5,1.5,1,1.5), oma=c(4,5.7,0,0), cex=1.5)
  y_limit<- max(xx$upper95CI)
  
  #Colour of plot dependant on species type, may have to add more species into ifelse statements
  sp_col<- ifelse(spp=="Chum","blue",ifelse(spp=="Coho","red",ifelse(spp=="Pink","black","green")))
  
  plot(escapement~year, data=xx, type="p", ylim=c(0, y_limit), xlim=c(min(xx$year)-0.5, max(xx$year)+0.5), pch=16, cex=2, lwd=2, col=sp_col, axes=FALSE); par(new=TRUE)

  #Add in se lines. Tried 95CI lines but they were to big for one of the estimates and didn't look good. Can change if need be.
  for(i in 1:dim(xx)[1]){
    lines(x=c(xx[i,1],xx[i,1]), y=c(xx[i,9], xx[i,10]), col=sp_col, pch=19, lwd=2, lty=1)}
  
  #Add TWN escapement estimates
  points(xx$year[1],est1, pch=15, col="#00000060", cex=1.5)
  points(xx$year[2],est2, pch=15, col="#00000060", cex=1.5)
  points(xx$year[3],est3, pch=15, col="#00000060", cex=1.5)
  
  
  #Add the axes and box.
  yy<- ifelse(spp=="Pink",2,1)
  axis(1, at=seq(min(xx$year), max(xx$year), by=yy), cex.axis=1.25)
  
  #Prevent y-axis from going into scientific notation
  myTicks = axTicks(2)
  axis(side=2, at = myTicks, labels = formatC(myTicks, format = 'd'), las=1, cex.axis=1.25)
  
  box()
  #Add the axes labels.
  mtext(paste(spp, "Population Estimate"), side=2, line=4, outer=TRUE, cex=2)
  mtext("Year", side=1, cex=2, outer=TRUE, line=2)
  
  #Add the legend.
  lx<- ifelse(spp=="Pink",yrs-2.5,yrs-1.75)
  legend(min(xx$year)+lx, y_limit-50, c("IFR AUC Estimate with 95 CI","TWN AUC Estimate"), pch=c(16,15), pt.cex=1, cex=0.8, col=c(sp_col,"#00000060"), bty="n", xpd=NA)

  
dev.off()
}

spp_plot(esc_dat,"Pink",740775,NA,NA)
spp_plot(esc_dat,"Coho",2454,4285,NA)
spp_plot(esc_dat,"Chum",49835,10623,NA)











