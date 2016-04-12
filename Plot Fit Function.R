##################################################################
#Function that creates plots based on model output and fish counts
##################################################################

plot_fit <- function(dat, sp){
  #Select the species to plot
  xx<- dplyr::filter(dat, species == sp)
  
  #Figure name will need to be changed based on desired file location
  fig_name<-paste("~/Dropbox (Instream)/Projects/TWN Indian River/4 - Figures & Tables/ir_auc_model_fit_",sp,".png")
  #Calculate number of plots based on number of years
  plot_no<-length(unique(xx$year))
  #Make height of plot based on number of years
  ht<- ifelse(plot_no==3,1800,ifelse(plot_no==2,1200,600))
  
  png(fig_name, height=ht, width=1500)
  
  #par is dependant on number of plots (plot_no) above
  par(mfrow=c(plot_no,1), mar=c(1.5,1.5,1,1.5), oma=c(4,5.2,0,0), cex=1.5)
  
  #Sort each data subet by day so plots properly
  vv<- dplyr::arrange(xx,year,day)
  
  d_ply(vv, c("year"), function(x){
    x$date <- as.POSIXct(as.character(x$date), format = "%Y-%m-%d")
    g <- glm(count~day+I(day^2), data=x, family=quasipoisson); g; c<-coef(g)
    days<-seq(from=range(x$day)[1], to=range(x$day)[2])
    #From equation 8, Millar et al. 2012
    y1=c[1]+(c[2]*days)+c[3]*(days^2)
    y<-exp(y1)
    y_limit<-as.numeric(ifelse(max(x$count)>=max(y), max(x$count), max(y)))
    
    #Colour of plot dependant on species type, may have to add more species into ifelse statements
    sp_col<- ifelse(sp=="chum","blue",ifelse(sp=="coho","red",ifelse(sp=="pink","black","green")))
    
    plot(x$count~as.Date(x$date), axes=FALSE, typ="b", xlab="", ylab="", pch=1, cex=2.5, lwd=3, col=sp_col,
         ylim=c(0, y_limit*1.2), xlim=c(min(as.Date(x$date)), max(as.Date(x$date))))
    
    x_range<- seq(min(as.Date(x$date)), max(as.Date(x$date)), by = "days")
    
    #Plot AUC model line
    lines(y=y, x=x_range, col="#00000050", lwd=2.5)
    
    #Create polygon to shade area under curve
    yzero<- rep(0, length(y))
    polygon(c(x_range, rev(x_range)), c(y, rev(yzero)), col = "#00000030", border = FALSE)
    
    r <- as.Date(range(xx$date))
    axis.Date(1, at = seq(r[1], r[2], by = "weeks"), format = "%b %d", cex.axis = 1.35)
    
    mtext(x$year[1], side=3, adj=0.02, line=-2, cex=2)
    myTicks = axTicks(2)
    axis(side=2, at = myTicks, labels = formatC(myTicks, format = 'd'), las=1, cex.axis=1.25)
    box()
  })
  
  mtext(paste(sp, "spawner count"), side=2, line=3.5, outer=TRUE, cex=2.5)
  mtext("Date", side=1, line=2, outer=TRUE, cex=2.5) #this provides the y-axis title.
  
  dev.off()
}

plot_fit(fish_dat, "chum")
plot_fit(fish_dat, "coho")
plot_fit(fish_dat, "pink")
