# Code for Anderson et al. (2021) 
# Marine Phytoplankton Functional Types Exhibit Diverse Responses to Thermal Change
# Stephanie I. Anderson updated 09/14/2020

# Contains:
## Phytoplankton functional type thermal reaction norms (Figure 1)
## Exponential curve fits
## Q10 calculations (Table 1)
## Comparison of thermal dependencies (Figure 2)
## Extended Figures 3 & 5

# Load packages
library(ggplot2)
library(quantreg)
library(lme4)
library(dplyr)
library(data.table)
library(cowplot)

# Load data
isolates <- read.csv("data/derived_traits.csv")
rates <- read.csv("data/growth_rates.csv")
source("scripts/nbcurve.R")
source("scripts/custom_theme.R")

#########################################################################
##### Thermal reaction norms #####
diatom<-subset(isolates, group=="diatoms")
cyano<-subset(isolates, group=="cyanobacteria")
dino<-subset(isolates, group=="dinoflagellates")
coccolithophores<-subset(isolates, group=="coccolithophores")
x<-seq(-2, 40, by=0.01) # temperature sequence

########## Diatoms ##########
# Fit 99th quantile regression
bissd <- rq(ln.r~temperature, data=subset(rates, group=="diatoms"),tau=0.99,ci=T)
cf.bd <- coef(bissd) #extract coefficients

# Calculate confidence intervals
## He & Hu (2002) method = "mcmb" uses the Markov chain marginal bootstrap
QR.b <- boot.rq(cbind(1,subset(rates, group=="diatoms")$temperature),
                subset(rates, group=="diatoms")$ln.r,tau=0.99, R=10000, method="mcmb")
ci_d <- t(apply(QR.b$B, 2, quantile, c(0.025,0.975)))

# plotting thermal performance curves
dev.off()
for(j in 1){
  pdf("figures/Diatom_TPC.pdf", width = 5.8, height = 4)
  plot.new()
  plot.window(c(-2,40),c(0,3))
  axis(1, 10*(-2:40), mgp=c(1,0.5,0))
  axis(2, 0.5*(0:6), mgp=c(1,0.5,0))
  box()
  for(i in 1:nrow(diatom)){
    o=diatom[i, "mu.c.opt.list"]
    w=diatom[i, "mu.wlist"]
    a=diatom[i, "mu.alist"]
    b=diatom[i, "mu.blist"]
    curve(nbcurve(x=x,opt=o,w=w,a=a,b=b),-2,40, col=alpha("black",alpha=0.6),ylim=c(0,3), 
        lty=1, add=T,xlab="", ylab="", cex=1.5)
  }
  # Add regression
  tempd <- seq(min(subset(rates, group=="diatoms")$temperature),
               max(subset(rates, group=="diatoms")$temperature), by=0.1)
  y1 <- c(exp(ci_d[1,2]+ci_d[2,2]*tempd))
  y2 <- c(exp(ci_d[1,1]+ci_d[2,1]*tempd))
  polygon(c(tempd, rev(tempd)),c(y1, rev(y2)),col=alpha("#026cb1",alpha=0.2), border=FALSE)
  curve(exp(cf.bd[[1]]+cf.bd[[2]]*x),min(subset(rates, group=="diatoms")$temperature),max(subset(rates, group=="diatoms")$temperature),add=T,col='#026cb1',lwd=2.5)
  
  # Eppley, 1972
  curve(0.59*exp(0.0633*x),-2,40,add=T,col='grey30',lwd=2.5, lty=2) # Eppley, 1972
  
  # add plot labels
  title(xlab=(expression(bold("Temperature (ºC)"))),
        ylab=(expression(bold("Specific Growth Rate (d"^"-1" *")"))), line=1.5, cex.lab=1)
  title(main=expression(bold("Diatoms")), line=-1, adj=0.05, cex=0.9)
  text(-1.6, 2.4, paste0("n=", length(diatom$isolate.code)), adj=c(0,0))
  text(-1.6, 2.1, paste0("N=", length(subset(rates, group=="diatoms")$isolate.code)), adj=c(0,0))
  dev.off()
}

########## Cyanobacteria ##########
# Fit 99th quantile regression
bissc_b<-rq(ln.r~temperature, data=subset(rates, group=="cyanobacteria"),tau=0.99,ci=T)
cf.bcb<-coef(bissc_b) #extract coefficients

# Calculate confidence intervals
QR.c <- boot.rq(cbind(1,subset(rates, group=="cyanobacteria")$temperature),
                subset(rates, group=="cyanobacteria")$ln.r,tau=0.99, R=10000, method="mcmb")
ci_cy<-t(apply(QR.c$B, 2, quantile, c(0.025,0.975)))

# plotting thermal performance curves
dev.off()
for(j in 1){
  pdf("figures/Cyanobacteria_TPC.pdf", width = 5.8, height = 4)
  plot.new()
  plot.window(c(-2,40),c(0,3))
  axis(1, 10*(-2:40), mgp=c(1,0.5,0))
  axis(2, 0.5*(0:6), mgp=c(1,0.5,0))
  box()
  for(i in 1:nrow(cyano)){
    o=cyano[i, "mu.c.opt.list"]
    w=cyano[i, "mu.wlist"]
    a=cyano[i, "mu.alist"]
    b=cyano[i, "mu.blist"]
    curve(nbcurve(x=x,opt=o,w=w,a=a,b=b),-2,40,ylim=c(0,3), col=alpha("black",alpha=0.6),lty=1, 
        xlab="", ylab="", add=T, cex=1.5)
  }
  # Add regression
  tempc <- seq(min(subset(rates, group=="cyanobacteria")$temperature),max(subset(rates, group=="cyanobacteria")$temperature), by=0.1)
  y1_cy <- c(exp(ci_cy[1,2]+ci_cy[2,2]*tempc))
  y2_cy <- c(exp(ci_cy[1,1]+ci_cy[2,1]*tempc))
  polygon(c(tempc, rev(tempc)),c(y1_cy, rev(y2_cy)),col=alpha("#ec3a25",alpha=0.2), border=FALSE)
  curve(exp(cf.bcb[[1]]+cf.bcb[[2]]*x),min(subset(rates, group=="cyanobacteria")$temperature),max(subset(rates, group=="cyanobacteria")$temperature),add=T,col='#ec3a25',lwd=2.5)
  
  # Eppley, 1972
  curve(0.59*exp(0.0633*x),-2,40,add=T,col='grey30',lwd=2.5, lty=2) # Eppley, 1972
  
  # add plot labels
  title(xlab=(expression(bold("Temperature (ºC)"))),
        ylab=(expression(bold("Specific Growth Rate (d"^"-1" *")"))), line=1.5, cex.lab=1)
  title(main=expression(bold("Cyanobacteria")), line=-1, adj=0.05, cex=0.9)
  text(-1.6, 2.4, paste0("n=", length(cyano$isolate.code)), adj=c(0,0))
  text(-1.6, 2.1, paste0("N=", length(subset(rates, group=="cyanobacteria")$isolate.code)),adj=c(0,0))
  dev.off()
}

########## Dinoflagellates ##########
# Fit 99th quantile regression
bissdi<-rq(ln.r~temperature, data=subset(rates, group=="dinoflagellates"),tau=0.99,ci=T) 
cf.di<-coef(bissdi) #extract coefficients

# Calculate confidence intervals
QR.d <- boot.rq(cbind(1,subset(rates, group=="dinoflagellates")$temperature),
                subset(rates, group=="dinoflagellates")$ln.r,tau=0.99, R=10000, method="mcmb")
ci_df<-t(apply(QR.d$B, 2, quantile, c(0.025,0.975)))

# plotting thermal performance curves
dev.off()
for(j in 1){
  pdf("figures/Dinoflagellate_TPC.pdf", width = 5.8, height = 4)
  plot.new()
  plot.window(c(-2,40),c(0,3))
  axis(1, 10*(-2:40), mgp=c(1,0.5,0))
  axis(2, 0.5*(0:6), mgp=c(1,0.5,0))
  box()
  for(i in 1:nrow(dino)){
    o=dino[i, "mu.c.opt.list"]
    w=dino[i, "mu.wlist"]
    a=dino[i, "mu.alist"]
    b=dino[i, "mu.blist"]
    curve(nbcurve(x=x,opt=o,w=w,a=a,b=b),-2,40,ylim=c(0,3),col=alpha("black",alpha=0.6), lty=1,
        xlab="", ylab="", add=T, cex=1.5)
  }
  # Add regression
  tempdi <- seq(min(subset(rates, group=="dinoflagellates")$temperature),
                max(subset(rates, group=="dinoflagellates")$temperature), by=0.1)
  y1_df <- c(exp(ci_df[1,2]+ci_df[2,2]*tempdi))
  y2_df <- c(exp(ci_df[1,1]+ci_df[2,1]*tempdi))
  polygon(c(tempdi, rev(tempdi)),c(y1_df, rev(y2_df)),col=alpha("#3ea127",alpha=0.2), border=FALSE)
  curve(exp(cf.di[[1]]+cf.di[[2]]*x),min(subset(rates, group=="dinoflagellates")$temperature),
        max(subset(rates, group=="dinoflagellates")$temperature),add=T,col='#3ea127',lwd=2.5) 
  
  # Eppley, 1972
  curve(0.59*exp(0.0633*x),-2,40,add=T,col='grey30',lwd=2.5, lty=2) # Eppley, 1972
  
  # add plot labels
  title(xlab=(expression(bold("Temperature (ºC)"))),
        ylab=(expression(bold("Specific Growth Rate (d"^"-1" *")"))), line=1.5, cex.lab=1)
  title(main=expression(bold("Dinoflagellates")), line=-1, adj=0.05,  cex=0.9)
  text(-1.6, 2.4, paste0("n=", length(dino$isolate.code)), adj=c(0,0))
  text(-1.6, 2.1, paste0("N=", length(subset(rates, group=="dinoflagellates")$isolate.code)), adj=c(0,0))
  dev.off()
}

########## Coccolithophores ##########
# Fit 99th quantile regression
bissco<-rq(ln.r~temperature, data=subset(rates, group=="coccolithophores"),tau=0.99,ci=T) #weights=wts 
cf.co<-coef(bissco) #extract coefficients

# Calculate confidence intervals
QR.c <- boot.rq(cbind(1,subset(rates, group=="coccolithophores")$temperature),
                subset(rates, group=="coccolithophores")$ln.r,tau=0.99, R=10000, method="mcmb")
ci_co<-t(apply(QR.c$B, 2, quantile, c(0.025,0.975)))

# plotting thermal performance curves
dev.off()
for(j in 1){
  pdf("figures/Coccolithophore_TPC.pdf", width = 5.8, height = 4)
  plot.new()
  plot.window(c(-2,40),c(0,3))
  axis(1, 10*(-2:40), mgp=c(1,0.5,0))
  axis(2, 0.5*(0:6), mgp=c(1,0.5,0))
  box()
  for(i in 1:nrow(coccolithophores)){
    o=coccolithophores[i, "mu.c.opt.list"]
    w=coccolithophores[i, "mu.wlist"]
    a=coccolithophores[i, "mu.alist"]
    b=coccolithophores[i, "mu.blist"]
    curve(nbcurve(x=x,opt=o,w=w,a=a,b=b),-2,40,ylim=c(0,3), col=alpha("black",alpha=0.6), lty=1, add=T,
        xlab="", ylab="", cex=1.5)
  }
  # Add regression
  tempco<-seq(min(subset(rates, group=="coccolithophores")$temperature),
              max(subset(rates, group=="coccolithophores")$temperature), by=0.1)
  y1_co <- c(exp(ci_co[1,2]+ci_co[2,2]*tempco))
  y2_co <- c(exp(ci_co[1,1]+ci_co[2,1]*tempco))
  polygon(c(tempco, rev(tempco)),c(y1_co, rev(y2_co)),col=alpha("orange",alpha=0.2), border=FALSE)
  curve(exp(cf.co[[1]]+cf.co[[2]]*x),min(subset(rates, group=="coccolithophores")$temperature),
        max(subset(rates, group=="coccolithophores")$temperature),add=T,col='orange',lwd=2.5)

  # Eppley, 1972
  curve(0.59*exp(0.0633*x),-2,40,add=T,col='grey30',lwd=2.5, lty=2) # Eppley, 1972

  # add plot labels
  title(xlab=(expression(bold("Temperature (ºC)"))),
        ylab=(expression(bold("Specific Growth Rate (d"^"-1" *")"))), line=1.5, cex.lab=1)
  title(main=expression(bold("Coccolithophores")), line=-1, adj=0.05, cex=0.9)
  text(-1.6, 2.4, paste0("n=", length(coccolithophores$isolate.code)), adj=c(0,0))
  text(-1.6, 2.1, paste0("N=", length(subset(rates, group=="coccolithophores")$isolate.code)), adj=c(0,0))
  dev.off()
}

########## All PFTs ##########
## (Extended Figure 3)
# Fit 99th quantile regression
biss<-rq(ln.r~temperature, data=rates, tau=0.99,ci=T)
cf.b<-coef(biss) #extract coefficients

# Calculate confidence intervals
QR.all <- boot.rq(cbind(1,rates$temperature),
                  rates$ln.r,tau=0.99, R=10000, method="mcmb")
ci<-t(apply(QR.all$B, 2, quantile, c(0.025,0.975)))

# plotting thermal performance curves
dev.off()
for(j in 1){
  pdf("figures/Extended_Figure3.pdf", width = 5.8, height = 4)
  plot.new()
  plot.window(c(-2,40),c(0,3))
  axis(1, 10*(-2:40), mgp=c(1,0.5,0))
  axis(2, 0.5*(0:6), mgp=c(1,0.5,0))
  box()
  for(i in 1:nrow(isolates)){
    o=isolates[i, "mu.c.opt.list"]
    w=isolates[i, "mu.wlist"]
    a=isolates[i, "mu.alist"]
    b=isolates[i, "mu.blist"]
    curve(nbcurve(x=x,opt=o,w=w,a=a,b=b),-2,40,ylim=c(0,3), col=alpha("black",alpha=0.6), lty=1, add=T,
          xlab="", ylab="", cex=1.5)
  }
  # Add regression
  temp<-seq(min(rates$temperature),
              max(rates$temperature), by=0.1)
  y1 <- c(exp(ci[1,2]+ci[2,2]*temp))
  y2 <- c(exp(ci[1,1]+ci[2,1]*temp))
  polygon(c(temp, rev(temp)),c(y1, rev(y2)),col=alpha("orangered1",alpha=0.2), border=FALSE)
  curve(exp(cf.b[[1]]+cf.b[[2]]*x),min(temp),max(temp),add=T,col='orangered1',lwd=2.5)
  
  # Eppley, 1972
  curve(0.59*exp(0.0633*x),-2,40,add=T,col='grey30',lwd=2.5, lty=2) # Eppley, 1972
  
  # add plot labels
  title(xlab=(expression(bold("Temperature (ºC)"))),
        ylab=(expression(bold("Specific Growth Rate (d"^"-1" *")"))), line=1.5)
  title(main=expression(bold("All PFTs")), line=-1, adj=0.05, cex=0.9)
  text(-1.6, 2.4, paste0("n=", length(isolates$isolate.code)), adj=c(0,0))
  text(-1.6, 2.1, paste0("N=", length(rates$isolate.code)), adj=c(0,0))
  dev.off()
}


########## Q10 Temperature Coefficient & Activation Energy ############
# Activation Energy
k = 8.617333262145*(10^(-5)) # boltzman constant
Ea<-function(b){ #slope (b)
  b*k*(273)^2
}

### Table 1 ###
group <- c('all.groups','coccolithophores','cyanobacteria', 'diatoms','dinoflagellates')
table1 <- as.data.frame(group)

# Number of unique isolates
table1$n <- rbind(length(isolates$isolate.code),length(coccolithophores$isolate.code),length(cyano$isolate.code),
          length(diatom$isolate.code),length(dino$isolate.code))

# Number of growth rate measurements
table1$N <- rbind(length(rates$isolate.code),length(subset(rates, group=="coccolithophores")$isolate.code),
               length(subset(rates, group=="cyanobacteria")$isolate.code),length(subset(rates, group=="diatoms")$isolate.code),
               length(subset(rates, group=="dinoflagellates")$isolate.code))

#Coefficients for exponential curves (calculated above)
table1$a <- rbind(round(cf.b[[1]],3), round(cf.co[[1]],3), round(cf.bcb[[1]],3), 
                                                           round(cf.bd[[1]],3), round(cf.di[[1]],3))

table1$a_ci <- rbind(paste0("[",round(ci[[1]],3),", ",round(ci[[3]],3),"]"),
                     paste0("[",round(ci_co[[1]],3),", ",signif(ci_co[[3]],3),"]"),
                     paste0("[",round(ci_cy[[1]],3),", ",round(ci_cy[[3]],3),"]"),
                     paste0("[",round(ci_d[[1]],3),", ",round(ci_d[[3]],3),"]"),
                     paste0("[",round(ci_df[[1]],3),", ",round(ci_df[[3]],3),"]"))

table1$b <- rbind(round(cf.b[[2]],3), round(cf.co[[2]],3), round(cf.bcb[[2]],3), 
                                                           round(cf.bd[[2]],3), round(cf.di[[2]],3))

table1$b_ci <- rbind(paste0("[",round(ci[[2]],3),", ",round(ci[[4]],3),"]"),
                     paste0("[",round(ci_co[[2]],3),", ",round(ci_co[[4]],3),"]"),
                     paste0("[",round(ci_cy[[2]],3),", ",round(ci_cy[[4]],3),"]"),
                     paste0("[",round(ci_d[[2]],3),", ",round(ci_d[[4]],3),"]"),
                     paste0("[",round(ci_df[[2]],3),", ",round(ci_df[[4]],3),"]"))

# calculated variables
table1$intercept = exp(table1$a) # y-intercept
table1$int_ci <- rbind(paste0("[",round(exp(ci[[1]]),3),", ",round(exp(ci[[3]]),3),"]"),
                     paste0("[",round(exp(ci_co[[1]]),3),", ",signif(exp(ci_co[[3]]),3),"]"),
                     paste0("[",round(exp(ci_cy[[1]]),3),", ",round(exp(ci_cy[[3]]),3),"]"),
                     paste0("[",round(exp(ci_d[[1]]),3),", ",round(exp(ci_d[[3]]),3),"]"),
                     paste0("[",round(exp(ci_df[[1]]),3),", ",round(exp(ci_df[[3]]),3),"]"))

table1$Q10 = exp(table1$b*10) # Q10
table1$Ea = Ea(table1$b) # acivation energy
table1$umax20 = exp(table1$a+table1$b*20) # maximum growth at 20ºC

write.csv(table1, "output/table1.csv")

###### Extended Data Figure 5 ######
cocco<-subset(rates, group =='coccolithophores')
cyano<-subset(rates, group =='cyanobacteria')
diatoms<-subset(rates, group =='diatoms')
dinos<-subset(rates, group =='dinoflagellates')

pdf("figures/Extended_Figure5.pdf", width = 7.2, height = 4.5)
x=tempco
par(mfrow=c(2,2), mar=c(0.5,3.5, 3.5, 0.5))
plot(cocco$temperature, cocco$r, xlim=c(-2, 40), ylim=c(0,3), xaxt='n',xlab='', ylab='', pch=20, col=alpha("black", 0.4))
axis(side=1,labels=F)
curve(exp(cf.co[[1]]+cf.co[[2]]*x),min(x), max(x),add=T,col='orange',lwd=2.5)
polygon(c(tempco, rev(tempco)),c(y1_co, rev(y2_co)),col=alpha("orange",alpha=0.2), border=FALSE)
curve(0.59*exp(0.0633*x),-2,40,add=T,col='grey30',lwd=2.5, lty=2)
title(main=expression(bold("Coccolithophores")), line=-1, adj=0.05, cex=1)
text(-1.6, 2.3, paste0("n=", length(unique(coccolithophores$isolate.code))), adj=c(0,0))
text(-1.6, 1.9, paste0("N=", length(subset(rates, group=="coccolithophores")$isolate.code)), adj=c(0,0))

x=tempc
par(mar=c(0.5,0.5,3.5,3.5))
plot(cyano$temperature, cyano$r, xlim=c(-2, 40), ylim=c(0,3), xaxt='n', yaxt='n',xlab='', ylab='', pch=20, col=alpha("black", 0.4))
axis(side=1,labels=F)
curve(exp(cf.bcb[[1]]+cf.bcb[[2]]*x),min(x), max(x),add=T,col=colors[2], lwd=2.5)
polygon(c(tempc, rev(tempc)),c(y1_cy, rev(y2_cy)),col=alpha(colors[2], alpha=0.2), border=FALSE)
curve(0.59*exp(0.0633*x),-2,40,add=T,col='grey30',lwd=2.5, lty=2)
title(main=expression(bold("Cyanobacteria")), line=-1, adj=0.05, cex=1)
text(-1.6, 2.3, paste0("n=", length(unique(cyano$isolate.code))), adj=c(0,0))
text(-1.6, 1.9, paste0("N=", length(subset(rates,group=="cyanobacteria")$isolate.code)), adj=c(0,0))

x=tempd
par(mar=c(3.5, 3.5,0.5,0.5))
plot(diatoms$temperature, diatoms$r, xlim=c(-2, 40), ylim=c(0,3),
     xlab='', ylab='', pch=20, col=alpha("black", 0.4))
curve(exp(cf.bd[[1]]+cf.bd[[2]]*x),min(x), max(x),add=T,col=colors[3],lwd=2.5)
polygon(c(tempd, rev(tempd)),c(y1, rev(y2)),col=alpha(colors[3],alpha=0.2), border=FALSE)
curve(0.59*exp(0.0633*x),-2,40,add=T,col='grey30',lwd=2.5, lty=2)
title(main=expression(bold("Diatoms")), line=-1, adj=0.05, cex=1)
text(-1.6, 2.3, paste0("n=", length(unique(diatoms$isolate.code))), adj=c(0,0))
text(-1.6, 1.9, paste0("N=", length(subset(rates, group=="diatoms")$isolate.code)), adj=c(0,0))

x=tempdi
par(mar=c(3.5, 0.5,0.5,3.5))
plot(dinos$temperature, dinos$r, xlim=c(-2, 40), ylim=c(0,3),
     xlab='', ylab='', yaxt='n',pch=20, col=alpha("black", 0.4))
curve(exp(cf.di[[1]]+cf.di[[2]]*x),min(x), max(x),add=T,col=colors[4],lwd=2.5)
polygon(c(tempdi, rev(tempdi)),c(y1_df, rev(y2_df)),col=alpha(colors[4],alpha=0.2), border=FALSE)
curve(0.59*exp(0.0633*x),-2,40,add=T,col='grey30',lwd=2.5, lty=2)
title(main=expression(bold("Dinoflagellates")), line=-1, adj=0.05, cex=1)
text(-1.6, 2.3, paste0("n=", length(unique(dinos$isolate.code))), adj=c(0,0))
text(-1.6, 1.9, paste0("N=", length(subset(rates, group=="dinoflagellates")$isolate.code)), adj=c(0,0))

mtext(expression(bold("Temperature (ºC)")), side = 1, outer = TRUE, line = -1.5)
mtext(expression(bold("Specific Growth Rate (d"^"-1" *")")), side = 2, outer = TRUE, line = -1.5)
dev.off()

#########################################################################
###### Calculate the Rate of Change for Reaction Norms #########################

# calculate the change in growth estimating change from lower 20% to max growth
growth.change.inc <- vector()
growth.change.dec <- vector()
lower <- vector()
upper <- vector()

for(i in 1:nrow(isolates)){
  o=isolates[i, "mu.c.opt.list"]
  w=isolates[i, "mu.wlist"]
  a=isolates[i, "mu.alist"]
  b=isolates[i, "mu.blist"]
  min=isolates[i, "tmin"]
  max=isolates[i, "tmax"]
  mumax = isolates[i, "mu.g.opt.val.list"]
  topt=isolates[i, "mu.g.opt.list"]
  
  #find 20% of µmax
  target = mumax * 0.20 
  x1=seq(min,topt, by=0.001)
  x2=seq(topt,max, by=0.001)
  
  #find the temperature values that result in the target rates
  lowerbound <- x1[which(abs(nbcurve(x1,o,w,a,b)-target)==min(abs(nbcurve(x1,o,w,a,b)-target)))]  
  upperbound <- x2[which(abs(nbcurve(x2,o,w,a,b)-target)==min(abs(nbcurve(x2,o,w,a,b)-target)))]
  lower<-append(lower, lowerbound)
  upper<-append(upper, upperbound)
  inc = (mumax-target)/(topt-lowerbound) 
  dec = abs((target-mumax)/(upperbound-topt))
  growth.change.inc<-rbind(growth.change.inc, inc)
  growth.change.dec<-rbind(growth.change.dec, dec)
}

isolates$lowerbound = lower
isolates$upperbound = upper
isolates$growth.change.inc = growth.change.inc
isolates$growth.change.dec = growth.change.dec

write.csv(isolates, "output/Isolate_growth_bounds.csv")

##### Figure 2A: Exponential curve comparison ####
inset<- ggplot(data.frame(x = c(-2, 25)), aes(x = x)) +
  stat_function(fun = nbcurve, args=list(8.19,31.2,0.21,0.1), color="grey30")+
  ylim(0, 0.8)+theme_classic()+
  geom_hline(yintercept=0.1566, linetype=2, color="grey60")+
  geom_hline(yintercept=0.1566*5, linetype=2, color="grey60")+
  theme(axis.title = element_blank(), 
        axis.text = element_blank(), axis.ticks = element_blank(),
        plot.margin = unit(c(0, 0.4, 0, 4), "lines"),
        plot.background = element_blank())+
  geom_text(x=-9, y=0.1566, label=expression(µ["20%max"]),
           color="black", size=3)+
  geom_text(x=-7, y=0.1566*5, label=expression(µ["max"]),
            color="black", size=3)+
  coord_cartesian(clip = "off")
inset

rofc<-
  ggplot(data=isolates)+
  geom_boxplot(aes(x=group, y=growth.change.inc), position=position_nudge(x=-0.22),width=0.4, color="black")+
  geom_boxplot(aes(x=group, y=growth.change.dec), fill="grey60",color="black", width=0.4, position=position_nudge(x=0.22))+
  labs(x="", y=expression(bold(paste( "Change in Performance (|", "µ", "|/ºC)"))))+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  guides(fill=FALSE, color=FALSE)+
  scale_x_discrete(labels=c("CO", "CY", "DT", "DF"))+
  annotation_custom(
    ggplotGrob(inset), 
    xmin = 1.75, xmax =Inf, ymin = 0.4, ymax = 0.6)+
  y
rofc


##### Figure 2B: Exponential curve comparison ####
# Color palette
colors <- c("orange", "#ec3a25","#026cb1","#3ea127","grey30", "#033175", "#84e04c", "#fd843d" ) # primary colors
x=seq(-2, 40, by=0.1)

envel<- ggplot(data.frame(x = c(-2, 35)), aes(x = x)) +
  stat_function(fun = dQ, aes(color="Diatoms", linetype="Diatoms"), lwd=0.8, xlim = c(min(tempd), max(tempd)))+
  geom_ribbon(data=data.frame(cbind(tempd, y1, y2)), aes(ymax=y2, ymin=y1, x=tempd), alpha=0.1, fill=colors[[3]])+
  stat_function(fun = cQ, aes(color="Cyanobacteria", linetype="Cyanobacteria"), lwd=0.8, xlim = c(min(tempc), max(tempc)))+
  geom_ribbon(data=data.frame(cbind(tempc, y1_cy, y2_cy)), aes(ymax=y2_cy, ymin=y1_cy, x=tempc), alpha=0.1, fill=colors[[2]])+
  stat_function(fun = coQ, aes(color="Coccolithophores", linetype="Coccolithophores"), lwd=0.8, xlim = c(min(tempco), max(tempco)))+
  geom_ribbon(data=data.frame(cbind(tempco, y1_co, y2_co)), aes(ymax=y2_co, ymin=y1_co, x=tempco), alpha=0.1, fill=colors[[1]])+
  stat_function(fun = diQ, aes(color="Dinoflagellates", linetype="Dinoflagellates"), lwd=0.8, xlim = c(min(tempdi), max(tempdi)))+
  geom_ribbon(data=data.frame(cbind(tempdi, y1_df, y2_df)), aes(ymax=y2_df, ymin=y1_df, x=tempdi), alpha=0.1, fill=colors[[4]])+
  stat_function(fun = ep, aes(color="Eppley (1972)", linetype="Eppley (1972)"), lwd=0.8)+
  labs(x="Temperature (ºC)", y=expression(bold("Specific Growth Rate (d"^"-1" *")")), color="")+
  scale_colour_manual("Groups", values=colors)+
  scale_linetype_manual(values=c(1,1,1,1,2), guide=FALSE)+
  guides(color=guide_legend(override.aes = list(linetype = c(1,1,1,1,2))))+ #overrides color so legend lines are dashed
  coord_cartesian(ylim = c(0.1,2.9))+
  y+theme(legend.position =c(0.20, 0.825), legend.text = element_text(size=9), legend.title = element_blank())
envel


### Save Figure 2 ###
plot_grid(rofc, envel, labels =letters[1:2])
ggsave("figures/Figure2.pdf", width = 8.3, height = 4.2)
