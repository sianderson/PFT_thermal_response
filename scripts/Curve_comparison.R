# Code for Anderson et al. (2021) 
# Marine Phytoplankton Functional Types Exhibit Diverse Responses to Thermal Change
# Stephanie I. Anderson updated 09/14/2021

# Contains:
## Statistical Analyses comparing thermal reaction norms
## Extended Figure 2

# Load packages
library(plyr)
library(ggplot2)
library(reshape2)
library(FSA)
library(car)
library(rstatix)
library(ggridges)

# Load data (Output from Thermal_reaction_norms.R)
isolates <- read.csv("output/Isolate_growth_bounds.csv")
source("scripts/custom_theme.R")

############################################################
###### Statistical Analyses ######
############################################################

############ Growth maxima ############
leveneTest(mu.g.opt.val.list ~ group, data=isolates) # different variances
setDT(isolates)[, list(GroupVariance=var(mu.g.opt.val.list)), by = group]

# not normally distributed
isolates %>%  
  group_by(group) %>%
  shapiro_test(mu.g.opt.val.list)

summary(aov(mu.g.opt.val.list ~ group, data=isolates))

# Post-hoc Dunn's test to account for non-normally distributed data
dunnTest(isolates$mu.g.opt.val.list, isolates$group)

# range of umax
ddply(isolates,. (group), summarize, range=max(mu.g.opt.val.list)-min(mu.g.opt.val.list))

############ Niche Width ############
# For entire group
ddply(isolates,. (group), summarize, niche=max(upperbound)-min(lowerbound))

# Average for each group
isolates$niche = isolates$upperbound - isolates$lowerbound
sum <- ddply(isolates,. (group), summarize, mean=mean(niche), n=length(unique(isolate.code)))
cor.test(sum$mean, sum$n) # corelation between niche and sample size

# Statistical Tests
leveneTest(niche ~ group, data=isolates) # different variances
isolates %>% #not normally distributed
  group_by(group) %>%
  shapiro_test(niche)
setDT(isolates)[, list(GroupVariance=var(niche)), by = group]
summary(aov(niche ~ group, data=isolates))
dunnTest(isolates$niche, isolates$group)

# Calculate the number of strains that exist at each temperature for each group using niche
bounds <- isolates[, c('group', 'lowerbound', 'upperbound')]
temp = as.data.frame(seq(-2,40, by=0.1))
temp = cbind(temp, temp)
temp = setNames(data.frame(t(temp[,-1])), temp[,1])
vr <- cbind(bounds, temp)

# For each isolate, mark whether survival would be possible (temperature within niche)
for(i in 4:ncol(vr)){
  vr[i]<-ifelse(vr[i] >= vr$lowerbound & vr[i] <= vr$upperbound, 1, 0)
}
v=vr[,c(1, 4:424)]
v <-melt(v, id='group')
v$variable<-as.numeric(as.character(v$variable))
v<-subset(v, value>0)
v$freq=v$variable*v$value


############################################################
###### Extended Figure 2 ######
############################################################
colors <- c("orange", "#ec3a25","#026cb1","#3ea127","grey30", "#033175", "#84e04c", "#fd843d" ) # primary colors

## Extended Figure 2A --> Viable range
viable<-ggplot(v, aes(x = freq, y=group, group=group)) +
  y+
  scale_y_discrete(limits = rev(c("coccolithophores","cyanobacteria","diatoms","dinoflagellates")), labels=c("DF","DT","CY","CO"))+
  scale_fill_manual("Groups", values = alpha(colors, 0.5))+
  theme(legend.position='bottom')+
  geom_density_ridges(aes(fill = group), scale = 2)+
  labs(x="Temperature (ºC)", y="Density")+
  guides(fill=FALSE)
viable

## Extended Figure 2B --> Growth maxima
mumax<-ggplot(isolates, aes(x=group, y=mu.g.opt.val.list, fill=group))+
  geom_boxplot()+
  scale_fill_manual("Groups", values = alpha(colors, 0.5))+
  guides(fill=FALSE, color=FALSE)+
  y + labs(x="", y=expression(bold("µ"["max"]*" (d"^"-1" *")")))+
  scale_x_discrete(labels=c("CO", "CY", "DT", "DF"))
mumax

# Extended Figure 2
plot_grid(viable, mumax, labels =letters[1:2], rel_widths=c(1, 1))
ggsave('figures/Extended_Figure2.pdf', width = 6.5, height = 3.5)

############ AIC Model Comparison ############
mod1 <- rq(ln.r~temperature+group, data=rates,tau=0.99)
mod2 <- rq(ln.r~temperature*group, data=rates,tau=0.99)
AICctab(mod1, mod2, nobs=length(rates$ln.r), delta=TRUE, sort=TRUE, weights=FALSE, base=TRUE)
