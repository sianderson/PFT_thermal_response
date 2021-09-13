# Code for Anderson et al. (2021) 
# Marine Phytoplankton Functional Types Exhibit Diverse Responses to Thermal Change
# Stephanie I. Anderson updated 09/14/2020

# Contains:
## Modeled temperatures for each strain based on isolation locations
## Calculation of thermal capacities
## Figure 3
## Extended Figure 6

# Load Libraries
library(raster)
library(ggplot2)

# Load Scripts
source("scripts/custom_theme.R")
source("scripts/nbcurve.R")

## temperature data
past_period<-raster("SST_data/CMIP_1950_period_mean.grd")
fut_period<-raster("SST_data/CMIP_2100_period_mean.grd")
diff_period<-raster("SST_data/CMIP_period_mean_diff.grd")

#######################################################
# Extract temperatures for each isolate of known origin
#######################################################
# Load data
isolates <- read.csv("data/derived_traits.csv")

# Remove strains of unknown origin
isolates <- subset(isolates, !is.na(isolation.latitude))

# Remove strains that had a poor quality modeled thermal maxima
isolates <- subset(isolates, maxqual=="good")

# correct for locations that are too far inside estuary
isolates$isolation.longitude.cor = ifelse(isolates$isolation.longitude==-123.17, 236.67, 
                                       ifelse(isolates$isolation.longitude<0, isolates$isolation.longitude+360,
                                              isolates$isolation.longitude))
isolates$isolation.latitude.cor = ifelse(isolates$isolation.latitude>59.4 & isolates$isolation.latitude< 59.8, 58.9, 
                                      ifelse(isolates$isolation.latitude==38.22, 37.8, 
                                             isolates$isolation.latitude))

# Extract modeled temperatures for each isolation location (1950-1970 & 2080-2100)
for(i in 1:nrow(isolates)) {
  location <- data.frame(isolates[,35],isolates[,36]) 
  Thab.past<-extract(x=past_period, y=location, method="bilinear")
  Thab.future<-extract(x=fut_period, y=location, method="bilinear")
}
isolates$Thab.past<-Thab.past
isolates$Thab.fut<-Thab.future

#######################################################
# Warming tolerance and thermal safety margin
#######################################################
# load libraries
library(plyr)
library(grid)
library(gtable)
library(gridExtra)
library(cowplot)
library(ggpubr)
library(tidyverse)
library(egg)

colors <- c("orange", "#ec3a25","#026cb1","#3ea127","grey30", "#033175", "#84e04c", "#fd843d" )
isolates$safety_margin = isolates$mu.g.opt.list - isolates$Thab.past
isolates$warming_tolerance = isolates$tmax-isolates$Thab.past

# Load data
temp_change<-read.csv("output/change_by_lat.csv") #output created in Thermal_reaction_norms.R

# Calculate the distance to the growth equivalence (DGE)
x=seq(0,50, by=0.01)
growth.eq <- data.frame()
for(i in 1:nrow(isolates)){
  o=isolates[i, "mu.c.opt.list"]
  w=isolates[i, "mu.wlist"]
  a=isolates[i, "mu.alist"]
  b=isolates[i, "mu.blist"]
  topt=isolates[i, "mu.g.opt.list"]
  thab=isolates[i, "Thab.past"]
 
   # Find the point above the thermal optimum (Topt) that results in the same growth rate
  ge<-list(uniroot(function(x) nbcurve(x,o,w,a,b)-nbcurve(thab,o,w,a,b),interval=c(topt, 50))$root)
  growth.eq<-rbind(growth.eq, ge)
}
names(growth.eq)<-"growth.eq"
isolates<-cbind(isolates,growth.eq)

# Remove strains that are already inhabiting temperatures above Topt
isolates$equivalence = ifelse(isolates$Thab.past<isolates$mu.g.opt.list, 
                           (isolates$growth.eq - isolates$Thab.past), 0)

# Summary statistics for each trait
sum_traits<-ddply(isolates, .(group), summarise, 
                  avg_wt=mean(warming_tolerance), avg_sm=mean(safety_margin), avg_ge=mean(equivalence))

# Predicted change in habitat temperature
isolates$Thab.change=isolates$Thab.fut-isolates$Thab.past

# Write csv to output folder
write.csv(isolates, file = "output/Thermal_capacity_by_group.csv")

#######################################################
# Figure 3
#######################################################
# Load data
isolates <- read.csv("output/Thermal_capacity_by_group.csv")

# Reorder groups
isolates$group <- factor(isolates$group, levels = c("coccolithophores", "cyanobacteria", "diatoms", "dinoflagellates"))

# Warming Tolerance
dist<-ggplot()+ 
  geom_hline(yintercept=0, linetype="dashed", color = "grey50")+
  geom_violin(data=isolates,aes(x=group, y=warming_tolerance, fill=group))+
  scale_fill_manual(values=colors)+
  labs(y="ºC", x= "")+ guides(fill=FALSE)+
  y+ theme(axis.title.y = element_text(size=9),
           axis.text.x=element_text(size=8),
           axis.text.y=element_text(size=8),
           rect = element_rect(fill = "transparent")) +
  scale_y_continuous(position = "right",breaks=c(0,10,20))+
  scale_x_discrete(labels=c("CO","CY", "DT", "DF"))+
  geom_point(data=sum_traits,
             aes(x=group, y=avg_wt, colour="grey30"), colour="grey20")
dist

wt<-ggplot()+
  geom_ribbon(data=temp_change, aes(x=lat, ymin = tos-tossd, ymax = tos+tossd), fill = "grey70")+
  geom_line(data=temp_change, aes(x=lat, y=tos), color="grey30", size=0.75)+
  geom_point(data=isolates, 
             aes(x=isolation.latitude.cor, y=warming_tolerance, color=group, shape=group), 
             size=2.5, alpha=0.7)+
  scale_shape_manual(values=c(15, 18, 17, 16))+
  scale_color_manual(values=colors)+
  scale_x_continuous(breaks=c(-50, 0, 50), labels=c("50ºS", "0º", "50ºN"))+
  y+labs(y="Warming Tolerance (ºC)", x= "Isolation Latitude",
         color="Functional \n Group", shape="Functional \n Group")+
  geom_hline(yintercept=0, linetype="dashed", color = "grey50")+
  annotation_custom(ggplotGrob(dist), xmin = -95, xmax = 8, ymin = 20, ymax = 38)+
  guides(color=FALSE, shape=FALSE)+
  theme(axis.title.x = element_blank())+
  ylim(-10,36)
wt

# Thermal Safety margin
sm_dist<-ggplot()+ 
  geom_hline(yintercept=0, linetype="dashed", color = "grey50")+
  geom_violin(data=isolates, aes(x=group, y=safety_margin, fill=group))+
  scale_fill_manual(values=colors)+
  labs(y="ºC", x= "")+ guides(fill=FALSE)+
  y+ theme(axis.title.y = element_text(size=9),
           axis.text.x=element_text(size=8),
           axis.text.y=element_text(size=8),
           rect = element_rect(fill = "transparent")) +
  scale_y_continuous(position = "right")+
  scale_x_discrete(labels=c("CO","CY", "DT", "DF"))+
  geom_point(data=sum_traits,
             aes(x=group, y=avg_sm, colour="grey30"), colour="grey20")
sm_dist

sm<-ggplot()+
  geom_ribbon(data=temp_change, aes(x=lat, ymin = tos-tossd, ymax = tos+tossd), fill = "grey70")+
  geom_line(data=temp_change, aes(x=lat, y=tos), color="grey30", size=0.75)+
  geom_point(data=isolates,
             aes(x=isolation.latitude.cor, y=safety_margin, color=group, shape=group), 
             size=2.5, alpha=0.7)+ ylim(-10,36)+
  scale_shape_manual(values=c(15, 18, 17, 16))+
  scale_color_manual(values=colors)+
  scale_x_continuous(breaks=c(-50, 0, 50), labels=c("50ºS", "0º", "50ºN"))+
  y+labs(y="Thermal Safety Margin (ºC)", x= "Isolation Latitude",
         color="", shape="")+
  geom_hline(yintercept=0, linetype="dashed", color = "grey50")+
  annotation_custom(ggplotGrob(sm_dist), xmin = -95, xmax = 8, ymin = 20, ymax = 38)+
  theme(axis.title.x = element_blank())+
  guides(color=FALSE, shape=FALSE)
sm

# Growth Equivalence
ge_dist<-ggplot()+ 
  geom_hline(yintercept=0, linetype="dashed", color = "grey50")+
  geom_violin(data=isolates,
              aes(x=group, y=equivalence, fill=group))+
  scale_fill_manual(values=colors)+
  labs(y="ºC", x= "")+ guides(fill=FALSE)+
  y+ theme(axis.title.y = element_text(size=9),
           axis.text.x=element_text(size=8),
           axis.text.y=element_text(size=8),
           rect = element_rect(fill = "transparent")) +
  scale_y_continuous(position = "right",breaks=c(0,10,20))+
  scale_x_discrete(labels=c("CO","CY", "DT", "DF"))+
  geom_point(data=sum_traits,
             aes(x=group, y=avg_ge, colour="grey30"), colour="grey20")
ge_dist

ge<-ggplot()+
  geom_ribbon(data=temp_change, aes(x=lat, ymin = tos-tossd, ymax = tos+tossd), fill = "grey70")+
  geom_line(data=temp_change, aes(x=lat, y=tos), color="grey30", size=0.75)+
  geom_point(data=isolates, 
             aes(x=isolation.latitude.cor, y=equivalence, color=group, shape=group), 
             size=2.5, alpha=0.7)+ 
  scale_shape_manual(guide=guide_legend(nrow=2),values=c(15, 18, 17, 16))+
  scale_color_manual(guide=guide_legend(nrow=2),values=colors)+
  scale_x_continuous(breaks=c(-50, 0, 50), labels=c("50ºS", "0º", "50ºN"))+
  y+labs(y="Distance to µ Equivalence (ºC)", x= "Isolation Latitude",
         color="", shape="")+ylim(-10,36)+
  geom_hline(yintercept=0, linetype="dashed", color = "grey50")+
   annotation_custom(
     ggplotGrob(ge_dist), 
     xmin = -95, xmax = 8, ymin = 20, ymax = 38)+
  theme(legend.position = c(0.5, 0.15),
        axis.title.x = element_blank(),
        legend.text = element_text(size=9),
        legend.spacing.x = unit(0.1, 'cm'))
ge

# Add space for labeled diagram
blank <- ggplot() + theme(panel.background = element_rect(fill = "white", colour="white"))
blank

gA <- ggplotGrob(sm)
gB <- ggplotGrob(ge)
gC <- ggplotGrob(wt)
maxHeight = grid::unit.pmax(gA$heights[2:5], gB$heights[2:5], gC$heights[2:5])
gA$heights[2:5] <- as.list(maxHeight)
gB$heights[2:5] <- as.list(maxHeight)
gC$heights[2:5] <- as.list(maxHeight)
g1 <- ggplot_gtable(ggplot_build(ge))
Tfinal <-grid.arrange(blank, gA,gB,gC, ncol=2, bottom = textGrob(expression(bold("Isolation Latitude"))))

p<-as_ggplot(Tfinal) + draw_plot_label(label = c("a", "b", "c", "d"), size = 15,
                                   x = c(0, 0.5, 0, 0.5), y = c(1, 1, 0.5, 0.5)) # Add labels
p

ggsave("figures/Figure3.pdf", width = 7, height = 7)

#######################################################
# Extended Figure 6
#######################################################

ggplot(data=isolates)+
  geom_smooth(aes(x=Thab.past, y=Thab.change+Thab.past), color="grey30", span=0.3)+
  geom_point(aes(x=Thab.past, y=tmax, color=group, shape=group, group=group), size=2.5, alpha=0.8)+
  scale_color_manual(values=alpha(colors,0.4))+
  scale_shape_manual(values=c(15, 18, 17, 16))+
  y+labs(y="Thermal Maximum (ºC)", x= "Habitat Temperature (ºC)",
         color="", shape="")+  theme(legend.position = c(0.77, 0.23),
                                     legend.text = element_text(size=10))+
  geom_abline(intercept = 0, slope = 1, linetype=2)+guides(line=FALSE)

ggsave("figures/Extended_Figure6.pdf", width = 4, height = 4)

#######################################################
# Significance Tests
###################################################
library(rstatix)
library(car)
library(onewaytests)
library(FSA)
library(data.table)

######## Warming Tolerance #########
# Variances not significantly different
leveneTest(warming_tolerance ~ group, data=isolates)

# not normally distributed
isolates %>% 
  group_by(group) %>%
  shapiro_test(warming_tolerance)

# With equal variance and NOT normal distribution of residuals, we used Kruskal-Wallis Test
kruskal.test(warming_tolerance ~ group, data=isolates)

######## Thermal Safety Margin #########
leveneTest(safety_margin ~ group, data=isolates)
setDT(isolates)[, list(GroupVariance=var(safety_margin)), by = group] #variance of highest < 4x lowest

subset(isolates) %>% 
  group_by(group) %>%
  shapiro_test(safety_margin)
# normally distributed

kruskal.test(safety_margin ~ group, data=isolates)
dunnTest(isolates$safety_margin, isolates$group) # Unadjusted p-values

######## Distance to Growth Equivalence #########
leveneTest(equivalence ~ group, data=isolates)
isolates %>% 
  group_by(group) %>%
  shapiro_test(equivalence)
# not normally distributed groups

# Kruskal-Wallis test for non-normal data
kruskal.test(equivalence ~ group, data=isolates)

