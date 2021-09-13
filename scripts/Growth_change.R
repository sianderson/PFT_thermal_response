# Code for Anderson et al. (2021) 
# Marine Phytoplankton Functional Types Exhibit Diverse Responses to Thermal Change
# Stephanie I. Anderson updated 09/16/2020

# Contains:
## Predicted Growth Change by 2100
## Figure 5
## Calculations for cyanobacteria range expansion

# Load Packages/scripts
library(raster)
library(rasterVis)
library(ggplot2)
library(RColorBrewer)
library(ncdf4)
library(sf)
library(data.table)
source("scripts/nbcurve.R")
source("scripts/custom_theme.R")

# Load data
isolates <- read.csv("data/derived_traits.csv")
diatom<-subset(isolates, group=="diatoms")
cyano<-subset(isolates, group=="cyanobacteria")
dino<-subset(isolates, group=="dinoflagellates")
coccolithophores<-subset(isolates, group=="coccolithophores")

## temperature data
past<-raster("SST_data/CMIP_1950_period_mean.grd")
future<-raster("SST_data/CMIP_2100_period_mean.grd")
change<-raster("SST_data/CMIP_period_mean_diff.grd")

#######################################################
# Make Projections for Proportional Growth Change
#######################################################
# Function to remove growth values less than 20% of max growth rate
fun <- function(x, z) {x[x<z] <- NA; return(x) } # where z is 20% of the max growth

# Function to make raster binary (present/absent based on whether growth is possible or not)
binary <- function(x) {
  x[is.na(x)] <- 0
  x[x>0] <- 1
  x} 

# Function to turn no growth into NA
# necessary for calculating extent
fun3 <- function(x) { x[x<=0] <- NA; return(x) }

# Need to run through this for each group separately
pft<-coccolithophores

prop_change<-stack()
gr_now<-stack()
gr_future<-stack()
expansion <- stack()

for(i in 1:nrow(pft)){            
    o=pft[i, "mu.c.opt.list"]
    w=pft[i, "mu.wlist"]
    a=pft[i, "mu.alist"]
    b=pft[i, "mu.blist"]
    mumax = pft[i, "mu.g.opt.val.list"]
    
    # 20% of the max growth rate
    z = mumax*0.2
    
    # past/future growth rates for each strain
    gr<-nbcurve(x=past,opt=o,w=w,a=a,b=b) 
    gr2<-nbcurve(x=future,opt=o,w=w,a=a,b=b) 
    
    # remove growth values less than 20% of max growth rate (where strains wouldn't be viable)
    gr <- calc(gr, function(x){fun(x, z = z)})
    gr2 <- calc(gr2, function(x){fun(x, z = z)})
    
    # calculate proportional change
    gr_change=(gr2-gr)/gr
    
    # create rasters that have just the extent information 
    extent <- calc(gr2, binary) - calc(gr, binary)
    extent <- calc(extent, fun3)
    
    # generate raster stacks
    expansion <- stack(expansion, extent)
    gr_now<-stack(gr_now, gr)
    gr_future<-stack(gr_future, gr2)
    prop_change<-stack(prop_change, gr_change)
  }
# calculate the median growth and extent change for each PFT
mean_change=calc(prop_change, median, na.rm=TRUE)
mean_expansion=calc(expansion, median, na.rm=TRUE)

# Write raster stack for each PFT
writeRaster(mean_expansion, filename = "output/Growth_change/cocco_expansion_1950_2100.grd", overwrite=TRUE)
writeRaster(mean_change, filename = "output/Growth_change/cocco_prop_change_1950_2100.grd", overwrite=TRUE)

# renaming the generated raster stacks (after run with each PFT)
diatom_change=mean_change
diatom_expansion = mean_expansion
cyano_change=mean_change
cyano_expansion = mean_expansion
dino_change=mean_change
dino_expansion = mean_expansion
cocco_change=mean_change
cocco_expansion = mean_expansion

# If no changes necessary, load here
diatom_change<-raster("output/Growth_change/diatom_prop_change_1950_2100.grd")
cyano_change<-raster("output/Growth_change/cyano_prop_change_1950_2100.grd")
dino_change<-raster("output/Growth_change/dino_prop_change_1950_2100.grd")
cocco_change<-raster("output/Growth_change/cocco_prop_change_1950_2100.grd")

diatom_expansion<-raster("output/Growth_change/diatom_expansion_1950_2100.grd")
cyano_expansion<-raster("output/Growth_change/cyano_expansion_1950_2100.grd")
dino_expansion<-raster("output/Growth_change/dino_expansion_1950_2100.grd")
cocco_expansion<-raster("output/Growth_change/cocco_expansion_1950_2100.grd")

# Stack results into new raster stack
mean_change_all<-stack(cocco_change, cyano_change,diatom_change,dino_change) 
mean_expansion_all<-stack(cocco_expansion, cyano_expansion,diatom_expansion,dino_expansion)
names<-c("Coccolithophores","Cyanobacteria", "Diatoms", "Dinoflagellates")
names(mean_change_all)<-names
names(mean_expansion_all)<-names

# Save
writeRaster(mean_change_all, filename = "output/Growth_change/all_prop_change_1950_2100.grd", overwrite=TRUE)

#############################################
# Average growth change by latitude (Figure 4A)
#############################################
# must compute this for each group
mean_change = diatom_change

# data.frame to hold results
diatom_change_dt <- data.table(lat=yFromRow(mean_change, 1:nrow(mean_change)))

# convert to matrices
diatom_change_m <- as.matrix(mean_change)

# average by lat
diatom_change_dt$tos <- apply(diatom_change_m, MARGIN=1, FUN=mean, na.rm=TRUE)

# SD by lat
diatom_change_dt$tossd <- apply(diatom_change_m, MARGIN=1, FUN=sd, na.rm=TRUE)

# Sample size by lat
diatom_change_dt$tos_n <- apply(diatom_change_m, MARGIN=1, FUN=function(x) length(which(!is.na(x))))

# confidence intervals
diatom_change_dt$ci_u<- diatom_change_dt$tos+1.645*(diatom_change_dt$tossd/sqrt(diatom_change_dt$tos_n)) #90%
diatom_change_dt$ci_l<- diatom_change_dt$tos-1.645*(diatom_change_dt$tossd/sqrt(diatom_change_dt$tos_n)) #90%

# save files
write.csv(diatom_change_dt, file='output/Growth_change/change_by_lat_diatom.csv')

# Load files if already made
dino_lat<-read.csv('output/Growth_change/change_by_lat_dino.csv')
diatom_lat<-read.csv('output/Growth_change/change_by_lat_diatom.csv')
cyano_lat<-read.csv('output/Growth_change/change_by_lat_cyano.csv')
cocco_lat<-read.csv('output/Growth_change/change_by_lat_cocco.csv')
diatom_lat$group="diatoms"
dino_lat$group="dinoflagellates"
cyano_lat$group="cyanobacteria"
cocco_lat$group="coccolithophores"
group_lat<-rbind(cocco_lat,cyano_lat, diatom_lat, dino_lat)

colors <- c("orange", "#ec3a25","#026cb1","#3ea127","grey30", "#033175", "#84e04c", "#fd843d" ) # primary colors
group_lat$group = factor(group_lat$group, levels=c('coccolithophores','cyanobacteria','diatoms','dinoflagellates'))

# Figure 4A
trend<-ggplot(data=group_lat, aes(x=lat, y=tos, group=group, color=group))+
  geom_hline(yintercept=0, linetype=2, color="grey70")+
  geom_line()+coord_flip()+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=alpha(colors, 0.4), labels = c("Coccolithophores", "Cyanobacteria", "Diatoms", "Dinoflagellates"))+
  geom_ribbon(aes(x=lat, ymin = tos-tossd, ymax = tos+tossd, fill=group),colour = NA)+y+
  guides(fill=FALSE)+
  scale_y_continuous(breaks=c(-0.5,0,0.5, 1.0), limits = c(-0.55,1.15))+
  scale_x_continuous(breaks=c(-50, 0, 50), labels=c("50ºS", "0º", "50ºN"))+
  labs(y="Proportional \n Growth Change", x="Latitude", fill="")+
  theme(legend.position = c(0.73,0.13),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        legend.text = element_text(size=8),
        legend.key.size = unit(0.4, 'cm'),
        legend.title = element_blank(),
        axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.text.x = element_text(size=10),
        axis.text.y = element_text(size=10))
trend
ggsave("figures/Figure4A.pdf", width = 3, height = 4.2)

#######################################################
# Figure 4B
#######################################################
# Set colors
cutpts <- seq(-0.6,1.2, by=0.1)
colourCount = length(unique(cutpts))
my.palette <- colorRampPalette(brewer.pal(n = 11, name = "RdBu"),bias=0.55)

# Change projection
proj4string(mean_expansion_all) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
proj4string(mean_change_all) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

# Load a world map
world<-st_as_sf(maps::map("world2", fill = TRUE, plot = FALSE, boundary=TRUE, interior=FALSE))
world_outline <- as(st_geometry(world), Class="Spatial")

# Plot
p0 <- levelplot(mean_expansion_all, par.settings=GrTheme)
p1 <- rasterVis::levelplot(mean_change_all, layout=c(2,2),
                    scales=list(),
                     col.regions=rev(my.palette(colourCount)),at=cutpts,
                     xlab=expression(bold("Longitude")),
                     ylab=expression(bold("Latitude")),
                     par.strip.text = list(fontface="bold",cex=0.9),
                     par.settings=list(axis.text=list(cex=0.8),
                                       par.xlab.text=list(cex=0.9),
                                       panel.background= list(col="transparent"),
                                       strip.background = list(col="transparent"),
                                       strip.border = list(col="transparent")))
pdf("figures/Figure4B.pdf", width = 6.5, height = 4.2)
maps <- p1 +as.layer(p0, under = TRUE) + latticeExtra::layer(sp.lines(world_outline, col="grey10",fill="grey10", lwd=0.5))
maps
dev.off()

#############################################
# Calculate range expansion for cyanos
#############################################
# functions to create binary structure, present/absent
fun3 <- function(x) {x[!(is.na(x))] <- 1; return(x) }
fun4 <- function(x) {x[is.na(x)] <- 0; return(x) }

# convert growth area to binary and find new range
b_cyano <- calc(cyano_change, fun3) 
b_cyano <- calc(b_cyano, fun4)
e_cyano <- calc(cyano_expansion, fun4)
new_range = e_cyano-b_cyano 

# Return back to original structure of NAs where strains are absent
fun5 <- function(x) {x[x<1] <- NA; return(x) }
new_range <-calc(new_range, fun5)

# Raster (r) of cell grid sizes (Created in Area_on_Earth.R)
load('output/Area_on_Earth.rdata')

# Calculate the area of cyano expansion
km = new_range * r 
cellStats(km, stat='sum')
km_old = b_cyano*r
cellStats(km_old, stat='sum') 
cellStats(km, stat='sum')/cellStats(km_old, stat='sum') 

