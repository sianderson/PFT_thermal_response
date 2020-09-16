# Code for Anderson et al. (in review, Nature Communications)
# Stephanie I. Anderson updated 09/14/2020
# Adapted from Pinsky et al. 2019 (https://doi.org/10.5281/zenodo.2576197)

# Script prepares SST multi-model ensemble data for later use 

# Load Packages
Packages <- c("sdmpredictors", "sp", "raster", "readr", "dplyr",
              "tidyr", "ggplot2", "ncdf4", "ncdf4.helpers", "rasterVis",
              "PCICt", "chron", "lattice", "latticeExtra", "RColorBrewer")
lapply(Packages, library, character.only = TRUE)

###### sst download #########
# Downloaded CMIP5 multi-model ensemble mean on 7/1/19 from http://climexp.knmi.nl/
sstfile <- nc_open("/Volumes/SABackup/SST_data/tos_Omon_modmean_rcp85_ave.nc")

# extract variables
lon <- ncvar_get(sstfile, "lon")
lat <- ncvar_get(sstfile, "lat")
nlon <- dim(lon)
nlat <- dim(lat)
print(c(nlon, nlat))
dlname <- ncatt_get(sstfile,"tos","long_name")
dunits <- ncatt_get(sstfile,"tos","units")
fillvalue <- ncatt_get(sstfile,"tos","_FillValue")

# get global attributes
title <- ncatt_get(sstfile,0,"title")
institution <- ncatt_get(sstfile,0,"institution")
datasource <- ncatt_get(sstfile,0,"source")
history <- ncatt_get(sstfile,0,"history")
Conventions <- ncatt_get(sstfile,0,"Conventions")

# get temperature
tmp_array <- ncvar_get(sstfile,"tos")
dlname <- ncatt_get(sstfile,"tos","long_name")
dunits <- ncatt_get(sstfile,"tos","units")
fillvalue <- ncatt_get(sstfile,"tos","_FillValue")
dim(tmp_array)

# extract time
time <- ncvar_get(sstfile,"time")
head(time)
tunits <- ncatt_get(sstfile,"time","units")
nt <- dim(time)
nt

# convert time -- split the time units string into fields
tustr <- strsplit(tunits$value, " ")
tdstr <- strsplit(unlist(tustr)[3], "-")
tmonth <- as.integer(unlist(tdstr)[2])
tyear <- as.integer(unlist(tdstr)[1])
tday <- as.integer(unlist(tdstr)[3])
chron(time,origin=c(tmonth, tday, tyear))

# replace netCDF fill values with NA's
tmp_array[tmp_array==fillvalue$value] <- NA
length(na.omit(as.vector(tmp_array[,,1])))

# get a single slice or layer (January)
tmp_slice <- tmp_array[,,m] # m=month

infile <- nc_open("/Volumes/SABackup/SST_data/tos_Omon_modmean_rcp85_ave.nc")
tos85 <- ncvar_get(infile, 'tos') 
dim(tos85)
dimnames(tos85) <- list(lon=1:288, lat=1:144, time=1:2880) # tos 288 x 144 x 2800 lon x lat x months 1861-2100
dimnames(tos85)[[1]] <- infile$var[['tos']]$dim[[1]]$vals # add lon
dimnames(tos85)[[2]] <- infile$var[['tos']]$dim[[2]]$vals # add lat
dimnames(tos85)[[3]] <- paste(rep(1861:2100, rep(12, length(1861:2100))), paste('_', 1:12, sep=''), sep='') # year and month of the CMIP5 time dimension
nc_close(infile)

# change dimension order
tos85 <- aperm(tos85, c(2,1,3)) 

 # tos85
 lons <- as.numeric(colnames(tos85))
 lats <- as.numeric(rownames(tos85))
 newlatord <- sort(lats, decreasing=TRUE)
 newrownms <- as.character(newlatord)
 newlatord <- as.character(newlatord)
 tos85 <- tos85[newlatord, ,]

#############################################
# select time periods to analyze
climtime <- as.character(1950:1970) # year of the climatology period
futtime <- as.character(2080:2100) # year the future period

# Reshape CMIP5 to lat x lon x month x year
# tosmax85
nms <- dimnames(tos85)
mos <- as.numeric(unlist(strsplit(nms[[3]], split='_'))[seq(2,length(nms[[3]])*2,by=2)]) # months
yrs <- as.numeric(unlist(strsplit(nms[[3]], split='_'))[seq(1,length(nms[[3]])*2,by=2)]) # years
dm <- dim(tos85)
tos85bymo <- array(tos85, dim=c(dm[1], dm[2], 12, dm[3]/12), dimnames=list(lat=nms[[1]], lon=nms[[2]], mo=1:12, year=sort(unique(yrs)))) # add a month dimension (3rd dimension)

# subset by time period
tospast<-tos85bymo[,,,climtime]
tosfut<-tos85bymo[,,,futtime]

# find the mean for each month (margins = lat, lon, month, year --> we're calculating across all years for each month and location)
meanpast <- apply(tospast, MARGIN=c(1,2,3), FUN=mean)
meanfut <- apply(tosfut, MARGIN=c(1,2,3), FUN=mean)

# find the mean for each time period (margins = lat, lon, month, year --> we're calculating across all years for each location)
meanpast_yr <- apply(tospast, MARGIN=c(1,2), FUN=mean)
meanfut_yr <- apply(tosfut, MARGIN=c(1,2), FUN=mean)

# now find the average distance
mean_diff<-meanfut-meanpast
mean_diff_yr<-meanfut_yr-meanpast_yr

#######################################################
# make some raster bricks and plots
#######################################################
past <- brick(meanpast,ymn=min(lats), ymx=max(lats), xmn=min(lons), xmx=max(lons))
future <- brick(meanfut, ymn=min(lats), ymx=max(lats), xmn=min(lons), xmx=max(lons))
change <- brick(mean_diff, ymn=min(lats), ymx=max(lats), xmn=min(lons), xmx=max(lons))
past_period <- raster(meanpast_yr, ymn=min(lats), ymx=max(lats), xmn=min(lons), xmx=max(lons))
fut_period <- raster(meanfut_yr, ymn=min(lats), ymx=max(lats), xmn=min(lons), xmx=max(lons))
diff_period <- raster(mean_diff_yr, ymn=min(lats), ymx=max(lats), xmn=min(lons), xmx=max(lons))

# Convert to celcius
past<-past-273.15
future<-future-273.15
past_period<-past_period-273.15
fut_period<-fut_period-273.15
change = fut_period-past_period

# Save 
writeRaster(past_period, "SST_data/CMIP_1950_period_mean.grd")
writeRaster(fut_period, "SST_data/CMIP_2100_period_mean.grd")
writeRaster(change, "SST_data/CMIP_period_mean_diff.grd")

#############################################
# Average temperature change by latitude
#############################################
library(data.table)
# data.frame to hold results
diff_period_dt <- data.table(lat=yFromRow(diff_period, 1:nrow(diff_period)))

# convert to matrices
diff_period_m <- as.matrix(diff_period)

# average by lat
diff_period_dt$tos <- apply(diff_period_m, MARGIN=1, FUN=mean, na.rm=TRUE)

# plot to check
diff_period_dt[,plot(lat, tos, type='l', ylim=c(0,10))]

# SD by lat
diff_period_dt$tossd <- apply(diff_period_m, MARGIN=1, FUN=sd, na.rm=TRUE)

# Sample size by lat
diff_period_dt$tos_n <- apply(diff_period_m, MARGIN=1, FUN=function(x) sum(!is.na(x)))

# write out
write.csv(diff_period_dt, file='output/change_by_lat.csv')


#############################################
# Extended Figure 4
#############################################
change<-raster("SST_data/CMIP_period_mean_diff.grd")

cutpts <- seq(-0.5,6, by=0.5)
colourCount = length(unique(cutpts))
my.palette <- colorRampPalette(brewer.pal(n = 11, name = "Reds"), bias=0.9)

crs(change)<-"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
fun <- function(x) {x[is.na(x)] <- 100; return(x) }
fun2 <- function(x) {x[x<100] <- NA; return(x) }

# raster of the continents
continents<- calc(diff_period, fun)
continents<- calc(continents, fun2)
crs(continents)<-"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

pdf("figures/Extended_Figure4.pdf", width = 6.2, height = 5)
p0<-rasterVis::levelplot(change, col.regions=my.palette(colourCount),at=cutpts,
                         xlab=expression(bold("Longitude")),
                         ylab=expression(bold("Latitude")),
                         par.strip.text = list(fontface="bold", cex=0.9),
                         par.settings=list(axis.text=list( cex=0.8),
                                           par.xlab.text=list( cex=0.9),
                                           panel.background= list(col="transparent"),
                                           strip.background = list(col="transparent"),
                                           strip.border = list(col="transparent")),
                         colorkey=list(title=expression("Change \n in SST \n (ºC)"), row=3, column=1, vjust=2))
p1 <- levelplot(continents, col.regions="black")
p0 + as.layer(p1, under = TRUE)
dev.off()
