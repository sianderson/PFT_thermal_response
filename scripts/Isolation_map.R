# Code for Anderson et al. (in review, Nature Communications)
# Stephanie I. Anderson updated 09/16/2020

# Contains:
# Map of Strain isolation locations (Extended Figure 1)

# Load Packages
library(rgdal)
library(sp)

# load shape files (downloaded from https://www.naturalearthdata.com/)
shape_path <- "data/maps/"
coast_shapefile <- paste(shape_path, "ne_50m_coastline/ne_50m_coastline.shp", sep="")
ocean_shapefile <- paste(shape_path, "ne_110m_ocean/ne_110m_ocean.shp", sep="")
lakes_shapefile <- paste(shape_path, "ne_110m_lakes/ne_110m_lakes.shp", sep="")
admin0_shapefile <- paste(shape_path, "ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp", sep="")
bb_shapefile <- paste(shape_path, "ne_110m_graticules_all/ne_110m_wgs84_bounding_box.shp", sep="")
grat30_shapefile <- paste(shape_path, "ne_110m_graticules_all/ne_110m_graticules_30.shp", sep="")

layer <- ogrListLayers(coast_shapefile)
coast_lines <- readOGR(coast_shapefile, layer=layer)
 
layer <- ogrListLayers(ocean_shapefile)
ocean_poly <- readOGR(ocean_shapefile, layer=layer)

layer <- ogrListLayers(admin0_shapefile)
admin0_poly <- readOGR(admin0_shapefile, layer=layer)

layer <- ogrListLayers(lakes_shapefile)
lakes_poly <- readOGR(lakes_shapefile, layer=layer)

layer <- ogrListLayers(grat30_shapefile)
grat30_lines <- readOGR(grat30_shapefile, layer=layer)

layer <- ogrListLayers(bb_shapefile)
bb_poly <- readOGR(bb_shapefile, layer=layer)
bb_lines <- as(bb_poly, "SpatialLines")

# set Robinson CRS
robin_crs <- CRS("+proj=robin +lon_0=0w")

bb_poly_proj <- spTransform(bb_poly, robin_crs)
coast_lines_proj <- spTransform(coast_lines, robin_crs)
admin0_poly_proj <- spTransform(admin0_poly, robin_crs)
lakes_poly_proj <- spTransform(lakes_poly, robin_crs)
grat30_lines_proj <- spTransform(grat30_lines, robin_crs)

# convert polygons to spatial lines
lakes_lines_proj <- as(lakes_poly_proj, "SpatialLines")
admin0_lines_proj <- as(admin0_poly_proj, "SpatialLines")
bb_lines_proj <- as(bb_poly_proj, "SpatialLines")

################################################
###### Isolation Data #######
isolates <- read.csv("data/derived_traits.csv")
isolates <- subset(isolates, !is.na(isolation.latitude))

diatom<-subset(isolates, group=="diatoms")
cyano<-subset(isolates, group=="cyanobacteria")
dino<-subset(isolates, group=="dinoflagellates")
coccolithophores<-subset(isolates, group=="coccolithophores")

# Transform Coordinates to Robinson Projection
coordinates(diatom) <- ~isolation.longitude+isolation.latitude
proj4string(diatom) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
robin.crs <- CRS("+proj=robin +lon_0=0w")
diatom_proj <- spTransform(diatom, robin.crs)

coordinates(cyano) <- ~isolation.longitude+isolation.latitude
proj4string(cyano) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
cyano_proj <- spTransform(cyano, robin.crs)

coordinates(coccolithophores) <- ~isolation.longitude+isolation.latitude
proj4string(coccolithophores) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
cocco_proj <- spTransform(coccolithophores, robin.crs)

coordinates(dino) <- ~isolation.longitude+isolation.latitude
proj4string(dino) <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
dino_proj <- spTransform(dino, robin.crs)

################################################
####### Extended Figure 1 ########
################################################
colors <- c("orange", "#ec3a25","#026cb1","#3ea127")

pdf("figures/Extended_Figure1.pdf", width = 4.5, height = 4)
par(mar=c(0,0,0,0))
plot(bb_poly_proj, col="gray90", bor="black", lwd=0.1)
plot(admin0_poly_proj, col="white", bor="white", lwd=0.4, add=TRUE)
plot(grat30_lines_proj, col="black", lwd=0.3, add=TRUE)
plot(coast_lines_proj, col="grey30", lwd=0.5, add=TRUE)
plot(bb_lines_proj, col="black", lwd=1.0, add=TRUE)
plot(cocco_proj, bg=alpha(colors[1], 0.7), col="black",lwd=1.0, add=TRUE, pch=21, cex=1.3)
plot(cyano_proj, bg=alpha(colors[2], 0.7), col="black",lwd=1.0, add=TRUE, pch=22, cex=1.3)
plot(diatom_proj, bg=alpha(colors[3], 0.7), col="black",lwd=1.0, add=TRUE, pch=23, cex=1.2)
plot(dino_proj, bg=alpha(colors[4], 0.7), col="black",lwd=1.0, add=TRUE, pch=24, cex=1.3)

legend(-18000000, -5000000, legend=c("Coccolithophores","Cyanobacteria","Diatoms","Dinoflagellates"), bg="white",
       title=expression(bold("Functional Group")), pch=c(21,22,23,24), pt.lwd=0.6, col="black", pt.bg=colors, cex=0.8)
dev.off()
