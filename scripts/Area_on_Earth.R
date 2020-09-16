# Code for Anderson et al. (in review, Nature Communications)
# Stephanie I. Anderson updated 09/14/2020

# function to calculate area on Earth (latitude-longitude quandrangle)
# lats and lons from SST_multimodel_ensemble.R

# input is in degrees
# output is in km^2

surface <-
  function(lat1, lon1, lat2, lon2){
    rad<-pi/180
    A1<-2*pi*(1-sin(lat1*rad))*(6371^2)
    A2<-2*pi*(1-sin(lat2*rad))*(6371^2)
    surface<-abs(A1-A2)*abs(lon1-lon2)/360
    surface
  }

surfacearea<-data.frame()
for (i in 1:(length(lats)+1)){
  lat1=lats[i]-0.625
  lat2=lats[(i-1)]-0.625
  lon1=lons[2]-0.625
  lon2=lons[1]-0.625
  sa=surface(lat1, lon1, lat2, lon2)
  surfacearea=rbind(surfacearea, sa)
}
surfacearea[144,1]=surfacearea[1,1]
surfacearea<-as.vector(surfacearea$X210.731826067799)
surfacearea<-matrix(surfacearea,nrow=length(surfacearea),ncol=288)

#generate a raster that can be used for future math
r <-raster(
  surfacearea,
  xmn=min(lons), xmx=max(lons),
  ymn=min(lats), ymx=max(lats)
)
plot(r)
save(r, file='output/Area_on_Earth.rdata')

