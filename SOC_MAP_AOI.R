# soc layer:

#Install all necessary packages
install.packages(c("raster","rgdal","SoilR","Formula","soilassessment","abind","ncdf4"))

#Load the packages into R
library(raster)
library(sf) # works!
library(terra)
library(rgdal)

# Set the path to GSOCmap and Area of interest (AOI) vector.
WD_AOI<-("H:/BisonRothC/INPUTS/AOI_POLYGON")
WD_GSOC<-("H:/BisonRothC/INPUTS/SOC_MAP")

# Open the shapefile of the AOI (region/country)
setwd(WD_AOI)
#  AOI<-readOGR("Departamento_Pergamino.shp") OUTDATED CODE
AOI <- read_sf('Bialowieza_Forest_3_2.shp')

#Open FAO GSOC MAP 
setwd(WD_GSOC)
SOC_MAP<-rast("GSOCmap1.5.0.tif")

SOC_MAP_AOI<-crop(SOC_MAP,AOI)
SOC_MAP_AOI<-mask(SOC_MAP_AOI,AOI)
writeRaster(SOC_MAP_AOI,filename="SOC_MAP_AOI.tif",filetype="GTiff")
