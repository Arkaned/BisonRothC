library(terra)
library(sf)
library(sp)

# https://files.isric.org/soilgrids/former/2017-03-10/data/ (clay layers collected from here)

WD_AOI<-("H:/BisonRothC/INPUTS/AOI_POLYGON")
WD_ISRIC<-("D:/CLAY_LAYERS")
WD_CLAY<-("H:/BisonRothC/INPUTS/CLAY")
# Open the shapefile of the region/country
setwd(WD_AOI)
AOI<-read_sf("Bialowieza_Forest_3_2.shp")
# Open Clay layers  (ISRIC)
setwd(WD_ISRIC)
Clay1<-rast(c("CLYPPT_M_sl1_250m_ll.tif")) # 0-1 cm
Clay2<-rast(c("CLYPPT_M_sl2_250m_ll.tif")) # 1-5 cm
Clay3<-rast(c("CLYPPT_M_sl3_250m_ll.tif")) # 5-15 cm
Clay4<-rast(c("CLYPPT_M_sl4_250m_ll.tif")) # 15-30 cm
Clay1_AR<-crop(Clay1,AOI)
Clay2_AR<-crop(Clay2,AOI)
Clay3_AR<-crop(Clay3,AOI)
Clay4_AR<-crop(Clay4,AOI)
# Average of four depths 
WeightedAverage<-function(r1,r2,r3,r4){return(r1*(1/30)+r2*(4/30)+r3*(10/30)+r4*(15/30))}
Clay_WA<-lapp(x=c(Clay1_AR,Clay2_AR,Clay3_AR,Clay4_AR),fun=WeightedAverage)
Clay_WA_AOI<-mask(Clay_WA,AOI)
setwd(WD_CLAY)
writeRaster(Clay_WA_AOI,filename="Clay_WA_AOI.tif")


# testing individual file/levels

writeRaster(Clay1_AR,filename="Clay_1_test.tif", overwrite = TRUE)
writeRaster(Clay2_AR,filename="Clay_2_test.tif", overwrite = TRUE)
writeRaster(Clay3_AR,filename="Clay_3_test.tif", overwrite = TRUE)
writeRaster(Clay4_AR,filename="Clay_4_test.tif", overwrite = TRUE)
