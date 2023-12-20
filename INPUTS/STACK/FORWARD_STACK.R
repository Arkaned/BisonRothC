library(terra)

rm(list = ls())
WD_AOI<-("H:/BisonRothC/INPUTS/AOI_POLYGON")
WD_SOC<-("H:/BisonRothC/INPUTS/SOC_MAP")
WD_CLAY<-("H:/BisonRothC/INPUTS/CLAY")
WD_CLIM<-("H:/BisonRothC/INPUTS/TERRA_CLIMATE_1")
WD_LU<-("H:/BisonRothC/INPUTS/LAND_USE")
WD_COV<-("H:/BisonRothC/INPUTS/COV")
WD_STACK<-("H:/BisonRothC/INPUTS/STACK")

# Open the shapefile of the region/country
setwd(WD_AOI)
AOI<-vect("Bialowieza_Forest_3_2.shp")
# Then, we will open the SOC layer and the clay layer.
#Open SOC MAP 
setwd(WD_SOC)
SOC_MAP_AOI<-rast(c("SOC_MAP_AOI.tif"))
# Open Clay layers  (ISRIC)
setwd(WD_CLAY)
Clay_WA_AOI<-rast(c("Clay_WA_AOI.tif"))
Clay_WA_AOI_res<-resample(Clay_WA_AOI,SOC_MAP_AOI,method='bilinear') 
#arkan's little plots:
plot(Clay_WA_AOI_res)

#Then we will open the 2000-2020 average climate layers created (as the one created in script number 2)
#Open Precipitation layer (CRU https://crudata.uea.ac.uk/cru/data/hrg/)
setwd(WD_CLIM)
PREC<-rast(c("Prec_Stack_01-19_TC.tif"))
PREC_AOI<-crop(PREC,AOI)
PREC_AOI<-resample(PREC_AOI,SOC_MAP_AOI)
PREC_AOI<-mask(PREC_AOI,AOI)
PREC_AOI<-c(PREC_AOI)
# arkan's little plots:
plot(PREC_AOI)

#Open Temperatures layer (CRU https://crudata.uea.ac.uk/cru/data/hrg/)
TEMP<-rast(c("Temp_Stack_01-19_TC.tif"))
TEMP_AOI<-crop(TEMP,AOI)
TEMP_AOI<-resample(TEMP_AOI,SOC_MAP_AOI)
TEMP_AOI<-mask(TEMP_AOI,AOI)
TEMP_AOI<-c(TEMP_AOI)
#arkan's little plots:
plot(TEMP_AOI)
#Open Potential Evapotranspiration layer (CRU https://crudata.uea.ac.uk/cru/data/hrg/)
PET<-rast(c("PET_Stack_01-19_TC.tif"))
PET_AOI<-crop(PET,AOI)
PET_AOI<-resample(PET_AOI,SOC_MAP_AOI)
PET_AOI<-mask(PET_AOI,AOI)
PET_AOI<-c(PET_AOI)
#arkan's little plots:
plot(PET_AOI)
# Then, we will open the land use layer (latest available year) created in script number 10.
setwd(WD_LU)
LU_AOI<-rast(c("ESA_Land_Cover_12clases_FAO_AOI.tif"))
plot(LU_AOI)
# Then, we will open the vegetation cover layer created in script number 7.
# Open Vegetation Cover 
setwd(WD_COV)
Cov_AOI<-rast(c('Cov_stack_AOI.tif'))
plot(Cov_AOI)
# As in the previous scripts, this script creates a DR layer (DPM/RPM ratio), assigning default DPM/RPM ratios to each FAO land cover class (See Table 9.13). Users can modify these ratios according to local expertise and available local information. 
# Open Land use layer and convert it to DR layer (mod 12 , 14 classes)
#DPM/RPM (decomplosable vs resistant plant material...como se divide los C inputs)
#(1) Most agricultural crops and improved grassland or tree crops 1.44 
#(2) Unimproved grassland and schrub 0.67
#(3) Deciduous and tropical woodland 0.25    
DR<-(LU_AOI==2 | LU_AOI==12| LU_AOI==13)*1.44+ (LU_AOI==4)*0.25 + (LU_AOI==3 | LU_AOI==5 | LU_AOI==6 | LU_AOI==8)*0.67




# STACK all layers
Stack_Set_AR<-c(SOC_MAP_AOI,Clay_WA_AOI_res,TEMP_AOI,PREC_AOI,PET_AOI,DR,LU_AOI,Cov_AOI)
setwd(WD_STACK)
writeRaster(Stack_Set_AR,filename=("Stack_Set_FOWARD.tif"), overwrite=TRUE)
