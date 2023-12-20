library(terra)

#### Prepare the layers for the SPIN UP process of the Roth C Model. 
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
plot(AOI)

#The second step is to load the latest version of FAO Soil Organic Carbon map layer (Master Layer), created in script number 0. 
#Open SOC MAP FAO
setwd(WD_SOC)
SOC_MAP_AOI<-rast(c("SOC_MAP_AOI.tif"))
# Next, we will open the clay content layer (from script number 8):
  # Open Clay layer
setwd(WD_CLAY)
Clay_WA_AOI<-rast(c("Clay_WA_AOI.tif"))
Clay_WA_AOI_res<-resample(Clay_WA_AOI,SOC_MAP_AOI,method='bilinear')
Clay_AR_Avg<-crop(Clay_WA_AOI_res,AOI)
Clay_AR_Avg<-mask(Clay_AR_Avg,AOI)
Clay_AR_Avg_res<-resample(Clay_AR_Avg,SOC_MAP_AOI,method='bilinear')
#Arkan's little plot
plot(Clay_AR_Avg_res)
# Next, we will open the climate raster layers (generated in script number 1).  These layers come from the CRU database, but the user can choose local layers if desired, as long as they match the arrangement and format needed for running the model.
#Open Precipitation layer 
setwd(WD_CLIM)
PREC<-rast(c("Prec_Stack_1981-2001_TC.tif"))
PREC_AOI<-crop(PREC,AOI)
PREC_AOI<-resample(PREC_AOI,SOC_MAP_AOI)
PREC_AOI<-mask(PREC_AOI,AOI)
PREC_AOI<-c(PREC_AOI)
#arkan's little plot
plot(PREC_AOI)

#Open Temperatures layer (CRU https://crudata.uea.ac.uk/cru/data/hrg/)
TEMP<-rast(c("Temp_Stack_1981-2001_TC.tif"))
TEMP_AOI<-crop(TEMP,AOI)
TEMP_AOI<-resample(TEMP_AOI,SOC_MAP_AOI)
TEMP_AOI<-mask(TEMP_AOI,AOI)
TEMP_AOI<-c(TEMP_AOI)
#arkan's little plot
plot(TEMP_AOI)

#Open Potential Evapotranspiration layer (CRU https://crudata.uea.ac.uk/cru/data/hrg/)
PET<-rast(c("PET_Stack_1981-2001_TC.tif"))
PET_AOI<-crop(PET,AOI)
PET_AOI<-resample(PET_AOI,SOC_MAP_AOI)
PET_AOI<-mask(PET_AOI,AOI)
PET_AOI<-c(PET_AOI)
#arkan's little plot
plot(PET_AOI)

# Next, we will open, resample and mask the land use raster layer to be used in the spin up phase (representative 1980-2000 period).  In this example we will use the ESA land used reclassified into FAO land use classes (script 9)
# OPen Land Use layer reclassify to FAO classes 
# 0 No Data
# 1 Artificial
# 2 Croplands
# 3 Grassland
# 4 Tree Covered
# 5 Shrubs Covered
# 6 Herbaceous vegetation flooded
# 7 Mangroves
# 8 Sparse Vegetation
# 9 Baresoil
# 10 Snow and Glaciers
# 11 Waterbodies
# 12 TreeCrops
# 13 Paddy fields
setwd(WD_LU)
LU_AOI<-rast(c("ESA_Land_Cover_12clases_FAO_AOI.tif"))
#arkan's little plot
plot(LU_AOI)
writeRaster(LU_AOI,filename=("LU_Stack_Modelling1.tif"))

# Then, we will open the vegetation cover layers (created in script number 7):
  # Open Vegetation Cover layer 
setwd(WD_COV)
Cov_AOI<-rast(c('Cov_stack_AOI.tif'))
#arkan's little plot
plot(Cov_AOI)




# Use Land use layer to convert it to DR layer 
#DPM/RPM (decomplosable vs resistant plant material)
#(1) Most agricultural crops and improved grassland and tree crops 1.44 
#(2) Unimproved grassland and schrub 0.67
#(3) Deciduous and tropical woodland 0.25    
DR<-(LU_AOI==2 | LU_AOI==12| LU_AOI==13)*1.44+ (LU_AOI==4)*0.25 + (LU_AOI==3 | LU_AOI==5 | LU_AOI==6 | LU_AOI==8)*0.67
# Finally, we will create a stack with all the raster layers that have been prepared.
# STACK all layers
Stack_Set_AOI<-c(SOC_MAP_AOI,Clay_WA_AOI_res,TEMP_AOI,PREC_AOI,PET_AOI,DR,LU_AOI,Cov_AOI)
setwd(WD_STACK)
writeRaster(Stack_Set_AOI,filename=("Stack_Set_SPIN_UP_AOI1.tif"))
