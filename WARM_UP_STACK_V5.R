rm(list = ls())

library(terra)

# Set the number of years of the warm up
nWUP <- 18

WD_AOI<-("H:/BisonRothC/INPUTS/AOI_POLYGON")
WD_SOC<-("H:/BisonRothC/INPUTS/SOC_MAP")
WD_CLAY<-("H:/BisonRothC/INPUTS/CLAY")
WD_CLIM<-("H:/BisonRothC/INPUTS/TERRA_CLIMATE_1")
WD_LU<-("H:/BisonRothC/INPUTS/LAND_USE")
WD_COV<-("H:/BisonRothC/INPUTS/COV")
WD_STACK<-("H:/BisonRothC/INPUTS/STACK")
WD_NPP<-("H:/BisonRothC/INPUTS/NPP")

setwd(WD_AOI)
AOI<-vect("Bialowieza_Forest_3_2.shp")

#Open SOC MAP 
setwd(WD_SOC)
SOC_MAP_AOI<-rast(c("SOC_MAP_AOI.tif"))
# Then we will open the clay layer created in script number 8:
# Open Clay layers  (ISRIC)
setwd(WD_CLAY)
Clay_WA_AOI<-rast(c("Clay_WA_AOI.tif"))
Clay_WA_AOI_res<-resample(Clay_WA_AOI,SOC_MAP_AOI,method='bilinear') 
# arkan's little plot
plot(Clay_WA_AOI_res)


# OPen Land Use layer (ESA)
setwd(WD_LU)
LU_AOI<-rast(c("ESA_Land_Cover_12clases_FAO_AOI.tif"))
#arkan's little plot
plot(LU_AOI)
# Open Vegetation Cover layer 
setwd(WD_COV)
Cov_AOI<-rast(c('Cov_stack_AOI.tif'))

LU_Stack <- rast(c(replicate(nWUP, LU_AOI)))
# Create DR Layer from LU layer (ESA land use , 14 classes)
#DPM/RPM (decomposable vs resistant plant material)
#(1) Most agricultural crops and improved grassland or tree crops 1.44 
#(2) Unimproved grassland and shrub 0.67
#(3) Deciduous and tropical woodland 0.25    
DR<-(LU_AOI==2 | LU_AOI==12| LU_AOI==13)*1.44+ (LU_AOI==4)*0.25 + (LU_AOI==3 | LU_AOI==5 | LU_AOI==6 | LU_AOI==8)*0.67
DR_Stack<-LU_Stack
for (i in 1:length(LU_Stack)){
  DR_Stack[[i]]<-(LU_Stack[[i]]==2 | LU_Stack[[i]]==12)*1.44+ (LU_Stack[[i]]==4)*0.25 + (LU_Stack[[i]]==3 | LU_Stack[[i]]==5 | LU_Stack[[i]]==6 | LU_Stack[[i]]==8)*0.67
}
plot(DR_Stack)
setwd(WD_LU)
writeRaster(LU_AOI,filename=("LU_Stack_Modelling_WU1.tif"))

# STACK all layers
Stack_Set_AOI<-c(SOC_MAP_AOI,Clay_WA_AOI_res,Cov_AOI,LU_Stack[[1]],DR_Stack[[1]]) # Arkan: I put LU_Stack[1] and DR_Stack[1] because they don't change over the years.
# arkan's little plots:
plot(Stack_Set_AOI)
setwd(WD_STACK)
writeRaster(Stack_Set_AOI,filename=("Stack_Set_WARM_UP_AOI1.tif"), overwrite=TRUE)
