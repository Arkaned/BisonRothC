library(terra)
library(tmap)
library(ggplot2)

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
plot(SOC_MAP_AOI, col=gray((1000:1)/1000), plg=list(cex.legend=2), pax = list(cex.axis = 2, cex.lab=2))
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

# colours for plot:
colors <- data.frame(value= c(1:13),
                     col= c("grey",
                            "darkorange",
                            "green",
                            "darkgreen",
                            "brown",
                            "cyan",
                            "blue",
                            "green",
                            "beige",
                            "grey",
                            "darkblue",
                            "forestgreen",
                            "pink"))




plot(LU_AOI, box=FALSE, axes=FALSE, legend=TRUE, col=colors)

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

plot_labels <- rast(c("ESA_Land_Cover_12clases_FAO_AOI.tif"))
label_values <- c(1:13)
label_List <- c('Artificial', 'Croplands', 'Grassland', 'Tree Covered', 'Shrubs Covered', 'Herbaceous flood vegetation', 'Mangroves', 'Sparse vegetation', 'Bare soil', 'Snow and Glaciers', 'Water bodies', 'Tree crops', 'Paddy Fields')

class_Dat <- data.frame(label_values, label_List)

classify(x=plot_labels, rcl=class_Dat)







# Then, we will open the vegetation cover layer created in script number 7.
# Open Vegetation Cover 
setwd(WD_COV)
Cov_AOI<-rast(c('Cov_stack_AOI.tif'))
plot(Cov_AOI)

##### attempt to visualize this nicely:

## first would need to be melted:
Cov_AOI_plot <- terra::rast(c('Cov_stack_AOI.tif'))
cov_Wrangled <- as.data.frame(Cov_AOI_plot, xy=TRUE) %>%
  #--- remove cells with NA for any of the layers ---#
  na.omit()
  #--- change the variable names --
#-- take a look --#

head(cov_Wrangled)

for (i in 3:14) {
  colnames(cov_Wrangled)[i] <- paste("Month", as.character(i-2))
}

melted_cov <- melt(cov_Wrangled, id=c("x", "y"))

## now map it (one map):
(
  g_tmax_map <- ggplot(data = melted_cov) +
    geom_raster(aes(x = x, y = y, fill=value )) +
    
    scale_fill_viridis_c() +
    theme_void() +
    theme(
      legend.position = "bottom"
    )
)

# ## now map it (multiplemap):
ggplot(data=melted_cov) +
  geom_raster(data = melted_cov, aes(x = x, y = y, fill = value)) +
  facet_wrap(variable ~ .) +
  coord_equal() +
  scale_fill_viridis_c() +
  theme_void() +
  theme(
    legend.position = "bottom"
  )


#############
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
terra::writeRaster(Stack_Set_AR,filename=("Stack_Set_FOWARD1.tif"))
