library(terra)
library(sf)
library(sp)
library(raster)

# Points to Raster

# MSc Ing Agr Luciano E Di Paolo
# Dr Ing Agr Guillermo E Peralta

rm(list=ls()) 


WD_F<-("H:/BisonRothC/OUTPUTS/3_FORWARD/FOWARD_AOI_final")

WD_SOC<-("H:/BisonRothC/INPUTS/SOC_MAP")

WD_AOI<-("H:/BisonRothC/INPUTS/AOI_POLYGON")

WD_MAPS<-("H:/BisonRothC/OUTPUTS/4_MAPS")

#Define the name of the Country ("ISO3CountryCode")
name<-"Bialowieza" 

#Open FORWARD vector

setwd(WD_F)
FOWARD<-read_sf("FOWARD_AOI_final.shp")


#####################################################################################################################################################

foward_vis <- terra::vect("FOWARD_AOI_final.shp")


Max_BAU2 <- c()
Min_BAU2 <- c()
Med_BAU2 <- c()

Max_SSM2 <- c()
Med_SSM2 <- c()
Min_SSM2 <- c()


SUM_Max_BAU2 <- 0
SUM_Min_BAU2 <- 0
SUM_Med_BAU2 <- 0
SUM_Max_SSM2 <- 0
SUM_Med_SSM2 <- 0
SUM_Min_SSM2 <- 0




for (val in 1:length(foward_vis)){
  if(is.na(foward_vis[,2][[1]][[1]][val]) | foward_vis[,2][[1]][[1]][val] == 0) {next}
  
  Max_BAU2[val] <- foward_vis[, 14][[1]][[1]][val]
  Min_BAU2[val] <- foward_vis[, 13][[1]][[1]][val]
  Med_BAU2[val] <- foward_vis[, 3][[1]][[1]][val]
  Max_SSM2[val] <- foward_vis[, 16][[1]][[1]][val]
  Med_SSM2[val] <- foward_vis[, 11][[1]][[1]][val]
  Min_SSM2[val] <- foward_vis[, 15][[1]][[1]][val]
  
  
  # now calculate the difference between everything and such without bison...
  SUM_Max_BAU2 <- SUM_Max_BAU + (Max_BAU[val] - Med_BAU[val])
  SUM_Min_BAU2 <- SUM_Min_BAU + (Min_BAU[val] - Med_BAU[val])
  SUM_Max_SSM2 <- SUM_Max_SSM + (Max_SSM[val] - Med_BAU[val])
  #SUM_Med_SSM <- SUM_Med_SSM + (Med_SSM[val] - Med_BAU[val])
  SUM_Min_SSM2 <- SUM_Min_SSM + (Min_SSM[val] - Med_BAU[val])
  
  if (val%%100 == 0) {print(val)}
}


sum_Diff2 <- data.frame(SUM_Max_BAU2,
                       
                       SUM_Max_SSM2,
                       SUM_Min_BAU2,
                       #SUM_Med_SSM,
                       SUM_Min_SSM2
                       #Max_BAU,
                       #Max_SSM,
                       
                       #Min_BAU,
                       #Min_SSM
)
diff_Dat2 <- data.frame(Max_BAU2,
                       Max_SSM2,
                       Med_BAU,
                       Med_SSM,
                       Min_BAU2,
                       Min_SSM2)


melted_sum2 <- melt(sum_Diff2)

ggplot(melted_sum2, aes(x=variable, y=value)) +
  geom_bar(stat="identity")


melted_diff2 <- melt(diff_Dat2)

ggplot(melted_diff2) + 
  geom_boxplot(aes(x=variable, y=value, color=variable, fill=variable)) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(mapping=aes(x=variable, y=value), color="black", size=0.1, alpha=0.1) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("A boxplot with jitter") +
  xlab("")





#####################################################################################################################################################

#Open SOC MAP (master layer)

setwd(WD_SOC)
SOC_MAP<-raster::raster(c("SOC_MAP_AOI.tif"))

#Creates emtpy raster 

empty_raster<-SOC_MAP*0

# Open the country vector boundaries

setwd(WD_AOI)
Country<-read_sf(c("Bialowieza_Forest_3_2.shp"))

# Cut the raster with the country vector
Country_raster<-crop(empty_raster,Country)

# Arkan: maybe problem is that raster has NA's?
#Country_raster[!is.na(Country_raster)]<-1
#Country_raster[is.na(Country_raster)]<-0
# Replace Na values for zero values

FOWARD[is.na(FOWARD)] <- 0

#FOWARD$SOC_BAU_20[is.na(FOWARD$SOC_BAU_20)]<-0
#FOWARD$Low_Scena0[is.na(FOWARD$Low_Scena0)]<-0
#FOWARD$Med_Scena0[is.na(FOWARD$Med_Scena0)]<-0
#FOWARD$High_Scen0[is.na(FOWARD$High_Scen0)]<-0
#FOWARD$SOC_t0[is.na(FOWARD$SOC_t0)]<-0
#FOWARD$UNC_t0[is.na(FOWARD$UNC_t0)]<-0
#FOWARD$UNC_BAU[is.na(FOWARD$UNC_BAU)]<-0
#FOWARD$UNC_SSM[is.na(FOWARD$UNC_SSM)]<-0


setwd(WD_MAPS)


# ARKAN: WHY DO SOME OF THESE VALUES TURN INTO ZERO's???!!!!

###################### ATTEMPT AT A SOLUTION FOR ERROR: "Error in update && (!hasValues(y)) : invalid 'x' type in 'x && y'" ##############

#  wopt=list(gdal=c("COMPRESS=NONE", "TFW=YES"), datatype="FL8S"



################# Attempt to create files with a forloop ##############

# attempt <- mask(Country_raster[], FOWARD$SOC_BAU_20


###############################################################################

Country_BAU_2040_Map<-rasterize(FOWARD, Country_raster, field=FOWARD$SOC_BAU_20, update=TRUE)
writeRaster(Country_BAU_2040_Map,filename=paste0(name,"_GSOCseq_finalSOC_BAU_Map030.tiff"),overwrite=TRUE)

# Points to Raster Low Scenario
Country_Lwr_2040_Map<-rasterize(FOWARD, Country_raster ,FOWARD$Low_Scena0, update=TRUE)
writeRaster(Country_Lwr_2040_Map,filename="_GSOCseq_finalSOC_SSM1_Map030.tiff",overwrite=TRUE)

# Points to Raster Med Scenario
Country_Med_2040_Map<-rasterize(FOWARD, Country_raster ,FOWARD$Med_Scena0, update=TRUE)
writeRaster(Country_Med_2040_Map,filename="_GSOCseq_finalSOC_SSM2_Map030.tiff",overwrite=TRUE)

# Points to Raster High Scenario
Country_Hgh_2040_Map<-rasterize(FOWARD, Country_raster ,FOWARD$High_Scen0, update=TRUE)
writeRaster(Country_Hgh_2040_Map,filename="_GSOCseq_finalSOC_SSM3_Map030.tiff",overwrite=TRUE)

# Points to Raster initial SOC (t0) 2018/2020
Country_SOC_2018_Map<-rasterize(FOWARD, Country_raster ,FOWARD$SOC_t0, update=TRUE)
writeRaster(Country_SOC_2018_Map,filename="_GSOCseq_T0_Map030.tiff",overwrite=TRUE)
plot(Country_SOC_2018_Map)

# Difference BAU 2040 - SOC 2018

Diff_BAU_SOC_2018<-Country_BAU_2040_Map-Country_SOC_2018_Map
writeRaster(Diff_BAU_SOC_2018,filename=paste0(name,"_GSOCseq_AbsDiff_BAU_Map030.tiff"),overwrite=TRUE)
writeRaster(Diff_BAU_SOC_2018/20,filename=paste0(name,"_GSOCseq_ASR_BAU_Map030.tiff"),overwrite=TRUE)

# Difference Low Scenario - SOC 2018

Diff_Lw_SOC_2018<-Country_Lwr_2040_Map-Country_SOC_2018_Map
writeRaster(Diff_Lw_SOC_2018,filename=paste0(name,"_GSOCseq_AbsDiff_SSM1_Map030.tiff"),overwrite=TRUE)
writeRaster(Diff_Lw_SOC_2018/20,filename=paste0(name,"_GSOCseq_ASR_SSM1_Map030.tiff"),overwrite=TRUE)

# Difference Med Scenario - SOC 2018

Diff_Md_SOC_2018<-Country_Med_2040_Map-Country_SOC_2018_Map
writeRaster(Diff_Md_SOC_2018,filename=paste0(name,"_GSOCseq_AbsDiff_SSM2_Map030.tiff"),overwrite=TRUE)
writeRaster(Diff_Md_SOC_2018/20,filename=paste0(name,"_GSOCseq_ASR_SSM2_Map030.tiff"),overwrite=TRUE)

# Difference High Scenario - SOC 2018

Diff_Hg_SOC_2018<-Country_Hgh_2040_Map-Country_SOC_2018_Map
writeRaster(Diff_Hg_SOC_2018,filename=paste0(name,"_GSOCseq_AbsDiff_SSM3_Map030"),format="GTiff",overwrite=TRUE)
writeRaster(Diff_Hg_SOC_2018/20,filename=paste0(name,"_GSOCseq_ASR_SSM3_Map030"),format="GTiff",overwrite=TRUE)

# Difference Low Scenario - BAU 2040

Diff_Lw_BAU_2040<-Country_Lwr_2040_Map-Country_BAU_2040_Map
writeRaster(Diff_Lw_BAU_2040,filename=paste0(name,"_GSOCseq_RelDiff_SSM1_Map030"),format="GTiff",overwrite=TRUE)
writeRaster(Diff_Lw_BAU_2040/20,filename=paste0(name,"_GSOCseq_RSR_SSM1_Map030"),format="GTiff",overwrite=TRUE)
plot(Diff_Lw_BAU_2040,  zlim=c(-0.01, 0.01))
# Difference Med Scenario - BAU 2040

Diff_Md_BAU_2040<-Country_Med_2040_Map-Country_BAU_2040_Map
writeRaster(Diff_Md_BAU_2040,filename=paste0(name,"_GSOCseq_RelDiff_SSM2_Map030"),format="GTiff",overwrite=TRUE)
writeRaster(Diff_Md_BAU_2040/20,filename=paste0(name,"_GSOCseq_RSR_SSM2_Map030"),format="GTiff",overwrite=TRUE)

terra::plot(terra::rast(Diff_Md_BAU_2040),
            col=colorRampPalette(c("blue", "cyan", "green"))(4000), 
            axes=FALSE,
            box=FALSE,
            plg=list(cex = 4))
            #plg=list(title="Difference in C per Ha",    # parameters for drawing legend)
             #        title.cex = 2,      # legend title size
              #       cex = 2)) # legend text size


# col=colorRampPalette(c("red", "purple", "blue", "cyan", "green"))(1000)) # zlim=c(-0.0001, 0.0001)
options(scipen=999)
# if error: "Error in plot.new() : figure margins too large


# Difference High Scenario - BAU 2040

Diff_Hg_BAU_2040<-Country_Hgh_2040_Map-Country_BAU_2040_Map
writeRaster(Diff_Hg_BAU_2040,filename=paste0(name,"_GSOCseq_RelDiff_SSM3_Map030"),format="GTiff",overwrite=TRUE)
writeRaster(Diff_Hg_BAU_2040/20,filename=paste0(name,"_GSOCseq_RSR_SSM3_Map030"),format="GTiff",overwrite=TRUE)
plot(Diff_Hg_BAU_2040,  zlim=c(-0.01, 0.01))
# Uncertainties SOC 2018

UNC_2018<-rasterize(FOWARD, Country_raster ,FOWARD$UNC_t0, update=TRUE)
writeRaster(UNC_2018,filename=paste0(name,"_GSOCseq_T0_UncertaintyMap030.tiff"),overwrite=TRUE)
UNC_2018
# Uncertainties SOC BAU 2038

UNC_BAU<-rasterize(FOWARD, Country_raster ,FOWARD$UNC_BAU, update=TRUE)
writeRaster(UNC_BAU,filename=paste0(name,"_GSOCseq_BAU_UncertaintyMap030"),format="GTiff",overwrite=TRUE)

# Uncertainties SOC SSM 

UNC_SSM<-rasterize(FOWARD, Country_raster ,FOWARD$UNC_SSM, update=TRUE)
writeRaster(UNC_SSM,filename=paste0(name,"_GSOCseq_SSM_UncertaintyMap030"),format="GTiff",overwrite=TRUE)

# Uncertainties for the Absolute difference SSM_ - SOC2018


UNC_abs_rate_BAU<- sqrt((FOWARD$UNC_BAU*FOWARD$SOC_BAU_20)^2 + ((FOWARD$UNC_t0*FOWARD$SOC_t0)^2)/abs(FOWARD$SOC_t0+FOWARD$SOC_BAU_20))
UNC_abs_rate_BAU[is.na(UNC_abs_rate_BAU)]<-0
UNC_abs_rate_BAU_Map<-rasterize(FOWARD, Country_raster ,UNC_abs_rate_BAU, update=TRUE)
writeRaster(UNC_abs_rate_BAU_Map,filename=paste0(name,"_GSOCseq_ASR_BAU_UncertaintyMap030"),format="GTiff",overwrite=TRUE)

UNC_abs_rate_Lw<- sqrt((FOWARD$UNC_SSM*FOWARD$Low_Scena0)^2 + (FOWARD$UNC_t0*FOWARD$SOC_t0)^2)/abs(FOWARD$SOC_t0+FOWARD$Low_Scena0)
UNC_abs_rate_Lw[is.na(UNC_abs_rate_Lw)] <- 0
UNC_abs_rate_Lw_Map<-rasterize(FOWARD, Country_raster ,UNC_abs_rate_Lw, update=TRUE)
writeRaster(UNC_abs_rate_Lw_Map,filename=paste0(name,"_GSOCseq_ASR_SSM1_UncertaintyMap030"),format="GTiff",overwrite=TRUE)

UNC_abs_rate_Md<- sqrt((FOWARD$UNC_SSM*FOWARD$Med_Scena0)^2 + (FOWARD$UNC_t0*FOWARD$SOC_t0)^2)/abs(FOWARD$SOC_t0+FOWARD$Med_Scena0)
UNC_abs_rate_Md[is.na(UNC_abs_rate_Md)]<-0
UNC_abs_rate_Md_Map<-rasterize(FOWARD, Country_raster ,UNC_abs_rate_Md, update=TRUE)
writeRaster(UNC_abs_rate_Md_Map,filename=paste0(name,"_GSOCseq_ASR_SSM2_UncertaintyMap030"),format="GTiff",overwrite=TRUE)

UNC_abs_rate_Hg<- sqrt((FOWARD$UNC_SSM*FOWARD$High_Scen0)^2 + (FOWARD$UNC_t0*FOWARD$SOC_t0)^2)/abs(FOWARD$SOC_t0+FOWARD$High_Scen0)
UNC_abs_rate_Hg[is.na(UNC_abs_rate_Hg)]<-0
UNC_abs_rate_Hg_Map<-rasterize(FOWARD, Country_raster ,UNC_abs_rate_Hg, update=TRUE)
writeRaster(UNC_abs_rate_Hg_Map,filename=paste0(name,"_GSOCseq_ASR_SSM3_UncertaintyMap030"),format="GTiff",overwrite=TRUE)

# Uncertainties for the Relative difference  SSM_ - SOCBAU

UNC_Rel_rate_Lw<- sqrt((FOWARD$UNC_SSM*FOWARD$Low_Scena0)^2 + (FOWARD$UNC_BAU*FOWARD$SOC_BAU_20)^2)/abs(FOWARD$SOC_BAU_20+FOWARD$Low_Scena0)
UNC_Rel_rate_Lw[is.na(UNC_Rel_rate_Lw)]<-0
UNC_Rel_rate_Lw_Map<-rasterize(FOWARD, Country_raster ,UNC_Rel_rate_Lw, update=TRUE)
writeRaster(UNC_Rel_rate_Lw_Map,filename=paste0(name,"_GSOCseq_RSR_SSM1_UncertaintyMap030"),format="GTiff",overwrite=TRUE)

UNC_Rel_rate_Md<- sqrt((FOWARD$UNC_SSM*FOWARD$Med_Scena0)^2 + (FOWARD$UNC_BAU*FOWARD$SOC_BAU_20)^2)/abs(FOWARD$SOC_BAU_20+FOWARD$Med_Scen_0)
UNC_Rel_rate_Md[is.na(UNC_Rel_rate_Md)]<-0
UNC_Rel_rate_Md_Map<-rasterize(FOWARD, Country_raster ,UNC_Rel_rate_Md, update=TRUE)
writeRaster(UNC_Rel_rate_Md_Map,filename=paste0(name,"_GSOCseq_RSR_SSM2_UncertaintyMap030"),format="GTiff",overwrite=TRUE)

UNC_Rel_rate_Hg<- sqrt((FOWARD$UNC_SSM*FOWARD$High_Scen0)^2 + (FOWARD$UNC_BAU*FOWARD$SOC_BAU_20)^2)/abs(FOWARD$SOC_BAU_20+FOWARD$High_Scen0)
UNC_Rel_rate_Hg[is.na(UNC_Rel_rate_Hg)]<-0
UNC_Rel_rate_Hg_Map<-rasterize(FOWARD, Country_raster ,UNC_Rel_rate_Hg, update=TRUE)
writeRaster(UNC_Rel_rate_Hg_Map,filename=paste0(name,"_GSOCseq_RSR_SSM3_UncertaintyMap030"),format="GTiff",overwrite=TRUE)


