library(terra)
library(SoilR)
library(soilassessment)

rm(list=ls()) 

WD_OUT<-("H:/BisonRothC/OUTPUTS/3_FORWARD")
working_dir<-setwd("H:/BisonRothC")
# OPEN THE VECTOR OF POINTS
Vector<-vect("INPUTS/TARGET_POINTS/target_points.shp")
# OPEN THE RESULT VECTOR FROM THE WARM UP PROCESS
WARM_UP<-vect("OUTPUTS/2_WARM_UP/WARM_UP_AOI_final/WARM_UP_AOI_final.shp")
# OPEN THE STACK WITH THE VARIABLES FOR THE FOWARD PROCESS
Stack_Set_1<- rast("INPUTS/STACK/Stack_Set_FOWARD1.tif")




# Set the increase in Carbon input for each land use and each scenario
#Crops and Crop trees
Low_Crops<-0.9
Med_Crops<-1
High_Crops<-1.1
#Shrublands, Grasslands , Herbaceous vegetation flooded & Sparse Vegetation
Low_Grass<-0.9
Med_Grass<-1
High_Grass<-1.1
#Forest
Low_Forest<-0.9
Med_Forest<-1
High_Forest<-1.1
#Paddy Fields
Low_PaddyFields<- 0.9
Med_PaddyFields<-1
High_PaddyFields<-1.1

# extract variables to points
Variables<-extract(Stack_Set_1,Vector,bind=TRUE)

# Creates an empty vector
FOWARD<-Vector
check_Points <- Vector

# Extract the layers from the Vector
SOC_im<-WARM_UP[[4]]
clay_im<-Variables[[3]] 
Cinputs_im<-WARM_UP[[10]] 
DR_im<-Variables[[40]]
LU_im<-Variables[[41]]

# Define the years to run the model
years=seq(1/12,100,by=1/12)


##################### NPP FOR BISON ##############

NPP<-rast(c("H:/BisonRothC/INPUTS/NPP/NPP_MIAMI_MEAN_81-00_AOI.tif"))
NPP_MEAN_MIN<-rast(c("H:/BisonRothC/INPUTS/NPP/NPP_MIAMI_MEAN_81-00_AOI_MIN.tif"))
NPP_MEAN_MAX<-rast(c("H:/BisonRothC/INPUTS/NPP/NPP_MIAMI_MEAN_81-00_AOI_MAX.tif"))


Variables<- extract(NPP, Variables, bind=TRUE) #NPP_mean  = Variables[i, 54]
Variables<- extract(NPP_MEAN_MIN, Variables, bind=TRUE)
Variables<- extract(NPP_MEAN_MAX, Variables, bind=TRUE)

#################################################



#############function set up starts###############
Roth_C<-function(Cinputs,years,DPMptf, RPMptf, BIOptf, HUMptf, FallIOM,Temp,Precip,Evp,Cov,Cov1,Cov2,soil.thick,SOC,clay,DR,bare1,LU)
{
  # Paddy Fields coefficent fPR = 0.4 if the target point is class = 13 , else fPR=1
  # From Shirato and Yukozawa 2004
  
  fPR=(LU == 13)*0.4 + (LU!=13)*1
  
  #Temperature effects per month
  fT=fT.RothC(Temp[,2]) 
  
  #Moisture effects per month . 
  
  fw1func<-function(P, E, S.Thick = 30, pClay = 32.0213, pE = 1, bare) 
  {
    
    M = P - E * pE
    Acc.TSMD = NULL
    for (i in 2:length(M)) {
      B = ifelse(bare[i] == FALSE, 1, 1.8)
      Max.TSMD = -(20 + 1.3 * pClay - 0.01 * (pClay^2)) * (S.Thick/23) * (1/B)
      Acc.TSMD[1] = ifelse(M[1] > 0, 0, M[1])
      if (Acc.TSMD[i - 1] + M[i] < 0) {
        Acc.TSMD[i] = Acc.TSMD[i - 1] + M[i]
      }
      else (Acc.TSMD[i] = 0)
      if (Acc.TSMD[i] <= Max.TSMD) {
        Acc.TSMD[i] = Max.TSMD
      }
    }
    b = ifelse(Acc.TSMD > 0.444 * Max.TSMD, 1, (0.2 + 0.8 * ((Max.TSMD - 
                                                                Acc.TSMD)/(Max.TSMD - 0.444 * Max.TSMD))))
    b<-clamp(b,lower=0.2)
    return(data.frame(b))
  }
  
  fW_2<- fw1func(P=(Precip[,2]), E=(Evp[,2]), S.Thick = soil.thick, pClay = clay, pE = 1, bare=bare1)$b 
  
  #Vegetation Cover effects  
  
  fC<-Cov2[,2]
  
  # Set the factors frame for Model calculations
  
  xi.frame=data.frame(years,rep(fT*fW_2*fC*fPR,length.out=length(years)))
  
  # RUN THE MODEL from SoilR
  #Loads the model 
  #Model3_spin=RothCModel(t=years,C0=c(DPMptf[[1]], RPMptf[[1]], BIOptf[[1]], HUMptf[[1]], FallIOM[[1]]),In=Cinputs,DR=DR,clay=clay,xi=xi.frame, pass=TRUE) 
  #Ct3_spin=getC(Model3_spin)
  
  # RUN THE MODEL from soilassesment
  # Arkan: I've made a mess with the c() or listing of the variables, hence the weird indexing, or maybe it's supposed to be like this. either way, it currently works
  Model3_spin=carbonTurnover(tt=years,C0=c(DPMptf[[1]][[1]], RPMptf[[1]][[1]], BIOptf[[1]][[1]], HUMptf[[1]][[1]], FallIOM[[1]][[1]]),In=Cinputs,Dr=DR,clay=clay,effcts=xi.frame, "euler") 
  
  Ct3_spin=Model3_spin[,2:6]
  
  # Get the final pools of the time series
  
  poolSize3_spin=as.numeric(tail(Ct3_spin,1))
  
  return(poolSize3_spin)
}
################function set up ends#############


############################################
#########YEAR SAMPLES######################
check_Years_10 <- years[1:(12*10)]
check_Years_20 <- years[1:(12*20)]
check_Years_30 <- years[1:(12*30)]
check_Years_40 <- years[1:(12*40)]
check_Years_50 <- years[1:(12*50)]
check_Years_60 <- years[1:(12*60)]
check_Years_70 <- years[1:(12*70)]
check_Years_80 <- years[1:(12*80)]
check_Years_90 <- years[1:(12*90)]
check_Years_100 <- years[1:(12*100)]




#############################################










# Iterates over the area of interest
##################for loop starts###############
for (i in 1:dim(Variables)[1]) {
  
  # Extract the variables 
  
  Vect<-as.data.frame(Variables[i,])
  
  Temp<-as.data.frame(t(Vect[4:15]))
  Temp<-data.frame(Month=1:12, Temp=Temp[,1])
  
  Precip<-as.data.frame(t(Vect[16:27]))
  Precip<-data.frame(Month=1:12, Precip=Precip[,1])
  
  Evp<-as.data.frame(t(Vect[28:39]))
  Evp<-data.frame(Month=1:12, Evp=Evp[,1])
  
  Cov<-as.data.frame(t(Vect[42:53]))
  Cov1<-data.frame(Cov=Cov[,1])
  Cov2<-data.frame(Month=1:12, Cov=Cov[,1])
  
  
  
  
  
  #Avoid calculus over Na values 
  
  if (any(is.na(Evp[,2])) | any(is.na(Temp[,2])) | any(is.na(SOC_im[i, ])) | any(is.na(clay_im[i, ])) | any(is.na(Precip[,2]))  |  any(is.na(Cov2[,2]))  |  any(is.na(Cov1[,1])) | any(is.na(Variables[i, 54])) | any(is.na(Cinputs_im[i, ])) | any(is.na(DR_im[i, ])) | (Cinputs_im[i, ]<0) |  (SOC_im[i, ]<0) | (clay_im[i, ]<0) ) {FOWARD[i,2]<-0}else{
    
    
    # Set the variables from the images
    
    soil.thick=30  #Soil thickness (organic layer topsoil), in cm
    SOC<-SOC_im[i, ]      #Soil organic carbon in Mg/ha 
    clay<-clay_im[i, ]        #Percent clay %
    Cinputs<-Cinputs_im[i, ]    #Annual C inputs to soil in Mg/ha/yr
    
    DR<-DR_im[i, ]              # DPM/RPM (decomplosable vs resistant plant material.)
    bare1<-(Cov1>0.8)           # If the surface is bare or vegetated
    LU<-LU_im[i, ]
    
    
    ######## BISON contributions to DR ##############
    
    # important to include NPP mean and max because that will effect how much of a contribution bison poop makes.
    
    b_R_1 <- 0.5 # bison waste estimated decomposability (bison ratio)
    b_R_2 <- 0.41
    b_R_3 <- 0.33
    
    
    
    B_c_mean <- (0.0687/Variables[1, 54][[1]][[1]]) # % of NPP that goes through bison (bison contribution)
    DR_mean_b1 <- B_c_mean*b_R_2 + DR*(1-B_c_mean) # o.o687 t ha-1 year-1, 0.5 DPM/RPM ratio
    
    
    B_c_mean <- (0.0687/Variables[1, 54][[1]][[1]]) # % of NPP that goes through bison (bison contribution)
    DR_mean_b2 <- B_c_mean*b_R_2 + DR*(1-B_c_mean)
    
    B_c_mean <- (0.0687/Variables[1, 54][[1]][[1]]) # % of NPP that goes through bison (bison contribution)
    DR_mean_b3 <- B_c_mean*b_R_3 + DR*(1-B_c_mean)
    
    #B_c_min <- (0.0687/NPP_M_MIN[w]) 
    #DR_min <- B_c_min*b_R + DR*(1-B_c_min)
                                                          ##### these two sections have been commented out 
    #B_c_max <- (0.0687/NPP_M_MAX[w]) 
    #DR_max <- B_c_max*b_R + DR*(1-B_c_max)
    
    
    
    
    
    
    
    
    
    # Final calculation of SOC  20 years in the future  (Business as usual)
    
    f_bau<-Roth_C(Cinputs=Cinputs,years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR,bare1=bare1,LU=LU)
    f_bau_t<-f_bau[1]+f_bau[2]+f_bau[3]+f_bau[4]+f_bau[5]
    
    #Unc BAU minimum 
    Cinputs_min<-WARM_UP[i,23][[1]][[1]]
    Cinputs_max<-WARM_UP[i,24][[1]][[1]]
    SOC_t0_min<-WARM_UP[i,11][[1]][[1]]
    SOC_t0_max<-WARM_UP[i,17][[1]][[1]]
    
    f_bau_min<-Roth_C(Cinputs=Cinputs_min,years=years,DPMptf=WARM_UP[i,12], RPMptf=WARM_UP[i,13], BIOptf=WARM_UP[i,14], HUMptf=WARM_UP[i,15], FallIOM=WARM_UP[i,16],Temp=Temp*1.02,Precip=Precip*0.95,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*0.8,clay=clay*0.9,DR=DR,bare1=bare1,LU=LU)
    f_bau_t_min<-f_bau_min[1]+f_bau_min[2]+f_bau_min[3]+f_bau_min[4]+f_bau_min[5]
    
    #Unc BAU maximum
    
    f_bau_max<-Roth_C(Cinputs=Cinputs_max,years=years,DPMptf=WARM_UP[i,18], RPMptf=WARM_UP[i,19], BIOptf=WARM_UP[i,20], HUMptf=WARM_UP[i,21], FallIOM=WARM_UP[i,22],Temp=Temp*0.98,Precip=Precip*1.05,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*1.2,clay=clay*1.1,DR=DR,bare1=bare1,LU=LU)
    f_bau_t_max<-f_bau_max[1]+f_bau_max[2]+f_bau_max[3]+f_bau_max[4]+f_bau_max[5]
      
      
      # Arkan's forest implementation:
      #Forest cover
    if (TRUE) {
      f_low<-Roth_C(Cinputs=(Cinputs*Low_Forest),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b1,bare1=bare1,LU=LU)
      f_low_t<-f_low[1]+f_low[2]+f_low[3]+f_low[4]+f_low[5]
        
      f_med<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
      f_med_t<-f_med[1]+f_med[2]+f_med[3]+f_med[4]+f_med[5]
        
      f_high<-Roth_C(Cinputs=(Cinputs*High_Forest),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b3,bare1=bare1,LU=LU)
      f_high_t<-f_high[1]+f_high[2]+f_high[3]+f_high[4]+f_high[5]
        
      #SSM forest unc min
        
      f_med_min<-Roth_C(Cinputs=(Cinputs_min*(Low_Forest-0.15)),years=years,DPMptf=WARM_UP[i,12], RPMptf=WARM_UP[i,13], BIOptf=WARM_UP[i,14], HUMptf=WARM_UP[i,15], FallIOM=WARM_UP[i,16],Temp=Temp*1.02,Precip=Precip*0.95,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*0.8,clay=clay*0.9,DR=DR_mean_b2,bare1=bare1,LU=LU)
      f_med_t_min<-f_med_min[1]+f_med_min[2]+f_med_min[3]+f_med_min[4]+f_med_min[5]
        
      #SSM forest unc max
        
      f_med_max<-Roth_C(Cinputs=(Cinputs_max*(Med_Forest+0.15)),years=years,DPMptf=WARM_UP[i,18], RPMptf=WARM_UP[i,19], BIOptf=WARM_UP[i,20], HUMptf=WARM_UP[i,21], FallIOM=WARM_UP[i,22],Temp=Temp*0.98,Precip=Precip*1.05,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*1.2,clay=clay*1.1,DR=DR_mean_b2,bare1=bare1,LU=LU)
      f_med_t_max<-f_med_max[1]+f_med_max[2]+f_med_max[3]+f_med_max[4]+f_med_max[5]
      
    }
    
    else{
      f_bau_t<-0
      f_low_t<-0
      f_med_t<-0
      f_high_t<-0
      f_bau_t_min<-0
      f_bau_t_max<-0
      f_med_t_min<-0
      f_med_t_max<-0
      SOC_t0_min<-0
      SOC_t0_max<-0
      
    }
    
    
    
    #############ARKANS EVIL VISUALIZATION PLAAAAAAAAAAAAN#################################################################################
   
    
    ########################## Year samples ##################
    f_bau_check_10 <- Roth_C(Cinputs=Cinputs,years=check_Years_10,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR,bare1=bare1,LU=LU)
    f_bau_check_10_t <- f_bau_check_10[1] + f_bau_check_10[2] + f_bau_check_10[3] + f_bau_check_10[4] + f_bau_check_10[5]
    f_bau_check_20 <- Roth_C(Cinputs=Cinputs,years=check_Years_20,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR,bare1=bare1,LU=LU)
    f_bau_check_20_t <- f_bau_check_20[1] + f_bau_check_20[2] + f_bau_check_20[3] + f_bau_check_20[4] + f_bau_check_20[5]
    f_bau_check_30 <- Roth_C(Cinputs=Cinputs,years=check_Years_30,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR,bare1=bare1,LU=LU)
    f_bau_check_30_t <- f_bau_check_30[1] + f_bau_check_30[2] + f_bau_check_30[3] + f_bau_check_30[4] + f_bau_check_30[5]
    f_bau_check_40 <- Roth_C(Cinputs=Cinputs,years=check_Years_40,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR,bare1=bare1,LU=LU)
    f_bau_check_40_t <- f_bau_check_40[1] + f_bau_check_40[2] + f_bau_check_40[3] + f_bau_check_40[4] + f_bau_check_40[5]
    f_bau_check_50 <- Roth_C(Cinputs=Cinputs,years=check_Years_50,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR,bare1=bare1,LU=LU)
    f_bau_check_50_t <- f_bau_check_50[1] + f_bau_check_50[2] + f_bau_check_50[3] + f_bau_check_50[4] + f_bau_check_50[5]
    f_bau_check_60 <- Roth_C(Cinputs=Cinputs,years=check_Years_60,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR,bare1=bare1,LU=LU)
    f_bau_check_60_t <- f_bau_check_60[1] + f_bau_check_60[2] + f_bau_check_60[3] + f_bau_check_60[4] + f_bau_check_60[5]
    f_bau_check_70 <- Roth_C(Cinputs=Cinputs,years=check_Years_70,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR,bare1=bare1,LU=LU)
    f_bau_check_70_t <- f_bau_check_70[1] + f_bau_check_70[2] + f_bau_check_70[3] + f_bau_check_70[4] + f_bau_check_70[5]
    f_bau_check_80 <- Roth_C(Cinputs=Cinputs,years=check_Years_80,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR,bare1=bare1,LU=LU)
    f_bau_check_80_t <- f_bau_check_80[1] + f_bau_check_80[2] + f_bau_check_80[3] + f_bau_check_80[4] + f_bau_check_80[5]
    f_bau_check_90 <- Roth_C(Cinputs=Cinputs,years=check_Years_90,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR,bare1=bare1,LU=LU)
    f_bau_check_90_t <- f_bau_check_90[1] + f_bau_check_90[2] + f_bau_check_90[3] + f_bau_check_90[4] + f_bau_check_90[5]
    f_bau_check_100 <- Roth_C(Cinputs=Cinputs,years=check_Years_100,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR,bare1=bare1,LU=LU)
    f_bau_check_100_t <- f_bau_check_100[1] + f_bau_check_100[2] + f_bau_check_100[3] + f_bau_check_100[4] + f_bau_check_100[5]
    
    f_bau_min_check_10 <-Roth_C(Cinputs=Cinputs_min,years=check_Years_10,DPMptf=WARM_UP[i,12], RPMptf=WARM_UP[i,13], BIOptf=WARM_UP[i,14], HUMptf=WARM_UP[i,15], FallIOM=WARM_UP[i,16],Temp=Temp*1.02,Precip=Precip*0.95,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*0.8,clay=clay*0.9,DR=DR,bare1=bare1,LU=LU)
    f_bau_min_check_10_t <- f_bau_min_check_10[1] + f_bau_min_check_10[2] + f_bau_min_check_10[3] + f_bau_min_check_10[4] + f_bau_min_check_10[5]
    f_bau_min_check_20 <-Roth_C(Cinputs=Cinputs_min,years=check_Years_20,DPMptf=WARM_UP[i,12], RPMptf=WARM_UP[i,13], BIOptf=WARM_UP[i,14], HUMptf=WARM_UP[i,15], FallIOM=WARM_UP[i,16],Temp=Temp*1.02,Precip=Precip*0.95,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*0.8,clay=clay*0.9,DR=DR,bare1=bare1,LU=LU)
    f_bau_min_check_20_t <- f_bau_min_check_20[1] + f_bau_min_check_20[2] + f_bau_min_check_20[3] + f_bau_min_check_20[4] + f_bau_min_check_20[5]
    f_bau_min_check_30 <-Roth_C(Cinputs=Cinputs_min,years=check_Years_30,DPMptf=WARM_UP[i,12], RPMptf=WARM_UP[i,13], BIOptf=WARM_UP[i,14], HUMptf=WARM_UP[i,15], FallIOM=WARM_UP[i,16],Temp=Temp*1.02,Precip=Precip*0.95,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*0.8,clay=clay*0.9,DR=DR,bare1=bare1,LU=LU)
    f_bau_min_check_30_t <- f_bau_min_check_30[1] + f_bau_min_check_30[2] + f_bau_min_check_30[3] + f_bau_min_check_30[4] + f_bau_min_check_30[5]
    f_bau_min_check_40 <-Roth_C(Cinputs=Cinputs_min,years=check_Years_40,DPMptf=WARM_UP[i,12], RPMptf=WARM_UP[i,13], BIOptf=WARM_UP[i,14], HUMptf=WARM_UP[i,15], FallIOM=WARM_UP[i,16],Temp=Temp*1.02,Precip=Precip*0.95,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*0.8,clay=clay*0.9,DR=DR,bare1=bare1,LU=LU)
    f_bau_min_check_40_t <- f_bau_min_check_40[1] + f_bau_min_check_40[2] + f_bau_min_check_40[3] + f_bau_min_check_40[4] + f_bau_min_check_40[5]
    f_bau_min_check_50 <-Roth_C(Cinputs=Cinputs_min,years=check_Years_50,DPMptf=WARM_UP[i,12], RPMptf=WARM_UP[i,13], BIOptf=WARM_UP[i,14], HUMptf=WARM_UP[i,15], FallIOM=WARM_UP[i,16],Temp=Temp*1.02,Precip=Precip*0.95,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*0.8,clay=clay*0.9,DR=DR,bare1=bare1,LU=LU)
    f_bau_min_check_50_t <- f_bau_min_check_50[1] + f_bau_min_check_50[2] + f_bau_min_check_50[3] + f_bau_min_check_50[4] + f_bau_min_check_50[5]
    f_bau_min_check_60 <-Roth_C(Cinputs=Cinputs_min,years=check_Years_60,DPMptf=WARM_UP[i,12], RPMptf=WARM_UP[i,13], BIOptf=WARM_UP[i,14], HUMptf=WARM_UP[i,15], FallIOM=WARM_UP[i,16],Temp=Temp*1.02,Precip=Precip*0.95,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*0.8,clay=clay*0.9,DR=DR,bare1=bare1,LU=LU)
    f_bau_min_check_60_t <- f_bau_min_check_60[1] + f_bau_min_check_60[2] + f_bau_min_check_60[3] + f_bau_min_check_60[4] + f_bau_min_check_60[5]
    f_bau_min_check_70 <-Roth_C(Cinputs=Cinputs_min,years=check_Years_70,DPMptf=WARM_UP[i,12], RPMptf=WARM_UP[i,13], BIOptf=WARM_UP[i,14], HUMptf=WARM_UP[i,15], FallIOM=WARM_UP[i,16],Temp=Temp*1.02,Precip=Precip*0.95,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*0.8,clay=clay*0.9,DR=DR,bare1=bare1,LU=LU)
    f_bau_min_check_70_t <- f_bau_min_check_70[1] + f_bau_min_check_70[2] + f_bau_min_check_70[3] + f_bau_min_check_70[4] + f_bau_min_check_70[5]
    f_bau_min_check_80 <-Roth_C(Cinputs=Cinputs_min,years=check_Years_80,DPMptf=WARM_UP[i,12], RPMptf=WARM_UP[i,13], BIOptf=WARM_UP[i,14], HUMptf=WARM_UP[i,15], FallIOM=WARM_UP[i,16],Temp=Temp*1.02,Precip=Precip*0.95,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*0.8,clay=clay*0.9,DR=DR,bare1=bare1,LU=LU)
    f_bau_min_check_80_t <- f_bau_min_check_80[1] + f_bau_min_check_80[2] + f_bau_min_check_80[3] + f_bau_min_check_80[4] + f_bau_min_check_80[5]
    f_bau_min_check_90 <-Roth_C(Cinputs=Cinputs_min,years=check_Years_90,DPMptf=WARM_UP[i,12], RPMptf=WARM_UP[i,13], BIOptf=WARM_UP[i,14], HUMptf=WARM_UP[i,15], FallIOM=WARM_UP[i,16],Temp=Temp*1.02,Precip=Precip*0.95,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*0.8,clay=clay*0.9,DR=DR,bare1=bare1,LU=LU)
    f_bau_min_check_90_t <- f_bau_min_check_90[1] + f_bau_min_check_90[2] + f_bau_min_check_90[3] + f_bau_min_check_90[4] + f_bau_min_check_90[5]
    f_bau_min_check_100 <-Roth_C(Cinputs=Cinputs_min,years=check_Years_100,DPMptf=WARM_UP[i,12], RPMptf=WARM_UP[i,13], BIOptf=WARM_UP[i,14], HUMptf=WARM_UP[i,15], FallIOM=WARM_UP[i,16],Temp=Temp*1.02,Precip=Precip*0.95,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*0.8,clay=clay*0.9,DR=DR,bare1=bare1,LU=LU)
    f_bau_min_check_100_t <- f_bau_min_check_100[1] + f_bau_min_check_100[2] + f_bau_min_check_100[3] + f_bau_min_check_100[4] + f_bau_min_check_100[5]
    
    f_bau_max_check_10<-Roth_C(Cinputs=Cinputs_max,years=check_Years_10,DPMptf=WARM_UP[i,18], RPMptf=WARM_UP[i,19], BIOptf=WARM_UP[i,20], HUMptf=WARM_UP[i,21], FallIOM=WARM_UP[i,22],Temp=Temp*0.98,Precip=Precip*1.05,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*1.2,clay=clay*1.1,DR=DR,bare1=bare1,LU=LU)
    f_bau_max_check_10_t <- f_bau_max_check_10[1] + f_bau_max_check_10[2] + f_bau_max_check_10[3] + f_bau_max_check_10[4] + f_bau_max_check_10[5]
    f_bau_max_check_20<-Roth_C(Cinputs=Cinputs_max,years=check_Years_20,DPMptf=WARM_UP[i,18], RPMptf=WARM_UP[i,19], BIOptf=WARM_UP[i,20], HUMptf=WARM_UP[i,21], FallIOM=WARM_UP[i,22],Temp=Temp*0.98,Precip=Precip*1.05,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*1.2,clay=clay*1.1,DR=DR,bare1=bare1,LU=LU)
    f_bau_max_check_20_t <- f_bau_max_check_20[1] + f_bau_max_check_20[2] + f_bau_max_check_20[3] + f_bau_max_check_20[4] + f_bau_max_check_20[5]
    f_bau_max_check_30<-Roth_C(Cinputs=Cinputs_max,years=check_Years_30,DPMptf=WARM_UP[i,18], RPMptf=WARM_UP[i,19], BIOptf=WARM_UP[i,20], HUMptf=WARM_UP[i,21], FallIOM=WARM_UP[i,22],Temp=Temp*0.98,Precip=Precip*1.05,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*1.2,clay=clay*1.1,DR=DR,bare1=bare1,LU=LU)
    f_bau_max_check_30_t <- f_bau_max_check_30[1] + f_bau_max_check_30[2] + f_bau_max_check_30[3] + f_bau_max_check_30[4] + f_bau_max_check_30[5]
    f_bau_max_check_40<-Roth_C(Cinputs=Cinputs_max,years=check_Years_40,DPMptf=WARM_UP[i,18], RPMptf=WARM_UP[i,19], BIOptf=WARM_UP[i,20], HUMptf=WARM_UP[i,21], FallIOM=WARM_UP[i,22],Temp=Temp*0.98,Precip=Precip*1.05,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*1.2,clay=clay*1.1,DR=DR,bare1=bare1,LU=LU)
    f_bau_max_check_40_t <- f_bau_max_check_40[1] + f_bau_max_check_40[2] + f_bau_max_check_40[3] + f_bau_max_check_40[4] + f_bau_max_check_40[5]
    f_bau_max_check_50<-Roth_C(Cinputs=Cinputs_max,years=check_Years_50,DPMptf=WARM_UP[i,18], RPMptf=WARM_UP[i,19], BIOptf=WARM_UP[i,20], HUMptf=WARM_UP[i,21], FallIOM=WARM_UP[i,22],Temp=Temp*0.98,Precip=Precip*1.05,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*1.2,clay=clay*1.1,DR=DR,bare1=bare1,LU=LU)
    f_bau_max_check_50_t <- f_bau_max_check_50[1] + f_bau_max_check_50[2] + f_bau_max_check_50[3] + f_bau_max_check_50[4] + f_bau_max_check_50[5]
    f_bau_max_check_60<-Roth_C(Cinputs=Cinputs_max,years=check_Years_60,DPMptf=WARM_UP[i,18], RPMptf=WARM_UP[i,19], BIOptf=WARM_UP[i,20], HUMptf=WARM_UP[i,21], FallIOM=WARM_UP[i,22],Temp=Temp*0.98,Precip=Precip*1.05,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*1.2,clay=clay*1.1,DR=DR,bare1=bare1,LU=LU)
    f_bau_max_check_60_t <- f_bau_max_check_60[1] + f_bau_max_check_60[2] + f_bau_max_check_60[3] + f_bau_max_check_60[4] + f_bau_max_check_60[5]
    f_bau_max_check_70<-Roth_C(Cinputs=Cinputs_max,years=check_Years_70,DPMptf=WARM_UP[i,18], RPMptf=WARM_UP[i,19], BIOptf=WARM_UP[i,20], HUMptf=WARM_UP[i,21], FallIOM=WARM_UP[i,22],Temp=Temp*0.98,Precip=Precip*1.05,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*1.2,clay=clay*1.1,DR=DR,bare1=bare1,LU=LU)
    f_bau_max_check_70_t <- f_bau_max_check_70[1] + f_bau_max_check_70[2] + f_bau_max_check_70[3] + f_bau_max_check_70[4] + f_bau_max_check_70[5]
    f_bau_max_check_80<-Roth_C(Cinputs=Cinputs_max,years=check_Years_80,DPMptf=WARM_UP[i,18], RPMptf=WARM_UP[i,19], BIOptf=WARM_UP[i,20], HUMptf=WARM_UP[i,21], FallIOM=WARM_UP[i,22],Temp=Temp*0.98,Precip=Precip*1.05,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*1.2,clay=clay*1.1,DR=DR,bare1=bare1,LU=LU)
    f_bau_max_check_80_t <- f_bau_max_check_80[1] + f_bau_max_check_80[2] + f_bau_max_check_80[3] + f_bau_max_check_80[4] + f_bau_max_check_80[5]
    f_bau_max_check_90<-Roth_C(Cinputs=Cinputs_max,years=check_Years_90,DPMptf=WARM_UP[i,18], RPMptf=WARM_UP[i,19], BIOptf=WARM_UP[i,20], HUMptf=WARM_UP[i,21], FallIOM=WARM_UP[i,22],Temp=Temp*0.98,Precip=Precip*1.05,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*1.2,clay=clay*1.1,DR=DR,bare1=bare1,LU=LU)
    f_bau_max_check_90_t <- f_bau_max_check_90[1] + f_bau_max_check_90[2] + f_bau_max_check_90[3] + f_bau_max_check_90[4] + f_bau_max_check_90[5]
    f_bau_max_check_100<-Roth_C(Cinputs=Cinputs_max,years=check_Years_100,DPMptf=WARM_UP[i,18], RPMptf=WARM_UP[i,19], BIOptf=WARM_UP[i,20], HUMptf=WARM_UP[i,21], FallIOM=WARM_UP[i,22],Temp=Temp*0.98,Precip=Precip*1.05,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*1.2,clay=clay*1.1,DR=DR,bare1=bare1,LU=LU)
    f_bau_max_check_100_t <- f_bau_max_check_100[1] + f_bau_max_check_100[2] + f_bau_max_check_100[3] + f_bau_max_check_100[4] + f_bau_max_check_100[5]
    
    
    
    
    
    F_SSM_low_check_10<-Roth_C(Cinputs=(Cinputs*Low_Forest),years=check_Years_10,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b1,bare1=bare1,LU=LU)
    F_SSM_low_check_10_t <- F_SSM_low_check_10[1] + F_SSM_low_check_10[2] + F_SSM_low_check_10[3] + F_SSM_low_check_10[4] + F_SSM_low_check_10[5]
    F_SSM_low_check_20<-Roth_C(Cinputs=(Cinputs*Low_Forest),years=check_Years_20,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b1,bare1=bare1,LU=LU)
    F_SSM_low_check_20_t <- F_SSM_low_check_20[1] + F_SSM_low_check_20[2] + F_SSM_low_check_20[3] + F_SSM_low_check_20[4] + F_SSM_low_check_20[5]
    F_SSM_low_check_30<-Roth_C(Cinputs=(Cinputs*Low_Forest),years=check_Years_30,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b1,bare1=bare1,LU=LU)
    F_SSM_low_check_30_t <- F_SSM_low_check_30[1] + F_SSM_low_check_30[2] + F_SSM_low_check_30[3] + F_SSM_low_check_30[4] + F_SSM_low_check_30[5]
    F_SSM_low_check_40<-Roth_C(Cinputs=(Cinputs*Low_Forest),years=check_Years_40,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b1,bare1=bare1,LU=LU)
    F_SSM_low_check_40_t <- F_SSM_low_check_40[1] + F_SSM_low_check_40[2] + F_SSM_low_check_40[3] + F_SSM_low_check_40[4] + F_SSM_low_check_40[5]
    F_SSM_low_check_50<-Roth_C(Cinputs=(Cinputs*Low_Forest),years=check_Years_50,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b1,bare1=bare1,LU=LU)
    F_SSM_low_check_50_t <- F_SSM_low_check_50[1] + F_SSM_low_check_50[2] + F_SSM_low_check_50[3] + F_SSM_low_check_50[4] + F_SSM_low_check_50[5]
    F_SSM_low_check_60<-Roth_C(Cinputs=(Cinputs*Low_Forest),years=check_Years_60,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b1,bare1=bare1,LU=LU)
    F_SSM_low_check_60_t <- F_SSM_low_check_60[1] + F_SSM_low_check_60[2] + F_SSM_low_check_60[3] + F_SSM_low_check_60[4] + F_SSM_low_check_60[5]
    F_SSM_low_check_70<-Roth_C(Cinputs=(Cinputs*Low_Forest),years=check_Years_70,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b1,bare1=bare1,LU=LU)
    F_SSM_low_check_70_t <- F_SSM_low_check_70[1] + F_SSM_low_check_70[2] + F_SSM_low_check_70[3] + F_SSM_low_check_70[4] + F_SSM_low_check_70[5]
    F_SSM_low_check_80<-Roth_C(Cinputs=(Cinputs*Low_Forest),years=check_Years_80,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b1,bare1=bare1,LU=LU)
    F_SSM_low_check_80_t <- F_SSM_low_check_80[1] + F_SSM_low_check_80[2] + F_SSM_low_check_80[3] + F_SSM_low_check_80[4] + F_SSM_low_check_80[5]
    F_SSM_low_check_90<-Roth_C(Cinputs=(Cinputs*Low_Forest),years=check_Years_90,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b1,bare1=bare1,LU=LU)
    F_SSM_low_check_90_t <- F_SSM_low_check_90[1] + F_SSM_low_check_90[2] + F_SSM_low_check_90[3] + F_SSM_low_check_90[4] + F_SSM_low_check_90[5]
    F_SSM_low_check_100<-Roth_C(Cinputs=(Cinputs*Low_Forest),years=check_Years_100,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b1,bare1=bare1,LU=LU)
    F_SSM_low_check_100_t <- F_SSM_low_check_100[1] + F_SSM_low_check_100[2] + F_SSM_low_check_100[3] + F_SSM_low_check_100[4] + F_SSM_low_check_100[5]
    
    
    F_SSM_med_check_10<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=check_Years_10,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    F_SSM_med_check_10_t <- F_SSM_med_check_10[1] + F_SSM_med_check_10[2] + F_SSM_med_check_10[3] + F_SSM_med_check_10[4] + F_SSM_med_check_10[5]
    F_SSM_med_check_20<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=check_Years_20,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    F_SSM_med_check_20_t <- F_SSM_med_check_20[1] + F_SSM_med_check_20[2] + F_SSM_med_check_20[3] + F_SSM_med_check_20[4] + F_SSM_med_check_20[5]
    F_SSM_med_check_30<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=check_Years_30,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    F_SSM_med_check_30_t <- F_SSM_med_check_30[1] + F_SSM_med_check_30[2] + F_SSM_med_check_30[3] + F_SSM_med_check_30[4] + F_SSM_med_check_30[5]
    F_SSM_med_check_40<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=check_Years_40,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    F_SSM_med_check_40_t <- F_SSM_med_check_40[1] + F_SSM_med_check_40[2] + F_SSM_med_check_40[3] + F_SSM_med_check_40[4] + F_SSM_med_check_40[5]
    F_SSM_med_check_50<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=check_Years_50,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    F_SSM_med_check_50_t <- F_SSM_med_check_50[1] + F_SSM_med_check_50[2] + F_SSM_med_check_50[3] + F_SSM_med_check_50[4] + F_SSM_med_check_50[5]
    F_SSM_med_check_60<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=check_Years_60,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    F_SSM_med_check_60_t <- F_SSM_med_check_60[1] + F_SSM_med_check_60[2] + F_SSM_med_check_60[3] + F_SSM_med_check_60[4] + F_SSM_med_check_60[5]
    F_SSM_med_check_70<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=check_Years_70,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    F_SSM_med_check_70_t <- F_SSM_med_check_70[1] + F_SSM_med_check_70[2] + F_SSM_med_check_70[3] + F_SSM_med_check_70[4] + F_SSM_med_check_70[5]
    F_SSM_med_check_80<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=check_Years_80,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    F_SSM_med_check_80_t <- F_SSM_med_check_80[1] + F_SSM_med_check_80[2] + F_SSM_med_check_80[3] + F_SSM_med_check_80[4] + F_SSM_med_check_80[5]
    F_SSM_med_check_90<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=check_Years_90,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    F_SSM_med_check_90_t <- F_SSM_med_check_90[1] + F_SSM_med_check_90[2] + F_SSM_med_check_90[3] + F_SSM_med_check_90[4] + F_SSM_med_check_90[5]
    F_SSM_med_check_100<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=check_Years_100,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    F_SSM_med_check_100_t <- F_SSM_med_check_100[1] + F_SSM_med_check_100[2] + F_SSM_med_check_100[3] + F_SSM_med_check_100[4] + F_SSM_med_check_100[5]
    
    
    F_SSM_high_check_10<-Roth_C(Cinputs=(Cinputs*High_Forest),years=check_Years_10,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b3,bare1=bare1,LU=LU)
    F_SSM_high_check_10_t <- F_SSM_high_check_10[1] + F_SSM_high_check_10[2] + F_SSM_high_check_10[3] + F_SSM_high_check_10[4] + F_SSM_high_check_10[5]
    F_SSM_high_check_20<-Roth_C(Cinputs=(Cinputs*High_Forest),years=check_Years_20,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b3,bare1=bare1,LU=LU)
    F_SSM_high_check_20_t <- F_SSM_high_check_20[1] + F_SSM_high_check_20[2] + F_SSM_high_check_20[3] + F_SSM_high_check_20[4] + F_SSM_high_check_20[5]
    F_SSM_high_check_30<-Roth_C(Cinputs=(Cinputs*High_Forest),years=check_Years_30,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b3,bare1=bare1,LU=LU)
    F_SSM_high_check_30_t <- F_SSM_high_check_30[1] + F_SSM_high_check_30[2] + F_SSM_high_check_30[3] + F_SSM_high_check_30[4] + F_SSM_high_check_30[5]
    F_SSM_high_check_40<-Roth_C(Cinputs=(Cinputs*High_Forest),years=check_Years_40,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b3,bare1=bare1,LU=LU)
    F_SSM_high_check_40_t <- F_SSM_high_check_40[1] + F_SSM_high_check_40[2] + F_SSM_high_check_40[3] + F_SSM_high_check_40[4] + F_SSM_high_check_40[5]
    F_SSM_high_check_50<-Roth_C(Cinputs=(Cinputs*High_Forest),years=check_Years_50,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b3,bare1=bare1,LU=LU)
    F_SSM_high_check_50_t <- F_SSM_high_check_50[1] + F_SSM_high_check_50[2] + F_SSM_high_check_50[3] + F_SSM_high_check_50[4] + F_SSM_high_check_50[5]
    F_SSM_high_check_60<-Roth_C(Cinputs=(Cinputs*High_Forest),years=check_Years_60,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b3,bare1=bare1,LU=LU)
    F_SSM_high_check_60_t <- F_SSM_high_check_60[1] + F_SSM_high_check_60[2] + F_SSM_high_check_60[3] + F_SSM_high_check_60[4] + F_SSM_high_check_60[5]
    F_SSM_high_check_70<-Roth_C(Cinputs=(Cinputs*High_Forest),years=check_Years_70,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b3,bare1=bare1,LU=LU)
    F_SSM_high_check_70_t <- F_SSM_high_check_70[1] + F_SSM_high_check_70[2] + F_SSM_high_check_70[3] + F_SSM_high_check_70[4] + F_SSM_high_check_70[5]
    F_SSM_high_check_80<-Roth_C(Cinputs=(Cinputs*High_Forest),years=check_Years_80,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b3,bare1=bare1,LU=LU)
    F_SSM_high_check_80_t <- F_SSM_high_check_80[1] + F_SSM_high_check_80[2] + F_SSM_high_check_80[3] + F_SSM_high_check_80[4] + F_SSM_high_check_80[5]
    F_SSM_high_check_90<-Roth_C(Cinputs=(Cinputs*High_Forest),years=check_Years_90,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b3,bare1=bare1,LU=LU)
    F_SSM_high_check_90_t <- F_SSM_high_check_90[1] + F_SSM_high_check_90[2] + F_SSM_high_check_90[3] + F_SSM_high_check_90[4] + F_SSM_high_check_90[5]
    F_SSM_high_check_100<-Roth_C(Cinputs=(Cinputs*High_Forest),years=check_Years_100,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b3,bare1=bare1,LU=LU)
    F_SSM_high_check_100_t <- F_SSM_high_check_100[1] + F_SSM_high_check_100[2] + F_SSM_high_check_100[3] + F_SSM_high_check_100[4] + F_SSM_high_check_100[5]
    
    
    
    
    #SSM forest unc min
    F_SSM_med_min_check_10 <-Roth_C(Cinputs=(Cinputs_min*(Low_Forest-0.15)),years=check_Years_10,DPMptf=WARM_UP[i,12], RPMptf=WARM_UP[i,13], BIOptf=WARM_UP[i,14], HUMptf=WARM_UP[i,15], FallIOM=WARM_UP[i,16],Temp=Temp*1.02,Precip=Precip*0.95,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*0.8,clay=clay*0.9,DR=DR_mean_b2,bare1=bare1,LU=LU)
    F_SSM_med_min_check_10_t <- F_SSM_med_min_check_10[1] + F_SSM_med_min_check_10[2] + F_SSM_med_min_check_10[3] + F_SSM_med_min_check_10[4] + F_SSM_med_min_check_10[5]
    F_SSM_med_min_check_20<-Roth_C(Cinputs=(Cinputs_min*(Low_Forest-0.15)),years=check_Years_20,DPMptf=WARM_UP[i,12], RPMptf=WARM_UP[i,13], BIOptf=WARM_UP[i,14], HUMptf=WARM_UP[i,15], FallIOM=WARM_UP[i,16],Temp=Temp*1.02,Precip=Precip*0.95,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*0.8,clay=clay*0.9,DR=DR_mean_b2,bare1=bare1,LU=LU)
    F_SSM_med_min_check_20_t <- F_SSM_med_min_check_20[1] + F_SSM_med_min_check_20[2] + F_SSM_med_min_check_20[3] + F_SSM_med_min_check_20[4] + F_SSM_med_min_check_20[5]
    F_SSM_med_min_check_30<-Roth_C(Cinputs=(Cinputs_min*(Low_Forest-0.15)),years=check_Years_30,DPMptf=WARM_UP[i,12], RPMptf=WARM_UP[i,13], BIOptf=WARM_UP[i,14], HUMptf=WARM_UP[i,15], FallIOM=WARM_UP[i,16],Temp=Temp*1.02,Precip=Precip*0.95,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*0.8,clay=clay*0.9,DR=DR_mean_b2,bare1=bare1,LU=LU)
    F_SSM_med_min_check_30_t <- F_SSM_med_min_check_30[1] + F_SSM_med_min_check_30[2] + F_SSM_med_min_check_30[3] + F_SSM_med_min_check_30[4] + F_SSM_med_min_check_30[5]
    F_SSM_med_min_check_40<-Roth_C(Cinputs=(Cinputs_min*(Low_Forest-0.15)),years=check_Years_40,DPMptf=WARM_UP[i,12], RPMptf=WARM_UP[i,13], BIOptf=WARM_UP[i,14], HUMptf=WARM_UP[i,15], FallIOM=WARM_UP[i,16],Temp=Temp*1.02,Precip=Precip*0.95,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*0.8,clay=clay*0.9,DR=DR_mean_b2,bare1=bare1,LU=LU)
    F_SSM_med_min_check_40_t <- F_SSM_med_min_check_40[1] + F_SSM_med_min_check_40[2] + F_SSM_med_min_check_40[3] + F_SSM_med_min_check_40[4] + F_SSM_med_min_check_40[5]
    F_SSM_med_min_check_50<-Roth_C(Cinputs=(Cinputs_min*(Low_Forest-0.15)),years=check_Years_50,DPMptf=WARM_UP[i,12], RPMptf=WARM_UP[i,13], BIOptf=WARM_UP[i,14], HUMptf=WARM_UP[i,15], FallIOM=WARM_UP[i,16],Temp=Temp*1.02,Precip=Precip*0.95,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*0.8,clay=clay*0.9,DR=DR_mean_b2,bare1=bare1,LU=LU)
    F_SSM_med_min_check_50_t <- F_SSM_med_min_check_50[1] + F_SSM_med_min_check_50[2] + F_SSM_med_min_check_50[3] + F_SSM_med_min_check_50[4] + F_SSM_med_min_check_50[5]
    F_SSM_med_min_check_60<-Roth_C(Cinputs=(Cinputs_min*(Low_Forest-0.15)),years=check_Years_60,DPMptf=WARM_UP[i,12], RPMptf=WARM_UP[i,13], BIOptf=WARM_UP[i,14], HUMptf=WARM_UP[i,15], FallIOM=WARM_UP[i,16],Temp=Temp*1.02,Precip=Precip*0.95,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*0.8,clay=clay*0.9,DR=DR_mean_b2,bare1=bare1,LU=LU)
    F_SSM_med_min_check_60_t <- F_SSM_med_min_check_60[1] + F_SSM_med_min_check_60[2] + F_SSM_med_min_check_60[3] + F_SSM_med_min_check_60[4] + F_SSM_med_min_check_60[5]
    F_SSM_med_min_check_70<-Roth_C(Cinputs=(Cinputs_min*(Low_Forest-0.15)),years=check_Years_70,DPMptf=WARM_UP[i,12], RPMptf=WARM_UP[i,13], BIOptf=WARM_UP[i,14], HUMptf=WARM_UP[i,15], FallIOM=WARM_UP[i,16],Temp=Temp*1.02,Precip=Precip*0.95,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*0.8,clay=clay*0.9,DR=DR_mean_b2,bare1=bare1,LU=LU)
    F_SSM_med_min_check_70_t <- F_SSM_med_min_check_70[1] + F_SSM_med_min_check_70[2] + F_SSM_med_min_check_70[3] + F_SSM_med_min_check_70[4] + F_SSM_med_min_check_70[5]
    F_SSM_med_min_check_80<-Roth_C(Cinputs=(Cinputs_min*(Low_Forest-0.15)),years=check_Years_80,DPMptf=WARM_UP[i,12], RPMptf=WARM_UP[i,13], BIOptf=WARM_UP[i,14], HUMptf=WARM_UP[i,15], FallIOM=WARM_UP[i,16],Temp=Temp*1.02,Precip=Precip*0.95,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*0.8,clay=clay*0.9,DR=DR_mean_b2,bare1=bare1,LU=LU)
    F_SSM_med_min_check_80_t <- F_SSM_med_min_check_80[1] + F_SSM_med_min_check_80[2] + F_SSM_med_min_check_80[3] + F_SSM_med_min_check_80[4] + F_SSM_med_min_check_80[5]
    F_SSM_med_min_check_90<-Roth_C(Cinputs=(Cinputs_min*(Low_Forest-0.15)),years=check_Years_90,DPMptf=WARM_UP[i,12], RPMptf=WARM_UP[i,13], BIOptf=WARM_UP[i,14], HUMptf=WARM_UP[i,15], FallIOM=WARM_UP[i,16],Temp=Temp*1.02,Precip=Precip*0.95,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*0.8,clay=clay*0.9,DR=DR_mean_b2,bare1=bare1,LU=LU)
    F_SSM_med_min_check_90_t <- F_SSM_med_min_check_90[1] + F_SSM_med_min_check_90[2] + F_SSM_med_min_check_90[3] + F_SSM_med_min_check_90[4] + F_SSM_med_min_check_90[5]
    F_SSM_med_min_check_100<-Roth_C(Cinputs=(Cinputs_min*(Low_Forest-0.15)),years=check_Years_100,DPMptf=WARM_UP[i,12], RPMptf=WARM_UP[i,13], BIOptf=WARM_UP[i,14], HUMptf=WARM_UP[i,15], FallIOM=WARM_UP[i,16],Temp=Temp*1.02,Precip=Precip*0.95,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*0.8,clay=clay*0.9,DR=DR_mean_b2,bare1=bare1,LU=LU)
    F_SSM_med_min_check_100_t <- F_SSM_med_min_check_100[1] + F_SSM_med_min_check_100[2] + F_SSM_med_min_check_100[3] + F_SSM_med_min_check_100[4] + F_SSM_med_min_check_100[5]
    
    
    #SSM forest unc max
    F_SSM_med_max_check_10<-Roth_C(Cinputs=(Cinputs_max*(Med_Forest+0.15)),years=check_Years_10,DPMptf=WARM_UP[i,18], RPMptf=WARM_UP[i,19], BIOptf=WARM_UP[i,20], HUMptf=WARM_UP[i,21], FallIOM=WARM_UP[i,22],Temp=Temp*0.98,Precip=Precip*1.05,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*1.2,clay=clay*1.1,DR=DR_mean_b2,bare1=bare1,LU=LU)
    F_SSM_med_max_check_10_t <- F_SSM_med_max_check_10[1] + F_SSM_med_max_check_10[2] + F_SSM_med_max_check_10[3] + F_SSM_med_max_check_10[4] + F_SSM_med_max_check_10[5]
    F_SSM_med_max_check_20<-Roth_C(Cinputs=(Cinputs_max*(Med_Forest+0.15)),years=check_Years_20,DPMptf=WARM_UP[i,18], RPMptf=WARM_UP[i,19], BIOptf=WARM_UP[i,20], HUMptf=WARM_UP[i,21], FallIOM=WARM_UP[i,22],Temp=Temp*0.98,Precip=Precip*1.05,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*1.2,clay=clay*1.1,DR=DR_mean_b2,bare1=bare1,LU=LU)
    F_SSM_med_max_check_20_t <- F_SSM_med_max_check_20[1] + F_SSM_med_max_check_20[2] + F_SSM_med_max_check_20[3] + F_SSM_med_max_check_20[4] + F_SSM_med_max_check_20[5]
    F_SSM_med_max_check_30<-Roth_C(Cinputs=(Cinputs_max*(Med_Forest+0.15)),years=check_Years_30,DPMptf=WARM_UP[i,18], RPMptf=WARM_UP[i,19], BIOptf=WARM_UP[i,20], HUMptf=WARM_UP[i,21], FallIOM=WARM_UP[i,22],Temp=Temp*0.98,Precip=Precip*1.05,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*1.2,clay=clay*1.1,DR=DR_mean_b2,bare1=bare1,LU=LU)
    F_SSM_med_max_check_30_t <- F_SSM_med_max_check_30[1] + F_SSM_med_max_check_30[2] + F_SSM_med_max_check_30[3] + F_SSM_med_max_check_30[4] + F_SSM_med_max_check_30[5]
    F_SSM_med_max_check_40<-Roth_C(Cinputs=(Cinputs_max*(Med_Forest+0.15)),years=check_Years_40,DPMptf=WARM_UP[i,18], RPMptf=WARM_UP[i,19], BIOptf=WARM_UP[i,20], HUMptf=WARM_UP[i,21], FallIOM=WARM_UP[i,22],Temp=Temp*0.98,Precip=Precip*1.05,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*1.2,clay=clay*1.1,DR=DR_mean_b2,bare1=bare1,LU=LU)
    F_SSM_med_max_check_40_t <- F_SSM_med_max_check_40[1] + F_SSM_med_max_check_40[2] + F_SSM_med_max_check_40[3] + F_SSM_med_max_check_40[4] + F_SSM_med_max_check_40[5]
    F_SSM_med_max_check_50<-Roth_C(Cinputs=(Cinputs_max*(Med_Forest+0.15)),years=check_Years_50,DPMptf=WARM_UP[i,18], RPMptf=WARM_UP[i,19], BIOptf=WARM_UP[i,20], HUMptf=WARM_UP[i,21], FallIOM=WARM_UP[i,22],Temp=Temp*0.98,Precip=Precip*1.05,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*1.2,clay=clay*1.1,DR=DR_mean_b2,bare1=bare1,LU=LU)
    F_SSM_med_max_check_50_t <- F_SSM_med_max_check_50[1] + F_SSM_med_max_check_50[2] + F_SSM_med_max_check_50[3] + F_SSM_med_max_check_50[4] + F_SSM_med_max_check_50[5]
    F_SSM_med_max_check_60<-Roth_C(Cinputs=(Cinputs_max*(Med_Forest+0.15)),years=check_Years_60,DPMptf=WARM_UP[i,18], RPMptf=WARM_UP[i,19], BIOptf=WARM_UP[i,20], HUMptf=WARM_UP[i,21], FallIOM=WARM_UP[i,22],Temp=Temp*0.98,Precip=Precip*1.05,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*1.2,clay=clay*1.1,DR=DR_mean_b2,bare1=bare1,LU=LU)
    F_SSM_med_max_check_60_t <- F_SSM_med_max_check_60[1] + F_SSM_med_max_check_60[2] + F_SSM_med_max_check_60[3] + F_SSM_med_max_check_60[4] + F_SSM_med_max_check_60[5]
    F_SSM_med_max_check_70<-Roth_C(Cinputs=(Cinputs_max*(Med_Forest+0.15)),years=check_Years_70,DPMptf=WARM_UP[i,18], RPMptf=WARM_UP[i,19], BIOptf=WARM_UP[i,20], HUMptf=WARM_UP[i,21], FallIOM=WARM_UP[i,22],Temp=Temp*0.98,Precip=Precip*1.05,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*1.2,clay=clay*1.1,DR=DR_mean_b2,bare1=bare1,LU=LU)
    F_SSM_med_max_check_70_t <- F_SSM_med_max_check_70[1] + F_SSM_med_max_check_70[2] + F_SSM_med_max_check_70[3] + F_SSM_med_max_check_70[4] + F_SSM_med_max_check_70[5]
    F_SSM_med_max_check_80<-Roth_C(Cinputs=(Cinputs_max*(Med_Forest+0.15)),years=check_Years_80,DPMptf=WARM_UP[i,18], RPMptf=WARM_UP[i,19], BIOptf=WARM_UP[i,20], HUMptf=WARM_UP[i,21], FallIOM=WARM_UP[i,22],Temp=Temp*0.98,Precip=Precip*1.05,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*1.2,clay=clay*1.1,DR=DR_mean_b2,bare1=bare1,LU=LU)
    F_SSM_med_max_check_80_t <- F_SSM_med_max_check_80[1] + F_SSM_med_max_check_80[2] + F_SSM_med_max_check_80[3] + F_SSM_med_max_check_80[4] + F_SSM_med_max_check_80[5]
    F_SSM_med_max_check_90<-Roth_C(Cinputs=(Cinputs_max*(Med_Forest+0.15)),years=check_Years_90,DPMptf=WARM_UP[i,18], RPMptf=WARM_UP[i,19], BIOptf=WARM_UP[i,20], HUMptf=WARM_UP[i,21], FallIOM=WARM_UP[i,22],Temp=Temp*0.98,Precip=Precip*1.05,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*1.2,clay=clay*1.1,DR=DR_mean_b2,bare1=bare1,LU=LU)
    F_SSM_med_max_check_90_t <- F_SSM_med_max_check_90[1] + F_SSM_med_max_check_90[2] + F_SSM_med_max_check_90[3] + F_SSM_med_max_check_90[4] + F_SSM_med_max_check_90[5]
    F_SSM_med_max_check_100<-Roth_C(Cinputs=(Cinputs_max*(Med_Forest+0.15)),years=check_Years_100,DPMptf=WARM_UP[i,18], RPMptf=WARM_UP[i,19], BIOptf=WARM_UP[i,20], HUMptf=WARM_UP[i,21], FallIOM=WARM_UP[i,22],Temp=Temp*0.98,Precip=Precip*1.05,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*1.2,clay=clay*1.1,DR=DR_mean_b2,bare1=bare1,LU=LU)
    F_SSM_med_max_check_100_t <- F_SSM_med_max_check_100[1] + F_SSM_med_max_check_100[2] + F_SSM_med_max_check_100[3] + F_SSM_med_max_check_100[4] + F_SSM_med_max_check_100[5]
    
    
    
    
    #############ARKANS EVIL VISUALIZATION PLAAAAAAAAAAAAN#################################################################################
    
    #########################################################
    #ARKAN: THIS IS MISSPELLED
    FOWARD[i,2]<-SOC
    FOWARD[i,3]<-f_bau_t
    FOWARD[i,4]<-f_bau[1]
    FOWARD[i,5]<-f_bau[2]
    FOWARD[i,6]<-f_bau[3]
    FOWARD[i,7]<-f_bau[4]
    FOWARD[i,8]<-f_bau[5]
    FOWARD[i,9]<-LU
    FOWARD[i,10]<-f_low_t
    FOWARD[i,11]<-f_med_t
    FOWARD[i,12]<-f_high_t
    FOWARD[i,13]<-f_bau_t_min
    FOWARD[i,14]<-f_bau_t_max
    FOWARD[i,15]<-f_med_t_min
    FOWARD[i,16]<-f_med_t_max
    FOWARD[i,17]<-SOC_t0_min
    FOWARD[i,18]<-SOC_t0_max
    
    
    
    
    check_points_list <- c(f_bau_check_10_t, f_bau_check_20_t, f_bau_check_30_t, f_bau_check_40_t, f_bau_check_50_t, f_bau_check_60_t, f_bau_check_70_t, f_bau_check_80_t, f_bau_check_90_t, f_bau_check_100_t,
                           f_bau_min_check_10_t, f_bau_min_check_20_t, f_bau_min_check_30_t, f_bau_min_check_40_t, f_bau_min_check_50_t, f_bau_min_check_60_t, f_bau_min_check_70_t, f_bau_min_check_80_t, f_bau_min_check_90_t, f_bau_min_check_100_t,
                           f_bau_max_check_10_t, f_bau_max_check_20_t, f_bau_max_check_30_t, f_bau_max_check_40_t, f_bau_max_check_50_t, f_bau_max_check_60_t, f_bau_max_check_70_t, f_bau_max_check_80_t, f_bau_max_check_90_t, f_bau_max_check_100_t,
                           F_SSM_low_check_10_t, F_SSM_low_check_20_t, F_SSM_low_check_30_t, F_SSM_low_check_40_t, F_SSM_low_check_50_t, F_SSM_low_check_60_t, F_SSM_low_check_70_t, F_SSM_low_check_80_t, F_SSM_low_check_90_t, F_SSM_low_check_100_t,
                           F_SSM_med_check_10_t, F_SSM_med_check_20_t, F_SSM_med_check_30_t, F_SSM_med_check_40_t, F_SSM_med_check_50_t, F_SSM_med_check_60_t, F_SSM_med_check_70_t, F_SSM_med_check_80_t, F_SSM_med_check_90_t, F_SSM_med_check_100_t,
                           F_SSM_high_check_10_t, F_SSM_high_check_20_t, F_SSM_high_check_30_t, F_SSM_high_check_40_t, F_SSM_high_check_50_t, F_SSM_high_check_60_t, F_SSM_high_check_70_t, F_SSM_high_check_80_t, F_SSM_high_check_90_t, F_SSM_high_check_100_t,
                           F_SSM_med_min_check_10_t, F_SSM_med_min_check_20_t, F_SSM_med_min_check_30_t, F_SSM_med_min_check_40_t, F_SSM_med_min_check_50_t, F_SSM_med_min_check_60_t, F_SSM_med_min_check_70_t, F_SSM_med_min_check_80_t, F_SSM_med_min_check_90_t, F_SSM_med_min_check_100_t,
                           F_SSM_med_max_check_10_t, F_SSM_med_max_check_20_t, F_SSM_med_max_check_30_t, F_SSM_med_max_check_40_t,F_SSM_med_max_check_50_t, F_SSM_med_max_check_60_t, F_SSM_med_max_check_70_t, F_SSM_med_max_check_80_t, F_SSM_med_max_check_90_t, F_SSM_med_max_check_100_t
    )
    for (x in 1:length(check_points_list+1)) {
      check_Points[i, x+1] <- check_points_list[x]
    }
    
    print(c(i,SOC,f_bau_t,f_low_t,f_med_t,f_high_t,f_bau_t_min,f_bau_t_max))
    
  }
}

############for loop ends##############

names(FOWARD)[2]="SOC_t0"
names(FOWARD)[3]="SOC_BAU_20"
names(FOWARD)[4]="DPM_BAU_20"
names(FOWARD)[5]="RPM_BAU_20"
names(FOWARD)[6]="BIO_BAU_20"
names(FOWARD)[7]="HUM_BAU_20"
names(FOWARD)[8]="IOM_BAU_20"
names(FOWARD)[9]="LandUse"
names(FOWARD)[10]="Low_Scenario"
names(FOWARD)[11]="Med_Scenario"
names(FOWARD)[12]="High_Scenario"
names(FOWARD)[13]="SOC_BAU_20_min"
names(FOWARD)[14]="SOC_BAU_20_max"
names(FOWARD)[15]="Med_Scen_min"
names(FOWARD)[16]="Med_Scen_max"
names(FOWARD)[17]="SOC_t0_min"
names(FOWARD)[18]="SOC_t0_max"




check_points_Names <- c('f_bau_check_10_t', 'f_bau_check_20_t', 'f_bau_check_30_t', 'f_bau_check_40_t', 'f_bau_check_50_t', 'f_bau_check_60_t', 'f_bau_check_70_t', 'f_bau_check_80_t', 'f_bau_check_90_t', 'f_bau_check_100_t',
                       'f_bau_min_check_10_t', 'f_bau_min_check_20_t', 'f_bau_min_check_30_t', 'f_bau_min_check_40_t', 'f_bau_min_check_50_t', 'f_bau_min_check_60_t', 'f_bau_min_check_70_t', 'f_bau_min_check_80_t', 'f_bau_min_check_90_t', 'f_bau_min_check_100_t',
                       'f_bau_max_check_10_t', 'f_bau_max_check_20_t', 'f_bau_max_check_30_t', 'f_bau_max_check_40_t', 'f_bau_max_check_50_t', 'f_bau_max_check_60_t', 'f_bau_max_check_70_t', 'f_bau_max_check_80_t', 'f_bau_max_check_90_t', 'f_bau_max_check_100_t',
                       'F_SSM_low_check_10_t', 'F_SSM_low_check_20_t', 'F_SSM_low_check_30_t', 'F_SSM_low_check_40_t', 'F_SSM_low_check_50_t', 'F_SSM_low_check_60_t', 'F_SSM_low_check_70_t', 'F_SSM_low_check_80_t', 'F_SSM_low_check_90_t', 'F_SSM_low_check_100_t',
                       'F_SSM_med_check_10_t', 'F_SSM_med_check_20_t', 'F_SSM_med_check_30_t', 'F_SSM_med_check_40_t', 'F_SSM_med_check_50_t', 'F_SSM_med_check_60_t', 'F_SSM_med_check_70_t', 'F_SSM_med_check_80_t', 'F_SSM_med_check_90_t', 'F_SSM_med_check_100_t',
                       'F_SSM_high_check_10_t', 'F_SSM_high_check_20_t', 'F_SSM_high_check_30_t', 'F_SSM_high_check_40_t', 'F_SSM_high_check_50_t', 'F_SSM_high_check_60_t', 'F_SSM_high_check_70_t', 'F_SSM_high_check_80_t', 'F_SSM_high_check_90_t', 'F_SSM_high_check_100_t',
                       'F_SSM_med_min_check_10_t', 'F_SSM_med_min_check_20_t', 'F_SSM_med_min_check_30_t', 'F_SSM_med_min_check_40_t', 'F_SSM_med_min_check_50_t', 'F_SSM_med_min_check_60_t', 'F_SSM_med_min_check_70_t', 'F_SSM_med_min_check_80_t', 'F_SSM_med_min_check_90_t', 'F_SSM_med_min_check_100_t',
                       'F_SSM_med_max_check_10_t', 'F_SSM_med_max_check_20_t', 'F_SSM_med_max_check_30_t', 'F_SSM_med_max_check_40_t','F_SSM_med_max_check_50_t', 'F_SSM_med_max_check_60_t', 'F_SSM_med_max_check_70_t', 'F_SSM_med_max_check_80_t', 'F_SSM_med_max_check_90_t', 'F_SSM_med_max_check_100_t')

for (x in 2:length(check_points_list+1)) {
  names(check_Points)[x] <- check_points_Names[x]
}


# Eliminate  values out of range
FOWARD$SOC_BAU_20[FOWARD$SOC_BAU_20<0]<-NA
FOWARD$Low_Scenario[FOWARD$Low_Scenario<0]<-NA
FOWARD$Med_Scenario[FOWARD$Med_Scenario<0]<-NA
FOWARD$High_Scenario[FOWARD$High_Scenario<0]<-NA
FOWARD$Med_Scen_min[FOWARD$Med_Scen_min<0]<-NA
FOWARD$Med_Scen_max[FOWARD$Med_Scen_max<0]<-NA

FOWARD$SOC_BAU_20[FOWARD$SOC_BAU_20>300]<-NA
FOWARD$Low_Scenario[FOWARD$Low_Scenario>300]<-NA
FOWARD$Med_Scenario[FOWARD$Med_Scenario>300]<-NA
FOWARD$High_Scenario[FOWARD$High_Scenario>300]<-NA
FOWARD$Med_Scen_min[FOWARD$Med_Scen_min>300]<-NA
FOWARD$Med_Scen_max[FOWARD$Med_Scen_max>300]<-NA

# Set the working directory 

setwd(WD_OUT)

# UNCERTAINTIES



UNC_SOC<-((FOWARD$SOC_BAU_20_max-FOWARD$SOC_BAU_20_min)/(2*FOWARD$SOC_BAU_20))*100

UNC_t0<-((FOWARD$SOC_t0_max-FOWARD$SOC_t0_min)/(2*FOWARD$SOC_t0))*100

UNC_SSM<-((FOWARD$Med_Scen_max-FOWARD$Med_Scen_min)/(2*FOWARD$Med_Scenario))*100

FOWARD[, 19]<-UNC_SOC
FOWARD[, 20]<-UNC_t0
FOWARD[, 21]<-UNC_SSM

names(FOWARD)[19]="UNC_BAU"
names(FOWARD)[20]="UNC_t0"
names(FOWARD)[21]="UNC_SSM"


# SAVE the Points (shapefile)

writeVector(check_Points, filename="check_Points_AOI_100y_divided", filetype="ESRI Shapefile") 
