library(terra)
library(SoilR)
library(soilassessment)

WD_OUT<-("H:/BisonRothC/OUTPUTS/3_FORWARD")
working_dir<-setwd("H:/BisonRothC")
# OPEN THE VECTOR OF POINTS
Vector<-vect("INPUTS/TARGET_POINTS/target_points.shp")
# OPEN THE RESULT VECTOR FROM THE WARM UP PROCESS
WARM_UP<-vect("OUTPUTS/2_WARM_UP/WARM_UP_AOI2/WARM_UP_AOI2.shp")
# OPEN THE STACK WITH THE VARIABLES FOR THE FOWARD PROCESS
Stack_Set_1<- rast("INPUTS/STACK/Stack_Set_FOWARD1.tif")

# Set the increase in Carbon input for each land use and each scenario
#Crops and Crop trees
Low_Crops<-1.05
Med_Crops<-1.10
High_Crops<-1.2
#Shrublands, Grasslands , Herbaceous vegetation flooded & Sparse Vegetation
Low_Grass<-1.05
Med_Grass<-1.10
High_Grass<-1.2
#Forest
Low_Forest<-1.05
Med_Forest<-1.10
High_Forest<-1.2
#Paddy Fields
Low_PaddyFields<-1.05
Med_PaddyFields<-1.10
High_PaddyFields<-1.2

# extract variables to points
Variables<-extract(Stack_Set_1,Vector,bind=TRUE)

# Creates an empty vector
FOWARD<-Vector

# Extract the layers from the Vector
SOC_im<-WARM_UP[[4]]
clay_im<-Variables[[3]] 
Cinputs_im<-WARM_UP[[10]] 
DR_im<-Variables[[40]]
LU_im<-Variables[[41]]

# Define the years to run the model
years=seq(1/12,20,by=1/12)



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
  
  if (any(is.na(Evp[,2])) | any(is.na(Temp[,2])) | any(is.na(SOC_im[i, ])) | any(is.na(clay_im[i, ])) | any(is.na(Precip[,2]))  |  any(is.na(Cov2[,2]))  |  any(is.na(Cov1[,1])) | any(is.na(Cinputs_im[i, ])) | any(is.na(DR_im[i, ])) | (Cinputs_im[i, ]<0) |  (SOC_im[i, ]<0) | (clay_im[i, ]<0) ) {FOWARD[i,2]<-0}else{
    
    
    # Set the variables from the images
    
    soil.thick=30  #Soil thickness (organic layer topsoil), in cm
    SOC<-SOC_im[i, ]      #Soil organic carbon in Mg/ha 
    clay<-clay_im[i, ]        #Percent clay %
    Cinputs<-Cinputs_im[i, ]    #Annual C inputs to soil in Mg/ha/yr
    
    DR<-DR_im[i, ]              # DPM/RPM (decomplosable vs resistant plant material.)
    bare1<-(Cov1>0.8)           # If the surface is bare or vegetated
    LU<-LU_im[i, ]
    
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
    
    # Crops and Tree crops
    if (LU==2 | LU==12){
      f_low<-Roth_C(Cinputs=(Cinputs*Low_Crops),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR,bare1=bare1,LU=LU)
      f_low_t<-f_low[1]+f_low[2]+f_low[3]+f_low[4]+f_low[5]
      
      f_med<-Roth_C(Cinputs=(Cinputs*Med_Crops),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR,bare1=bare1,LU=LU)
      f_med_t<-f_med[1]+f_med[2]+f_med[3]+f_med[4]+f_med[5]
      
      f_high<-Roth_C(Cinputs=(Cinputs*High_Crops),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR,bare1=bare1,LU=LU)
      f_high_t<-f_high[1]+f_high[2]+f_high[3]+f_high[4]+f_high[5]
      
      # SSM croplands unc min
      
      f_med_min<-Roth_C(Cinputs=(Cinputs_min*(Med_Crops-0.15)),years=years,DPMptf=WARM_UP[i,12], RPMptf=WARM_UP[i,13], BIOptf=WARM_UP[i,14], HUMptf=WARM_UP[i,15], FallIOM=WARM_UP[i,16],Temp=Temp*1.02,Precip=Precip*0.95,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*0.8,clay=clay*0.9,DR=DR,bare1=bare1,LU=LU)
      f_med_t_min<-f_med_min[1]+f_med_min[2]+f_med_min[3]+f_med_min[4]+f_med_min[5]
      
      # SSM croplands unc max
      
      f_med_max<-Roth_C(Cinputs=(Cinputs_max*(Med_Crops+0.15)),years=years,DPMptf=WARM_UP[i,18], RPMptf=WARM_UP[i,19], BIOptf=WARM_UP[i,20], HUMptf=WARM_UP[i,21], FallIOM=WARM_UP[i,22],Temp=Temp*0.98,Precip=Precip*1.05,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*1.2,clay=clay*1.1,DR=DR,bare1=bare1,LU=LU)
      f_med_t_max<-f_med_max[1]+f_med_max[2]+f_med_max[3]+f_med_max[4]+f_med_max[5]
      
    }
    #Shrublands, grasslands, and sparce vegetation
    else if (LU==3 | LU==5 | LU==6 | LU==8) {
      f_low<-Roth_C(Cinputs=(Cinputs*Low_Grass),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR,bare1=bare1,LU=LU)
      f_low_t<-f_low[1]+f_low[2]+f_low[3]+f_low[4]+f_low[5]
      
      f_med<-Roth_C(Cinputs=(Cinputs*Med_Grass),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR,bare1=bare1,LU=LU)
      f_med_t<-f_med[1]+f_med[2]+f_med[3]+f_med[4]+f_med[5]
      
      f_high<-Roth_C(Cinputs=(Cinputs*High_Grass),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR,bare1=bare1,LU=LU)
      f_high_t<-f_high[1]+f_high[2]+f_high[3]+f_high[4]+f_high[5]
      
      #SSM Shrublands unc min
      
      f_med_min<-Roth_C(Cinputs=(Cinputs_min*(Med_Grass-0.15)),years=years,DPMptf=WARM_UP[i,12], RPMptf=WARM_UP[i,13], BIOptf=WARM_UP[i,14], HUMptf=WARM_UP[i,15], FallIOM=WARM_UP[i,16],Temp=Temp*1.02,Precip=Precip*0.95,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*0.8,clay=clay*0.9,DR=DR,bare1=bare1,LU=LU)
      f_med_t_min<-f_med_min[1]+f_med_min[2]+f_med_min[3]+f_med_min[4]+f_med_min[5]
      
      #SSM Shrublands unc max
      
      f_med_max<-Roth_C(Cinputs=(Cinputs_max*(Med_Grass+0.15)),years=years,DPMptf=WARM_UP[i,18], RPMptf=WARM_UP[i,19], BIOptf=WARM_UP[i,20], HUMptf=WARM_UP[i,21], FallIOM=WARM_UP[i,22],Temp=Temp*0.98,Precip=Precip*1.05,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*1.2,clay=clay*1.1,DR=DR,bare1=bare1,LU=LU)
      f_med_t_max<-f_med_max[1]+f_med_max[2]+f_med_max[3]+f_med_max[4]+f_med_max[5]
      
    }
    # Paddy Fields 
    else if (LU==13) {
      f_low<-Roth_C(Cinputs=(Cinputs*Low_PaddyFields),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR,bare1=bare1,LU=LU)
      f_low_t<-f_low[1]+f_low[2]+f_low[3]+f_low[4]+f_low[5]
      
      f_med<-Roth_C(Cinputs=(Cinputs*Med_PaddyFields),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR,bare1=bare1,LU=LU)
      f_med_t<-f_med[1]+f_med[2]+f_med[3]+f_med[4]+f_med[5]
      
      f_high<-Roth_C(Cinputs=(Cinputs*High_PaddyFields),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR,bare1=bare1,LU=LU)
      f_high_t<-f_high[1]+f_high[2]+f_high[3]+f_high[4]+f_high[5]
      
      #SSM Forest unc min
      
      f_med_min<-Roth_C(Cinputs=(Cinputs_min*(Med_PaddyFields-0.15)),years=years,DPMptf=WARM_UP[i,12], RPMptf=WARM_UP[i,13], BIOptf=WARM_UP[i,14], HUMptf=WARM_UP[i,15], FallIOM=WARM_UP[i,16],Temp=Temp*1.02,Precip=Precip*0.95,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*0.8,clay=clay*0.9,DR=DR,bare1=bare1,LU=LU)
      f_med_t_min<-f_med_min[1]+f_med_min[2]+f_med_min[3]+f_med_min[4]+f_med_min[5]
      
      #SSM Forest unc max
      
      f_med_max<-Roth_C(Cinputs=(Cinputs_max*(Med_PaddyFields+0.15)),years=years,DPMptf=WARM_UP[i,18], RPMptf=WARM_UP[i,19], BIOptf=WARM_UP[i,20], HUMptf=WARM_UP[i,21], FallIOM=WARM_UP[i,22],Temp=Temp*0.98,Precip=Precip*1.05,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*1.2,clay=clay*1.1,DR=DR,bare1=bare1,LU=LU)
      f_med_t_max<-f_med_max[1]+f_med_max[2]+f_med_max[3]+f_med_max[4]+f_med_max[5]
    }
    
    
    
    # Arkan's forest implementation:
    #Forest cover
    else if (LU==4 | LU==12) {
      f_low<-Roth_C(Cinputs=(Cinputs*Low_Forest),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR,bare1=bare1,LU=LU)
      f_low_t<-f_low[1]+f_low[2]+f_low[3]+f_low[4]+f_low[5]
      
      f_med<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR,bare1=bare1,LU=LU)
      f_med_t<-f_med[1]+f_med[2]+f_med[3]+f_med[4]+f_med[5]
      
      f_high<-Roth_C(Cinputs=(Cinputs*High_Forest),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR,bare1=bare1,LU=LU)
      f_high_t<-f_high[1]+f_high[2]+f_high[3]+f_high[4]+f_high[5]
      
      #SSM forest unc min
      
      f_med_min<-Roth_C(Cinputs=(Cinputs_min*(Low_Forest-0.15)),years=years,DPMptf=WARM_UP[i,12], RPMptf=WARM_UP[i,13], BIOptf=WARM_UP[i,14], HUMptf=WARM_UP[i,15], FallIOM=WARM_UP[i,16],Temp=Temp*1.02,Precip=Precip*0.95,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*0.8,clay=clay*0.9,DR=DR,bare1=bare1,LU=LU)
      f_med_t_min<-f_med_min[1]+f_med_min[2]+f_med_min[3]+f_med_min[4]+f_med_min[5]
      
      #SSM forest unc max
      
      f_med_max<-Roth_C(Cinputs=(Cinputs_max*(Med_Forest+0.15)),years=years,DPMptf=WARM_UP[i,18], RPMptf=WARM_UP[i,19], BIOptf=WARM_UP[i,20], HUMptf=WARM_UP[i,21], FallIOM=WARM_UP[i,22],Temp=Temp*0.98,Precip=Precip*1.05,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*1.2,clay=clay*1.1,DR=DR,bare1=bare1,LU=LU)
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

writeVector(FOWARD, filename="FOWARD_AOI2", filetype="ESRI Shapefile") 
