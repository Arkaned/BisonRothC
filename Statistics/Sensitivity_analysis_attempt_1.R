library(terra)

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
sens_Values <- Vector

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
################function set up ends############









# Iterates over the area of interest
##################for loop starts###############
for (i in 2516) {
  
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
    
    
    # Arkan's forest implementation:
    #Forest cover
      
    stand_Bison<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    stand_bison_t<-stand_Bison[1]+stand_Bison[2]+stand_Bison[3]+stand_Bison[4]+stand_Bison[5]
      
    Cinp_min<-Roth_C(Cinputs=(Cinputs*Med_Forest*0.9),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    Cinp_min_t<-Cinp_min[1]+Cinp_min[2]+Cinp_min[3]+Cinp_min[4]+Cinp_min[5]
    
    Cinp_max<-Roth_C(Cinputs=(Cinputs*Med_Forest*1.1),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    Cinp_max_t<-Cinp_max[1]+Cinp_max[2]+Cinp_max[3]+Cinp_max[4]+Cinp_max[5]
    
    DPM_min<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=years,DPMptf=WARM_UP[i,5][[1]][[1]]*0.9, RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    DPM_min_t<-DPM_min[1]+DPM_min[2]+DPM_min[3]+DPM_min[4]+DPM_min[5]
    
    DPM_max<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=years,DPMptf=WARM_UP[i,5][[1]][[1]]*1.1, RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    DPM_max_t<-DPM_max[1]+DPM_max[2]+DPM_max[3]+DPM_max[4]+DPM_max[5]
    
    RPM_min<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6][[1]][[1]]*0.9, BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    RPM_min_t<-RPM_min[1]+RPM_min[2]+RPM_min[3]+RPM_min[4]+RPM_min[5]
    
    RPM_max<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6][[1]][[1]]*1.1, BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    RPM_max_t<-RPM_max[1]+RPM_max[2]+RPM_max[3]+RPM_max[4]+RPM_max[5]
    
    BIO_min<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7][[1]][[1]]*0.9, HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    BIO_min_t<-BIO_min[1]+BIO_min[2]+BIO_min[3]+BIO_min[4]+BIO_min[5]
    
    BIO_max<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7][[1]][[1]]*1.1, HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    BIO_max_t<-BIO_max[1]+BIO_max[2]+BIO_max[3]+BIO_max[4]+BIO_max[5]
    
    HUM_min<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8][[1]][[1]]*0.9, FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    HUM_min_t<-HUM_min[1]+HUM_min[2]+HUM_min[3]+HUM_min[4]+HUM_min[5]
    
    HUM_max<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8][[1]][[1]]*1.1, FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    HUM_max_t<-HUM_max[1]+HUM_max[2]+HUM_max[3]+HUM_max[4]+HUM_max[5]
    
    IOM_min<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9][[1]][[1]]*0.9,Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    IOM_min_t<-IOM_min[1]+IOM_min[2]+IOM_min[3]+IOM_min[4]+IOM_min[5]
    
    IOM_max<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9][[1]][[1]]*1.1,Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    IOM_max_t<-IOM_max[1]+IOM_max[2]+IOM_max[3]+IOM_max[4]+IOM_max[5]
    
    TEMP_min<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp*0.9,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    TEMP_min_t<-TEMP_min[1]+TEMP_min[2]+TEMP_min[3]+TEMP_min[4]+TEMP_min[5]
    
    TEMP_max<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp*1.1,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    TEMP_max_t<-TEMP_max[1]+TEMP_max[2]+TEMP_max[3]+TEMP_max[4]+TEMP_max[5]
    
    PREC_min<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip*0.9,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    PREC_min_t<-PREC_min[1]+PREC_min[2]+PREC_min[3]+PREC_min[4]+PREC_min[5]
    
    PREC_max<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip*1.1,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    PREC_max_t<-PREC_max[1]+PREC_max[2]+PREC_max[3]+PREC_max[4]+PREC_max[5]
    
    EVP_min<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp*0.9,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    EVP_min_t<-EVP_min[1]+EVP_min[2]+EVP_min[3]+EVP_min[4]+EVP_min[5]
    
    EVP_max<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp*1.1,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    EVP_max_t<-EVP_max[1]+EVP_max[2]+EVP_max[3]+EVP_max[4]+EVP_max[5]
    
    COV_min<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2*0.9,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    COV_min_t<-COV_min[1]+COV_min[2]+COV_min[3]+COV_min[4]+COV_min[5]
    
    COV_max<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2*1.1,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    COV_max_t<-COV_max[1]+COV_max[2]+COV_max[3]+COV_max[4]+COV_max[5]
    
    SOC_min<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*0.9,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    SOC_min_t<-SOC_min[1]+SOC_min[2]+SOC_min[3]+SOC_min[4]+SOC_min[5]
    
    SOC_max<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC*1.1,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    SOC_max_t<-SOC_min[1]+SOC_min[2]+SOC_min[3]+SOC_min[4]+SOC_min[5]
    
    CLY_min<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay*0.9,DR=DR_mean_b2,bare1=bare1,LU=LU)
    CLY_min_t<-CLY_min[1]+CLY_min[2]+CLY_min[3]+CLY_min[4]+CLY_min[5]
    
    CLY_max<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay*1.1,DR=DR_mean_b2,bare1=bare1,LU=LU)
    CLY_max_t<-CLY_min[1]+CLY_min[2]+CLY_min[3]+CLY_min[4]+CLY_min[5]
    
    DR_min<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2*0.9,bare1=bare1,LU=LU)
    DR_min_t<-DR_min[1]+DR_min[2]+DR_min[3]+DR_min[4]+DR_min[5]
    
    DR_max<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2*1.1,bare1=bare1,LU=LU)
    DR_max_t<-DR_max[1]+DR_max[2]+DR_max[3]+DR_max[4]+DR_max[5]
    
    
    bare1<-((Cov1*0.9)>0.8)
    
    BARE_mod_min<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    Bare_mod_min_t<-BARE_mod_min[1]+BARE_mod_min[2]+BARE_mod_min[3]+BARE_mod_min[4]+BARE_mod_min[5]
    
    bare1<-((Cov1*1.1)>0.8)
    
    Bare_mod_max<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    Bare_mod_max_t<-Bare_mod_max[1]+Bare_mod_max[2]+Bare_mod_max[3]+Bare_mod_max[4]+Bare_mod_max[5]
    
    bare1<-(Cov1>0.8)
    
    
    
    sens_Values[1,2]<-stand_bison_t
    sens_Values[1,3]<-Cinp_min_t
    sens_Values[1,4]<-Cinp_max_t
    sens_Values[1,5]<-DPM_min_t
    sens_Values[1,6]<-DPM_max_t
    sens_Values[1,7]<-RPM_min_t
    sens_Values[1,8]<-RPM_max_t
    sens_Values[1,9]<-BIO_min_t
    sens_Values[1,10]<-BIO_max_t
    sens_Values[1,11]<-HUM_min_t
    sens_Values[1,12]<-HUM_max_t
    sens_Values[1,13]<-IOM_min_t
    sens_Values[1,14]<-IOM_max_t
    sens_Values[1,15]<-TEMP_min_t
    sens_Values[1,16]<-TEMP_max_t
    sens_Values[1,17]<-PREC_min_t
    sens_Values[1,18]<-PREC_max_t
    sens_Values[1,19]<-EVP_min_t
    sens_Values[1,20]<-EVP_max_t
    sens_Values[1,21]<-COV_min_t
    sens_Values[1,22]<-COV_max_t
    sens_Values[1,23]<-SOC_min_t
    sens_Values[1,24]<-SOC_max_t
    sens_Values[1,25]<-CLY_min_t
    sens_Values[1,26]<-CLY_max_t
    sens_Values[1,27]<-DR_min_t
    sens_Values[1,28]<-DR_max_t
    sens_Values[1,29]<-Bare_mod_min_t
    sens_Values[1,30]<-Bare_mod_max_t
    
    
    
    #PLANT_DEC_min<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    #PLANT_DEC_min_t<-PLANT_DEC_min[1]+PLANT_DEC_min[2]+PLANT_DEC_min[3]+PLANT_DEC_min[4]+PLANT_DEC_min[5]
    
    #PLANT_DEC_max<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    #PLANT_DEC_max_t<-PLANT_DEC_max[1]+PLANT_DEC_max[2]+PLANT_DEC_max[3]+PLANT_DEC_max[4]+PLANT_DEC_max[5]
    
    #BISON_DEC_min<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    #BISON_DEC_min_t<-BISON_DEC_min[1]+BISON_DEC_min[2]+BISON_DEC_min[3]+BISON_DEC_min[4]+BISON_DEC_min[5]
    
    #BISON_DEC_max<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    #BISON_DEC_max_t<-BISON_DEC_max[1]+BISON_DEC_max[2]+BISON_DEC_max[3]+BISON_DEC_max[4]+BISON_DEC_max[5]
    
    #BISON_MASS_max<-Roth_C(Cinputs=(Cinputs*Med_Forest),years=years,DPMptf=WARM_UP[i,5], RPMptf=WARM_UP[i,6], BIOptf=WARM_UP[i,7], HUMptf=WARM_UP[i,8], FallIOM=WARM_UP[i,9],Temp=Temp,Precip=Precip,Evp=Evp,Cov=Cov,Cov1=Cov1,Cov2=Cov2,soil.thick=soil.thick,SOC=SOC,clay=clay,DR=DR_mean_b2,bare1=bare1,LU=LU)
    #BISON_MASS_t<-BISON_MASS_max[1]+BISON_MASS_max[2]+BISON_MASS_max[3]+BISON_MASS_max[4]+BISON_MASS_max[5]
    
    
    
    #############ARKANS EVIL VISUALIZATION PLAAAAAAAAAAAAN#################################################################################
    
    #########################################################
    #ARKAN: THIS IS MISSPELLED
    #FOWARD[i,2]<-SOC
    #FOWARD[i,3]<-f_bau_t
    #FOWARD[i,4]<-f_bau[1]
    #FOWARD[i,5]<-f_bau[2]
    #FOWARD[i,6]<-f_bau[3]
    #FOWARD[i,7]<-f_bau[4]
    #FOWARD[i,8]<-f_bau[5]
    #FOWARD[i,9]<-LU
    #FOWARD[i,10]<-f_low_t
    #FOWARD[i,11]<-f_med_t
    #FOWARD[i,12]<-f_high_t
    #FOWARD[i,13]<-f_bau_t_min
    #FOWARD[i,14]<-f_bau_t_max
    #FOWARD[i,15]<-f_med_t_min
    #FOWARD[i,16]<-f_med_t_max
    #FOWARD[i,17]<-SOC_t0_min
    #FOWARD[i,18]<-SOC_t0_max
    
    
    
    
    
    
    
    
    
    #print(c(i,SOC,f_bau_t,f_low_t,f_med_t,f_high_t,f_bau_t_min,f_bau_t_max))
    
  }
}

############for loop ends##############

names(sens_Values)[2]<-'stand_bison'
names(sens_Values)[3]<-'Cinp_min'
names(sens_Values)[4]<-'Cinp_max'
names(sens_Values)[5]<-'DPM_min'
names(sens_Values)[6]<-'DPM_max'
names(sens_Values)[7]<-'RPM_min'
names(sens_Values)[8]<-'RPM_max'
names(sens_Values)[9]<-'BIO_min'
names(sens_Values)[10]<-'BIO_max'
names(sens_Values)[11]<-'HUM_min'
names(sens_Values)[12]<-'HUM_max'
names(sens_Values)[13]<-'IOM_min'
names(sens_Values)[14]<-'IOM_max'
names(sens_Values)[15]<-'TEMP_min'
names(sens_Values)[16]<-'TEMP_max'
names(sens_Values)[17]<-'PREC_min'
names(sens_Values)[18]<-'PREC_max'
names(sens_Values)[19]<-'EVP_min'
names(sens_Values)[20]<-'EVP_max'
names(sens_Values)[21]<-'COV_min'
names(sens_Values)[22]<-'COV_max'
names(sens_Values)[23]<-'SOC_min'
names(sens_Values)[24]<-'SOC_max'
names(sens_Values)[25]<-'CLY_min'
names(sens_Values)[26]<-'CLY_max'
names(sens_Values)[27]<-'DR_min'
names(sens_Values)[28]<-'DR_max'
names(sens_Values)[29]<-'Bare_mod_min'
names(sens_Values)[30]<-'Bare_mod_max'



# Eliminate  values out of range
#FOWARD$SOC_BAU_20[FOWARD$SOC_BAU_20<0]<-NA
#FOWARD$Low_Scenario[FOWARD$Low_Scenario<0]<-NA
#FOWARD$Med_Scenario[FOWARD$Med_Scenario<0]<-NA
#FOWARD$High_Scenario[FOWARD$High_Scenario<0]<-NA
#FOWARD$Med_Scen_min[FOWARD$Med_Scen_min<0]<-NA
#FOWARD$Med_Scen_max[FOWARD$Med_Scen_max<0]<-NA

#FOWARD$SOC_BAU_20[FOWARD$SOC_BAU_20>300]<-NA
#FOWARD$Low_Scenario[FOWARD$Low_Scenario>300]<-NA
#FOWARD$Med_Scenario[FOWARD$Med_Scenario>300]<-NA
#FOWARD$High_Scenario[FOWARD$High_Scenario>300]<-NA
#FOWARD$Med_Scen_min[FOWARD$Med_Scen_min>300]<-NA
#FOWARD$Med_Scen_max[FOWARD$Med_Scen_max>300]<-NA

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

setwd("H:/BisonRothC/Statistics")





# SAVE the Points (shapefile)

writeVector(sens_Values, filename="Sensitivity_analysis_F_run", filetype="ESRI Shapefile", overwrite=TRUE) 
