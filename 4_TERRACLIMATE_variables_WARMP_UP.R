library(terra)

WD <- ("H:/BisonRothC/INPUTS/TERRA_CLIMATE_1")
setwd(WD)

# OPEN LAYERS


# Open the TerraClimate data from GEE

tmp<-rast(c("AverageTemperature_2001-2021.tif"))

pre_01_18<-rast(c("Precipitation_2001-2021.tif"))

pet_01_18<-rast(c("PET_2001-2021.tif"))

# TEMPERATURE


# Get one month temperature ( January)

tmp_Jan_1<-tmp[[1]]

dim(tmp_Jan_1)

# Create empty list
Rlist<-list()

# Average of 20 years (j)  and 12 months (i) 
##########for loop starts###############
for (i in 1:12) { 
  var_sum<-tmp_Jan_1*0
  k<-i
  
  for (j in 1:(dim(tmp)[3]/12)) {
    print(k)
    var_sum<-(var_sum + tmp[[k]])
    
    k<-k+12
    
  }
  #Save each month average. 
  
  var_avg<-var_sum/(dim(tmp)[3]/12)
  
  #writeRaster(ra,filename=name, format="GTiff")
  Rlist[[i]]<-var_avg
}
##########for loop ends#############
#save a stack of months averages

Temp_Stack<-rast(c(Rlist))
Temp_Stack<-Temp_Stack*0.1 # rescale to C
writeRaster(x=Temp_Stack,filename='Temp_Stack_01-19_TC.tif',overwrite=TRUE)



#############################################################################################################################

#PRECIPITATION


# Have one month Precipitation ( January)

pre_Jan_1<-pre_01_18[[1]]

dim(pre_Jan_1)

# Create empty list
Rlist<-list()


# Average of 20 years (j)  and 12 months (i) 

#########for loop starts############
for (i in 1:12) { 
  
  var_sum<-pre_Jan_1*0
  k<-i
  
  for (j in 1:(dim(pre_01_18)[3]/12)) {
    print(k)
    var_sum<-(var_sum + pre_01_18[[k]])
    
    k<-k+12
    
  }
  #Save each month average. 
  
  var_avg<-var_sum/(dim(pre_01_18)[3]/12)
  
  #writeRaster(ra,filename=name, format="GTiff",overwrite=TRUE)
  Rlist[[i]]<-var_avg
}
##########for loop ends##########

#save a stack of months averages

Prec_Stack<-rast(c(Rlist))
writeRaster(Prec_Stack,filename='Prec_Stack_01-19_TC.tif',overwrite=TRUE)


########################################################################

# POTENTIAL EVAPOTRANSPIRATION 

# Have one month ETP ( January)

pet_Jan_1<-pet_01_18[[1]]

dim(pet_Jan_1)

# Create empty list
Rlist<-list()

# Average of 18 years (j)  and 12 months (i) 
############for loop starts##############
for (i in 1:12) { 
  
  var_sum<-pet_Jan_1*0
  k<-i
  
  for (j in 1:(dim(pet_01_18)[3]/12)) {
    print(k)
    var_sum<-(var_sum + pet_01_18[[k]])
    
    k<-k+12
    
  }
  #Save each month average. 
  
  var_avg<-var_sum/(dim(pet_01_18)[3]/12)
  
  #writeRaster(ra,filename=name, format="GTiff",overwrite=TRUE)
  Rlist[[i]]<-var_avg
}
#########for loop ends############

#save a stack of months averages

PET_Stack<-rast(c(Rlist))
PET_Stack<-PET_Stack*0.1
writeRaster(PET_Stack,filename='PET_Stack_01-19_TC.tif',overwrite=TRUE)
