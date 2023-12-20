library(terra)
library(sf)
library(ggplot2)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(reshape2)
library(scales)

rm(list=ls())

divided_data <- "H:/BisonRothC/OUTPUTS/3_FORWARD/check_Points_AOI_100y_divided"

dFrame <- "H:/BisonRothC/OUTPUTS/3_FORWARD/FOWARD_AOI_100y_divided/FOWARD_AOI_100y_divided.shp"

WARMY <- "H:/BisonRothC/OUTPUTS/2_WARM_UP/WARM_UP_AOI_final/WARM_UP_AOI_final.shp"

WARMZ <- vect(WARMY)
foward_dat <- vect(dFrame)

setwd(divided_data)

div_dat <- vect("check_Points_AOI_100y_divided.shp")

p <- 1

yearsx <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)

f_bau_check_list <- c()
for (cheese in 2:11) {
  f_bau_check_list[cheese-1] <- div_dat[1, cheese][[1]][1]
}

plot(yearsx, f_bau_check_list)



########## For box plotting ###########

box_plot_f_bau_10 <- c()
box_plot_f_bau_20 <- c()
box_plot_f_bau_30 <- c()
box_plot_f_bau_40 <- c()
box_plot_f_bau_50 <- c()
box_plot_f_bau_60 <- c()
box_plot_f_bau_70 <- c()
box_plot_f_bau_80 <- c()
box_plot_f_bau_90 <- c()
box_plot_f_bau_100 <- c()


for (val in 1:length(div_dat)){
  box_plot_f_bau_10[val] <- div_dat[, 2][[1]][[1]][val]
  box_plot_f_bau_20[val] <- div_dat[, 3][[1]][[1]][val]
  box_plot_f_bau_30[val] <- div_dat[, 4][[1]][[1]][val]
  box_plot_f_bau_40[val] <- div_dat[, 5][[1]][[1]][val]
  box_plot_f_bau_50[val] <- div_dat[, 6][[1]][[1]][val]
  box_plot_f_bau_60[val] <- div_dat[, 7][[1]][[1]][val]
  box_plot_f_bau_70[val] <- div_dat[, 8][[1]][[1]][val]
  box_plot_f_bau_80[val] <- div_dat[, 9][[1]][[1]][val]
  box_plot_f_bau_90[val] <- div_dat[, 10][[1]][[1]][val]
  box_plot_f_bau_100[val] <- div_dat[, 11][[1]][[1]][val]
  
  
}


boxing_data <- data.frame(box_plot_f_bau_10,
                         box_plot_f_bau_20,
                         box_plot_f_bau_30,
                         box_plot_f_bau_40,
                         box_plot_f_bau_50,
                         box_plot_f_bau_60,
                         box_plot_f_bau_70,
                         box_plot_f_bau_80,
                         box_plot_f_bau_90,
                         box_plot_f_bau_100)
boxplot(boxing_data)





######### COMPARING MEDS, MAXES AND MINS (average plus std#####################


####### maybe first make a sum of the differences between SSM and and BAU:


SOC_SSM_LOW_diff <- c()
SOC_SSM_MED_diff <- c()
SOC_SSM_HIGH_diff <- c()

box_plot_f_SSM_M_10 <- c()
box_plot_f_SSM_M_20 <- c()
box_plot_f_SSM_M_30 <- c()
box_plot_f_SSM_M_40 <- c()
box_plot_f_SSM_M_50 <- c()
box_plot_f_SSM_M_60 <- c()
box_plot_f_SSM_M_70 <- c()
box_plot_f_SSM_M_80 <- c()
box_plot_f_SSM_M_90 <- c()
box_plot_f_SSM_M_100 <- c()

for (val in 1:length(div_dat)){
  if(is.na(foward_dat[,2][[1]][[1]][val]) | foward_dat[,2][[1]][[1]][val] == 0) {next}

  box_plot_f_SSM_M_10[val] <- div_dat[, 42][[1]][[1]][val]- div_dat[, 2][[1]][[1]][val]
  box_plot_f_SSM_M_20[val] <- div_dat[, 43][[1]][[1]][val] - div_dat[, 3][[1]][[1]][val]
  box_plot_f_SSM_M_30[val] <- div_dat[, 44][[1]][[1]][val] - div_dat[, 4][[1]][[1]][val]
  box_plot_f_SSM_M_40[val] <- div_dat[, 45][[1]][[1]][val] - div_dat[, 5][[1]][[1]][val]
  box_plot_f_SSM_M_50[val] <- div_dat[, 46][[1]][[1]][val] - div_dat[, 6][[1]][[1]][val]
  box_plot_f_SSM_M_60[val] <- div_dat[, 47][[1]][[1]][val] - div_dat[, 7][[1]][[1]][val]
  box_plot_f_SSM_M_70[val] <- div_dat[, 48][[1]][[1]][val] - div_dat[, 8][[1]][[1]][val]
  box_plot_f_SSM_M_80[val] <- div_dat[, 49][[1]][[1]][val] - div_dat[, 9][[1]][[1]][val]
  box_plot_f_SSM_M_90[val] <- div_dat[, 50][[1]][[1]][val] - div_dat[, 10][[1]][[1]][val]
  box_plot_f_SSM_M_100[val] <- div_dat[, 51][[1]][[1]][val] - div_dat[, 11][[1]][[1]][val]
  if (val%%100 == 0) {print(val)}
}



boxing_data <- data.frame(box_plot_f_SSM_M_10,
                          box_plot_f_SSM_M_20,
                          box_plot_f_SSM_M_30,
                          box_plot_f_SSM_M_40,
                          box_plot_f_SSM_M_50,
                          box_plot_f_SSM_M_60,
                          box_plot_f_SSM_M_70,
                          box_plot_f_SSM_M_80,
                          box_plot_f_SSM_M_90,
                          box_plot_f_SSM_M_100)
boxplot(boxing_data)

with_line <- melt(boxing_data)

ggplot(data=with_line, aes(x=variable, y=value)) + 
  geom_boxplot() +
  geom_hline(yintercept=0)





#################### bISON WASTE DECOMPOSABILITY OPTIONS INTO GROUPED BOX PLOTS???

box_plot_f_SSM_L_10 <- c()
box_plot_f_SSM_L_20 <- c()
box_plot_f_SSM_L_30 <- c()
box_plot_f_SSM_L_40 <- c()
box_plot_f_SSM_L_50 <- c()
box_plot_f_SSM_L_60 <- c()
box_plot_f_SSM_L_70 <- c()
box_plot_f_SSM_L_80 <- c()
box_plot_f_SSM_L_90 <- c()
box_plot_f_SSM_L_100 <- c()

box_plot_f_SSM_M_10 <- c()
box_plot_f_SSM_M_20 <- c()
box_plot_f_SSM_M_30 <- c()
box_plot_f_SSM_M_40 <- c()
box_plot_f_SSM_M_50 <- c()
box_plot_f_SSM_M_60 <- c()
box_plot_f_SSM_M_70 <- c()
box_plot_f_SSM_M_80 <- c()
box_plot_f_SSM_M_90 <- c()
box_plot_f_SSM_M_100 <- c()

box_plot_f_SSM_H_10 <- c()
box_plot_f_SSM_H_20 <- c()
box_plot_f_SSM_H_30 <- c()
box_plot_f_SSM_H_40 <- c()
box_plot_f_SSM_H_50 <- c()
box_plot_f_SSM_H_60 <- c()
box_plot_f_SSM_H_70 <- c()
box_plot_f_SSM_H_80 <- c()
box_plot_f_SSM_H_90 <- c()
box_plot_f_SSM_H_100 <- c()


for (val in 1:length(div_dat)){
  if(is.na(foward_dat[,2][[1]][[1]][val]) | foward_dat[,2][[1]][[1]][val] == 0) {next}
  
  box_plot_f_SSM_L_10[val] <- div_dat[, 32][[1]][[1]][val]- div_dat[, 2][[1]][[1]][val]
  box_plot_f_SSM_L_20[val] <- div_dat[, 33][[1]][[1]][val] - div_dat[, 3][[1]][[1]][val]
  box_plot_f_SSM_L_30[val] <- div_dat[, 34][[1]][[1]][val] - div_dat[, 4][[1]][[1]][val]
  box_plot_f_SSM_L_40[val] <- div_dat[, 35][[1]][[1]][val] - div_dat[, 5][[1]][[1]][val]
  box_plot_f_SSM_L_50[val] <- div_dat[, 36][[1]][[1]][val] - div_dat[, 6][[1]][[1]][val]
  box_plot_f_SSM_L_60[val] <- div_dat[, 37][[1]][[1]][val] - div_dat[, 7][[1]][[1]][val]
  box_plot_f_SSM_L_70[val] <- div_dat[, 38][[1]][[1]][val] - div_dat[, 8][[1]][[1]][val]
  box_plot_f_SSM_L_80[val] <- div_dat[, 39][[1]][[1]][val] - div_dat[, 9][[1]][[1]][val]
  box_plot_f_SSM_L_90[val] <- div_dat[, 40][[1]][[1]][val] - div_dat[, 10][[1]][[1]][val]
  box_plot_f_SSM_L_100[val] <- div_dat[, 41][[1]][[1]][val] - div_dat[, 11][[1]][[1]][val]
  
  
  box_plot_f_SSM_M_10[val] <- div_dat[, 42][[1]][[1]][val]- div_dat[, 2][[1]][[1]][val]
  box_plot_f_SSM_M_20[val] <- div_dat[, 43][[1]][[1]][val] - div_dat[, 3][[1]][[1]][val]
  box_plot_f_SSM_M_30[val] <- div_dat[, 44][[1]][[1]][val] - div_dat[, 4][[1]][[1]][val]
  box_plot_f_SSM_M_40[val] <- div_dat[, 45][[1]][[1]][val] - div_dat[, 5][[1]][[1]][val]
  box_plot_f_SSM_M_50[val] <- div_dat[, 46][[1]][[1]][val] - div_dat[, 6][[1]][[1]][val]
  box_plot_f_SSM_M_60[val] <- div_dat[, 47][[1]][[1]][val] - div_dat[, 7][[1]][[1]][val]
  box_plot_f_SSM_M_70[val] <- div_dat[, 48][[1]][[1]][val] - div_dat[, 8][[1]][[1]][val]
  box_plot_f_SSM_M_80[val] <- div_dat[, 49][[1]][[1]][val] - div_dat[, 9][[1]][[1]][val]
  box_plot_f_SSM_M_90[val] <- div_dat[, 50][[1]][[1]][val] - div_dat[, 10][[1]][[1]][val]
  box_plot_f_SSM_M_100[val] <- div_dat[, 51][[1]][[1]][val] - div_dat[, 11][[1]][[1]][val]
  
  box_plot_f_SSM_H_10[val] <- div_dat[, 52][[1]][[1]][val]- div_dat[, 2][[1]][[1]][val]
  box_plot_f_SSM_H_20[val] <- div_dat[, 53][[1]][[1]][val] - div_dat[, 3][[1]][[1]][val]
  box_plot_f_SSM_H_30[val] <- div_dat[, 54][[1]][[1]][val] - div_dat[, 4][[1]][[1]][val]
  box_plot_f_SSM_H_40[val] <- div_dat[, 55][[1]][[1]][val] - div_dat[, 5][[1]][[1]][val]
  box_plot_f_SSM_H_50[val] <- div_dat[, 56][[1]][[1]][val] - div_dat[, 6][[1]][[1]][val]
  box_plot_f_SSM_H_60[val] <- div_dat[, 57][[1]][[1]][val] - div_dat[, 7][[1]][[1]][val]
  box_plot_f_SSM_H_70[val] <- div_dat[, 58][[1]][[1]][val] - div_dat[, 8][[1]][[1]][val]
  box_plot_f_SSM_H_80[val] <- div_dat[, 59][[1]][[1]][val] - div_dat[, 9][[1]][[1]][val]
  box_plot_f_SSM_H_90[val] <- div_dat[, 60][[1]][[1]][val] - div_dat[, 10][[1]][[1]][val]
  box_plot_f_SSM_H_100[val] <- div_dat[, 61][[1]][[1]][val] - div_dat[, 11][[1]][[1]][val]
  
  
  if (val%%100 == 0) {print(val)}
  
  
}


s_b_s_SSM <- data.frame(box_plot_f_SSM_L_10,
                        box_plot_f_SSM_L_20,
                        box_plot_f_SSM_L_30,
                        box_plot_f_SSM_L_40,
                        box_plot_f_SSM_L_50,
                        box_plot_f_SSM_L_60,
                        box_plot_f_SSM_L_70,
                        box_plot_f_SSM_L_80,
                        box_plot_f_SSM_L_90,
                        box_plot_f_SSM_L_100,
                        
                        box_plot_f_SSM_M_10,
                        box_plot_f_SSM_M_20,
                        box_plot_f_SSM_M_30,
                        box_plot_f_SSM_M_40,
                        box_plot_f_SSM_M_50,
                        box_plot_f_SSM_M_60,
                        box_plot_f_SSM_M_70,
                        box_plot_f_SSM_M_80,
                        box_plot_f_SSM_M_90,
                        box_plot_f_SSM_M_100,
                        
                        box_plot_f_SSM_H_10,
                        box_plot_f_SSM_H_20,
                        box_plot_f_SSM_H_30,
                        box_plot_f_SSM_H_40,
                        box_plot_f_SSM_H_50,
                        box_plot_f_SSM_H_60,
                        box_plot_f_SSM_H_70,
                        box_plot_f_SSM_H_80,
                        box_plot_f_SSM_H_90,
                        box_plot_f_SSM_H_100
                        )

### separating into three dataframes, to be merged again melted. 

s_b_s_SSM_L <- data.frame(box_plot_f_SSM_L_10,
                          box_plot_f_SSM_L_20,
                          box_plot_f_SSM_L_30,
                          box_plot_f_SSM_L_40,
                          box_plot_f_SSM_L_50,
                          box_plot_f_SSM_L_60,
                          box_plot_f_SSM_L_70,
                          box_plot_f_SSM_L_80,
                          box_plot_f_SSM_L_90,
                          box_plot_f_SSM_L_100)

melted_L_SSM <- melt(s_b_s_SSM_L)

ggplot(data=melted_L_SSM) +
  geom_boxplot(aes(x=variable, y=value), color="black", fill="red3") + # default colour options can be listed out by calling grDevices::colors()
  geom_hline(yintercept=0, linetype="dashed") +
  scale_fill_viridis(discrete = TRUE, alpha=0.8) +
  geom_jitter(mapping=aes(x=variable, y=value), color="black", size=0.1, alpha=0.1) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=40),
    axis.title.x = element_text(size=25, face="bold"),
    axis.text.x = element_text(size=20),
    
    axis.title.y = element_text(size=25, face="bold"),
    axis.text.y = element_text(size=20),
    
    axis.line.y = element_line(),
    axis.line.x = element_line(),
    
    
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
    
  ) +
  ggtitle("SSM High decomposition") +
  scale_x_discrete(labels=c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) + 
  labs(x = "Years", y="Difference in t C per Ha") +
  coord_cartesian(ylim = c(-0.02, 0))
  
  


s_b_s_SSM_M <- data.frame(box_plot_f_SSM_M_10,
                          box_plot_f_SSM_M_20,
                          box_plot_f_SSM_M_30,
                          box_plot_f_SSM_M_40,
                          box_plot_f_SSM_M_50,
                          box_plot_f_SSM_M_60,
                          box_plot_f_SSM_M_70,
                          box_plot_f_SSM_M_80,
                          box_plot_f_SSM_M_90,
                          box_plot_f_SSM_M_100)

melted_M_SSM <- melt(s_b_s_SSM_M)

ggplot(data=melted_M_SSM) +
  geom_boxplot(aes(x=variable, y=value), color="black", fill="darkorange") + # default colour options can be listed out by calling grDevices::colors()
  geom_hline(yintercept=0, linetype="dashed") +
  scale_fill_viridis(discrete = TRUE, alpha=0.8) +
  geom_jitter(mapping=aes(x=variable, y=value), color="black", size=0.1, alpha=0.1) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=40),
    
    axis.title.x = element_text(size=25, face="bold"),
    axis.text.x = element_text(size=20),
    
    axis.title.y = element_text(size=25, face="bold"),
    axis.text.y = element_text(size=20),
    
    axis.line.y = element_line(),
    axis.line.x = element_line(),
    
    
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  ) +
  ggtitle("SSM Medium decomposition") +
  scale_x_discrete(labels=c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) + 
  labs(x = "Years", y="Difference in t C per Ha per Yr") + 
  scale_y_continuous(labels = label_comma()) +
  coord_cartesian(ylim = c(-0.0002, 0.0002))


s_b_s_SSM_H <- data.frame(box_plot_f_SSM_H_10,
                         box_plot_f_SSM_H_20,
                         box_plot_f_SSM_H_30,
                         box_plot_f_SSM_H_40,
                         box_plot_f_SSM_H_50,
                         box_plot_f_SSM_H_60,
                         box_plot_f_SSM_H_70,
                         box_plot_f_SSM_H_80,
                         box_plot_f_SSM_H_90,
                         box_plot_f_SSM_H_100)

melted_H_SSM <- melt(s_b_s_SSM_H)

ggplot(data=melted_H_SSM) +
  geom_boxplot(aes(x=variable, y=value), color="black", fill="forestgreen") + # default colour options can be listed out by calling grDevices::colors()
  geom_hline(yintercept=0, linetype="dashed") +
  scale_fill_viridis(discrete = TRUE, alpha=0.8) +
  geom_jitter(mapping=aes(x=variable, y=value), color="black", size=0.1, alpha=0.1) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=40),
    
    axis.title.x = element_text(size=25, face="bold"),
    axis.text.x = element_text(size=20),
    
    axis.title.y = element_text(size=25, face="bold", margin=margin(10, 10, 10, 10)),
    axis.text.y = element_text(size=20),
    
    axis.line.y = element_line(),
    axis.line.x = element_line(),
    
    
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
    
  ) +
  ggtitle("SSM Low decomposition") +
  scale_x_discrete(labels=c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) + 
  labs(x = "Years", y="Difference in t C per Ha per") + 
  coord_cartesian(ylim = c(0, 0.02))
  


#s_b_s_SSM_H < data.frame(box_plot_f_SSM_H_10,
 #                        box_plot_f_SSM_H_20,
  #                       box_plot_f_SSM_H_30,
   #                      box_plot_f_SSM_H_40,
    #                     box_plot_f_SSM_H_50,
     #                    box_plot_f_SSM_H_60,
      #                   box_plot_f_SSM_H_70,
       #                  box_plot_f_SSM_H_80,
          #               box_plot_f_SSM_H_90,
        ##                 box_plot_f_SSM_H_100)
#

###################### Difference between H and L:

max_Min_Diff <- melted_H_SSM

max_Min_Diff[2] <- melted_H_SSM[2] + melted_L_SSM[2]

ggplot(data=max_Min_Diff) +
  geom_boxplot(aes(x=variable, y=value), color="black", fill="forestgreen") + # default colour options can be listed out by calling grDevices::colors()
  geom_hline(yintercept=0) +
  scale_fill_viridis(discrete = TRUE, alpha=0.8) +
  geom_jitter(mapping=aes(x=variable, y=value), color="black", size=0.1, alpha=0.1) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=40)
  ) +
  ggtitle("Difference between low and high decomposition") +
  scale_x_discrete(labels=c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)) + 
  labs(x = "Years", y="t C per Ha per Yr") +
  scale_y_continuous(labels = label_comma())




#melted_total <- melt(c(s_b_s_SSM_L, s_b_s_SSM_M, s_b_s_SSM_H), id.vars= c(s_b_s_SSM_L, s_b_s_SSM_M, s_b_s_SSM_H))




#### Check for difference between



names <- c('box_plot_f_SSM_L_10',
           'box_plot_f_SSM_L_20',
           'box_plot_f_SSM_L_30',
           'box_plot_f_SSM_L_40',
           'box_plot_f_SSM_L_50',
           'box_plot_f_SSM_L_60',
           'box_plot_f_SSM_L_70',
           'box_plot_f_SSM_L_80',
           'box_plot_f_SSM_L_90',
           'box_plot_f_SSM_L_100',
           
          'box_plot_f_SSM_M_10',
          'box_plot_f_SSM_M_20',
          'box_plot_f_SSM_M_30',
          'box_plot_f_SSM_M_40',
          'box_plot_f_SSM_M_50',
          'box_plot_f_SSM_M_60',
          'box_plot_f_SSM_M_70',
          'box_plot_f_SSM_M_80',
          'box_plot_f_SSM_M_90',
          'box_plot_f_SSM_M_100',
           
           'box_plot_f_SSM_H_10',
           'box_plot_f_SSM_H_20',
           'box_plot_f_SSM_H_30',
           'box_plot_f_SSM_H_40',
           'box_plot_f_SSM_H_50',
           'box_plot_f_SSM_H_60',
           'box_plot_f_SSM_H_70',
           'box_plot_f_SSM_H_80',
           'box_plot_f_SSM_H_90',
           'box_plot_f_SSM_H_100')


melted_SSM <- melt(s_b_s_SSM, id.vars=names)

ggplot(data=melted_SSM, aes(x=variable, y=value, fill)) + 
  geom_boxplot() +
  geom_hline(yintercept=0)

ggplot(data=melted_SSM) +
  geom_boxplot(aes(x=variable, y=value, color=variable, fill=variable)) +
  geom_hline(yintercept=0) +
  scale_fill_viridis(discrete = TRUE, alpha=0.8) +
  geom_jitter(mapping=aes(x=variable, y=value), color="black", size=0.1, alpha=0.1) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("A boxplot with jitter") +
  xlab("")

##### modified to be in different facets #############

s_b_s_SSM_L <- data.frame(c(box_plot_f_SSM_L_10,
                          box_plot_f_SSM_L_20,
                          box_plot_f_SSM_L_30,
                          box_plot_f_SSM_L_40,
                          box_plot_f_SSM_L_50,
                          box_plot_f_SSM_L_60,
                          box_plot_f_SSM_L_70,
                          box_plot_f_SSM_L_80,
                          box_plot_f_SSM_L_90,
                          box_plot_f_SSM_L_100,
                          
                          box_plot_f_SSM_M_10,
                          box_plot_f_SSM_M_20,
                          box_plot_f_SSM_M_30,
                          box_plot_f_SSM_M_40,
                          box_plot_f_SSM_M_50,
                          box_plot_f_SSM_M_60,
                          box_plot_f_SSM_M_70,
                          box_plot_f_SSM_M_80,
                          box_plot_f_SSM_M_90,
                          box_plot_f_SSM_M_100,
                          
                          box_plot_f_SSM_H_10,
                          box_plot_f_SSM_H_20,
                          box_plot_f_SSM_H_30,
                          box_plot_f_SSM_H_40,
                          box_plot_f_SSM_H_50,
                          box_plot_f_SSM_H_60,
                          box_plot_f_SSM_H_70,
                          box_plot_f_SSM_H_80,
                          box_plot_f_SSM_H_90,
                          box_plot_f_SSM_H_100), 
)






ggplot(data=melted_SSM, aes(x=variable, y=value, fill)) + 
  geom_boxplot() +
  geom_hline(yintercept=0)


########################################################################

SOC_BAU_2020 <- c()
SOC_BAU_2120 <- c()
SOC_SSM_LOW <- c()
SOC_SSM_MED <- c()
SOC_SSM_HIGH <- c()

for (val in 1:length(div_dat)){
  if(is.na(foward_dat[,2][[1]][[1]][val]) | foward_dat[,2][[1]][[1]][val] == 0) {next}
  SOC_BAU_2020[val] <- foward_dat[,2][[1]][[1]][val]
  SOC_BAU_2120[val] <- foward_dat[,3][[1]][[1]][val]
  
  SOC_SSM_LOW[val] <- foward_dat[,10][[1]][[1]][val]
  SOC_SSM_MED[val] <- foward_dat[,11][[1]][[1]][val]
  SOC_SSM_HIGH[val] <- foward_dat[,12][[1]][[1]][val]
  if (val%%100 == 0) {print(val)}
}


final_box_data <- data.frame(SOC_BAU_2120, SOC_SSM_LOW, SOC_SSM_MED, SOC_SSM_HIGH)
names(final_box_data)[1] <- "SOC_BAU_2120"
names(final_box_data)[2] <- "SOC_SSM_LOW"
names(final_box_data)[3] <- "SOC_SSM_MED"
names(final_box_data)[4] <- "SOC_SSM_HIGH"

boxplot(final_box_data)
## fun fact: apparently lowest value in BAU is 8.56
names = data.frame('SOC_BAU_2120', 'SOC_SSM_LOW', 'SOC_SSM_MED', 'SOC_SSM_HIGH')


final_mod <- melt(final_box_data)


ggplot(data=final_mod) +
  geom_boxplot(aes(x=variable, y=value, color=variable, fill=variable)) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(mapping=aes(x=variable, y=value), color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("A boxplot with jitter") +
  xlab("")

################### EXPONENTIATED DATA ###################

expo_data <- final_box_data^3

final_mod2 <- melt(expo_data)

ggplot(data=final_mod2) +
  geom_boxplot(aes(x=variable, y=value, color=variable, fill=variable)) +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(mapping=aes(x=variable, y=value), color="black", size=0.4, alpha=0.9) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("A boxplot with jitter") +
  xlab("")


#############################

#### graphing uncertainty analysis:
Max_BAU <- c()
Min_BAU <- c()
Med_BAU <- c()

Max_SSM <- c()
Med_SSM <- c()
Min_SSM <- c()


SUM_Max_BAU <- 0
SUM_Min_BAU <- 0
SUM_Med_BAU <- 0
SUM_Max_SSM <- 0
SUM_Med_SSM <- 0
SUM_Min_SSM <- 0




for (val in 1:length(div_dat)){
  if(is.na(foward_dat[,2][[1]][[1]][val]) | foward_dat[,2][[1]][[1]][val] == 0) {next}
  
  Max_BAU[val] <- foward_dat[, 14][[1]][[1]][val]
  Min_BAU[val] <- foward_dat[, 13][[1]][[1]][val]
  Med_BAU[val] <- foward_dat[, 3][[1]][[1]][val]
  Max_SSM[val] <- foward_dat[, 16][[1]][[1]][val]
  Med_SSM[val] <- foward_dat[, 11][[1]][[1]][val]
  Min_SSM[val] <- foward_dat[, 15][[1]][[1]][val]
  
  
  # now calculate the difference between everything and such without bison...
  SUM_Max_BAU <- SUM_Max_BAU + (Max_BAU[val] - Med_BAU[val])
  SUM_Min_BAU <- SUM_Min_BAU + (Min_BAU[val] - Med_BAU[val])
  SUM_Max_SSM <- SUM_Max_SSM + (Max_SSM[val] - Med_BAU[val])
  SUM_Med_SSM <- SUM_Med_SSM + (Med_SSM[val] - Med_BAU[val])
  SUM_Min_SSM <- SUM_Min_SSM + (Min_SSM[val] - Med_BAU[val])
  
  if (val%%100 == 0) {print(val)}
}

sum_Diff <- data.frame(SUM_Max_BAU,
                       
                       SUM_Max_SSM,
                       SUM_Min_BAU,
                       #SUM_Med_SSM,
                       SUM_Min_SSM
                       #Max_BAU,
                       #Max_SSM,
                       #Med_BAU,
                       #Med_SSM,
                       #Min_BAU,
                       #Min_SSM
                       )
diff_Dat <- data.frame(Max_BAU,
                       Max_SSM,
                       Min_BAU,
                       Min_SSM)


melted_sum <- melt(sum_Diff)

ggplot(melted_sum, aes(x=variable, y=value)) +
  geom_bar(stat="identity")


melted_diff <- melt(diff_Dat)

ggplot(melted_diff) + 
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
  
### calculate mean and standard deviation for values:
avg_val = mean(Med_BAU[!is.na(Med_BAU)])
print(avg_val)
N = length(Med_BAU[!is.na(Med_BAU)])
deviations <- Med_BAU[!is.na(Med_BAU)] - mean(Min_BAU[!is.na(Med_BAU)])
s <- deviations^2 # the "s" in r.m.s.
m <- sum(s)/N # the mean; the "m" in r.m.s.
sd <- sqrt(m) # the sqare root, the "r" in r.m.s.
print(sd)

#get_Sum <- sum(Min_BAU[!is.na(Min_BAU)])
#print(get_Sum)


