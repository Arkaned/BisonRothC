# Sensitivity plot 

library(terra)

library(grid)
library(tidyverse)
library(shadowtext)
library(extrafont)

sens_data <- vect("H:/BisonRothC/Statistics/Sensitivity_analysis_F_run/Sensitivity_analysis_F_run.shp")





sens_Data <- (terra::as.data.frame(sens_data))
sens_Data <- na.omit(sens_Data)

names_List <- c()


for (i in 3:length(names(sens_Data))) {
  names_List[(i-2)] <- names(sens_Data)[i]
}

####### % calculations
values_perc <- c()
values_C <- c()

for (i in 3:length(sens_Data)) {
  values_perc[(i-2)] <- abs(1-(sens_Data[1, i]/sens_Data[1, 2]))
}

for (i in 3:length(sens_Data)) {
  values_C[(i-2)] <- abs(sens_Data[1, i] - sens_Data[1, 2])
}

plotting_data <- data.frame(
  count = values_perc,
  name = factor(names_List, levels=names_List),
  y = seq(length(names_List)) * 0.9
)


plotting_C <- data.frame(
  count = values_C,
  name = factor(names_List, levels=names_List),
  y = seq(length(names_List)) * 0.5
)

cleaned_insignificant <- plotting_C
cleaned_insignificant <- cleaned_insignificant[-c(1, 2, 3, 4, 7, 8, 21, 22, 25, 26, 27, 28), ]

colors_List <- c(1, 2, 1, 2, 1, 2, 3, 4, 3, 4, 3, 4, 5, 6, 7, 8)
color_vals <- c("#933705", "#7D3006", "#933705", "#7D3006", "#933705", "#7D3006",
                "#0882A6", "#096F8D", "#0882A6", "#096F8D", "#0882A6", "#096F8D",
                "#278D09", "#257F0A", 
                "#7F220A", "#91280D")
cleaned_insignificant[, 4] <- colors_List
cleaned_insignificant[, 5] <- color_vals
colnames(cleaned_insignificant)[4] <- "COLORS"
colnames(cleaned_insignificant)[5] <- "COLOR_VALS"


plt <- ggplot(cleaned_insignificant) +
  geom_col(aes(count, name), width = 0.8, fill=cleaned_insignificant$COLOR_VALS) +
  scale_y_discrete(name="Inputs +/- 10%") + 
  scale_x_continuous(name="SOC (t C per Ha)", expand=c(0,0), limits = c(0, 4)) +
  theme(axis.title=element_text(face="bold", size=20),
        axis.text.x = element_text(size=15, angle=0),
        axis.text.y = element_text(size=15, angle=0),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
  
  



# text = element_text(family = "Calibri"))












significant_vals <- 

clean_data <- plotting_data


plt <- ggplot(plotting_data) +
  geom_col(aes(count, names_List), width = 0.6) 




plt3 <- plt2 + 
  geom_shadowtext(
    data = subset(data, count < 8),
    aes(count, y = name, label = name),
    hjust = 0,
    nudge_x = 0.3,
    colour = "blue",
    bg.colour = "white",
    bg.r = 0.2,
    family = "Econ Sans Cnd",
    size = 7
  ) + 
  geom_text(
    data = subset(data, count >= 8),
    aes(0, y = name, label = name),
    hjust = 0,
    nudge_x = 0.3,
    colour = "white",
    family = "Econ Sans Cnd",
    size = 7
  )
