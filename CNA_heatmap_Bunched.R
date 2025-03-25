library(ggplot2)
library(readxl)
library(dplyr)
library(patchwork)
setwd("/home/carina/Documents/CNV_analysis/plots")


list_path <- c("m1" = "/home/carina/Documents/CNV_analysis/TAPS_summary/model1/model1_static.xlsx", 
               "m2" = "/home/carina/Documents/CNV_analysis/TAPS_summary/Model2/model2_all_query_static.xlsx", 
               "m3" = "/home/carina/Documents/CNV_analysis/TAPS_summary/Model3/model3_all_query_static.xlsx",
               "m4" = "/home/carina/Documents/CNV_analysis/TAPS_summary/Model4/Model4_query.xlsm",
               "m5" ="/home/carina/Documents/CNV_analysis/TAPS_summary/Model5/Model5_query.xlsm",
               "m6" = "/home/carina/Documents/CNV_analysis/TAPS_summary/Model6/model6_query.xlsx",
               "m7" = "/home/carina/Documents/CNV_analysis/TAPS_summary/Model7/model7_query.xlsm"
               
)


list_samples = list(
  m1 = c("1xD","1xD_rec", "1xP1_1E",  "1xP2c_1D" ,"1P5_2E", "1P3cFP_1", "1P5FP_3",   "1P4cKd_1", "1P5Kd_1"),
  m2 = c("2D", "2P3_1D", "2P3_rec", "2xP3_5D", "2xP6_2D",  "2P8_2D", "2P4FP_1", "2P6FP_2", "2P4Kd_2", 
         "2P6Kd_1", "2P6_3D_1IC-Lv", "2P6_3D_1IC-Ov", "2P6_3D_1IC-SR", "2P6_3D_2IC-Lv", "2P6_3D_2IC-SR", "2P6_3D_3IC-Lv"),
  m3 = c("3D", "3D_rec","3P3_1D", "3P4c_2D", "3P4c_3E", "3P10_1D","3P5cFP_1", "3P7FP_2", "3P5cKd_2", 
         "3P7Kd_3", "3P4IC_3-Lg","3P5IC_1-SR", "3P5IC_2-Lg", "3P5IC_3-Lg", "3P5TV_1-Lg", "3P5TV_2-Lg", "3P10TV_1", "3P10TV_2"),
  m4 = c("4D", "4P3_1E", "4P4SC_1E", "4xP6_1E", "4P7c_1D","4P4FP_2", "4P6FP_1", "4P4Kd_2", "4P6Kd_2"), #"4P6_1E", ,"4P2C_1E"
  
  m5 = c("5D", "5P1", "5xP3_1D","5P5_1E",  "5P6_1D", "5P3FP_1","5xP5FP_1",  "5P3Kd_2","5P5Kd_2"),
  m6 = c("6D", "6P2","6P3", "6P5SC_2E", "6P6SC_1E", "6P4FP_1", "6P6FP_1", "6P5Kd_2", "6P7Kd_2",
         "6P2IC_2-SR", "6P2TV_1-SR", "6P3IC_4", "6P3TV_1"),
  m7 =  c("7D", "7P4c_1D",  "7P6_2E", "7P8_3D","7P4cFP_1", "7P6FP_1", "7P5Kd_1", "7P7Kd_1", "7P4IC_2-Lv") # "7P3_1E"     "7P3_3E"  
)


df <- data.frame()

models <-c("m1", "m4", "m5")

for (m in models){
  for (i in 1:length(list_samples[[m]])){
    data <- read_xlsx(list_path[[m]],
                      sheet=list_samples[[m]][i]) %>% 
      mutate(log2=as.numeric(log2), imba=as.numeric(imba),
             value2 = case_when(log2 > limit ~ "Gain",log2 < -limit ~ "Loss",
                                log2 > -(limit-half.limit) & log2 <  (limit-half.limit) & (imba < 0.45 | imba > 0.55) ~ "Copy-Neutral LOH",
                                log2 > -(limit+half.limit) & log2 < -(limit-half.limit) & (imba < 0.45 | imba > 0.55) ~ "Subtle Loss with LOH",
                                TRUE ~ "Balanced"), start2 = Start / 100000, end2 = End / 100000)
    data$sample <- list_samples[[m]][i]
    data$model <- m
    
    df <- bind_rows(df, data)
    
  }
  
}


ggplot(df, aes(x = as.factor(start2), y=sample, fill=value2)) + 
geom_tile() + 
     facet_grid(model~Chromosome, space="free", scale="free") +
  scale_fill_manual(values = custom_colors) + 
  #scale_fill_manual(values = color_palette, name = "", labels = c(0, 1, 2, 3, 4, 5, 6, 7, 8)) +
  theme(axis.line.x=element_blank(),
        axis.text.x=element_blank(),
       axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank(),
      legend.position = "none")

