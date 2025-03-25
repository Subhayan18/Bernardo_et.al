# Load necessary library
library(ggplot2)
library(readxl)
library(dplyr)
library(patchwork)
setwd("C:/Users/Subhayan/Box Sync/Helping_with_analysis/Carina/Chr_summary_plot")
file_path<-c("C:/Users/Subhayan/Box Sync/Helping_with_analysis/Carina/TAPS_summary/model1/model1_static.xlsx")
  
plot.func <- function(data=data, sample_name){
  # Create your dataset
  plot.data <- data %>% select(
    Chr=Chromosome,
    Start = start2,
    End = end2,
    Value = value2
  )
# # Calculate the maximum End for each chromosome
# max_end_per_chr <- plot.data %>%
#   group_by(Chr) %>%
#   summarize(MaxEnd = max(End))
# 
# # Create the background track for each chromosome
# background <- max_end_per_chr %>%
#   mutate(Start = 0, Value = "Balanced") %>%
#   select(Chr, Start, End = MaxEnd, Value)

# Combine the background with the original data
plot.data <- bind_rows(background, plot.data)


# Define the desired order of chromosomes
chr_order <- c(paste0("chr", 1:22), "chrX", "chrY")

# Convert the Chr column to a factor with the specified order
plot.data$Chr <- factor(plot.data$Chr, levels = chr_order)

# Define your custom colors
custom_colors <- c("Balanced" = "white", 
                   "Gain" = "dodgerblue", 
                   "Loss" = "darkorange2",
                   "Subtle Loss with LOH" = "lightpink",
                   "Copy-Neutral LOH" = "cornsilk"
                   )

# Create the plot
ggplot(plot.data, aes(xmin = Start, xmax = End, ymin = 0, ymax = 1, fill = Value)) +
  geom_rect() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_minimal() +  # Use minimal theme
  theme(
    legend.position = "none",  # Removes the legend
    strip.text = element_blank(),  # Removes facet labels
    axis.title.x = element_blank(),  # Removes x-axis title
    axis.text.x = element_blank(),  # Removes x-axis text
    axis.title.y = element_text(size = 8),  # Removes y-axis title
    axis.text.y = element_blank(),  # Removes y-axis text
    axis.ticks.y = element_blank(),  # Removes y-axis ticks
    panel.grid = element_blank(),  # Removes grid lines
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0),  # Removes plot margins
    panel.spacing = unit(0, "lines")  # Removes spacing between panels
  ) +
  scale_fill_manual(values = custom_colors) +  # Apply custom colors 
  facet_grid(. ~ Chr, space = "free_x", scales = "free_x") +
  labs(y = sample_name) + annotation_custom(
    grob = grid::rectGrob(gp = grid::gpar(col = "black", fill = NA, lwd = 1)),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) # Add the sample name as y-axis label 

}
background<-read.csv("background.csv")

sample_names <- excel_sheets(file_path)
limit=0.15
half.limit=0.1

plot_list <- list()
for (i in 1:length(sample_names)){
data <- read_xlsx(file_path,
                  sheet=sample_names[i]) %>% 
  mutate(log2=as.numeric(log2), imba=as.numeric(imba),
    value2 = case_when(log2 > limit ~ "Gain",log2 < -limit ~ "Loss",
    log2 > -(limit-half.limit) & log2 <  (limit-half.limit) & (imba < 0.45 | imba > 0.55) ~ "Copy-Neutral LOH",
    log2 > -(limit+half.limit) & log2 < -(limit-half.limit) & (imba < 0.45 | imba > 0.55) ~ "Subtle Loss with LOH",
    TRUE ~ "Balanced"), start2 = Start / 100000, end2 = End / 100000)
plot_list[[i]] <- plot.func(data, sample_names[i])
}

#case_when(log2 > limit ~ "Gain", log2 < -limit ~ "Loss", TRUE ~ "Balanced")

# Combine plots using patchwork
combined_plot <- wrap_plots(plot_list, ncol = 1,widths=2)

ggsave("model1_combined_limit_0.15_new.pdf", plot = combined_plot, width = 32, height = 18, units = "in", dpi = 600)

# data <- read_xlsx("C:/Users/Subhayan/Box Sync/Helping_with_analysis/Carina/TAPS_summary/model1/model1_static.xlsx",
#                   sheet=sample_names[4])
# data <- data %>% filter(Chromosome=="chr3")
