# Load the necessary libraries
library(ggplot2)
library(dplyr)
library(ggrepel)

# Upload data
report.pg_matrix <- read.csv("DE_results_human_WTvsDPP4.csv", sep = ";", header = TRUE, dec = ".", check.names = FALSE)
experiment_annotation <- read.csv("experiment_annotation_WTvsDPP4.csv", sep = ";", header = TRUE, dec = ".", check.names = FALSE)

# Filtering and processing relevant data
data_volcano <- report.pg_matrix %>%
  select(Protein_ID, log2_FC, p_adj) %>%
  filter(!is.na(log2_FC) & !is.na(p_adj)) %>% # Delete NA
  mutate(
    log_p_value = -log10(p_adj),
    significance = ifelse(p_adj < 0.05 & abs(log2_FC) > 1, "Significant", "Not significant")
  )

# Create Volcano Plot
png("VolcanoPlot.png", width = 2000, height = 3000, res = 300)  # Export with good resolution
ggplot(data_volcano, aes(x = log2_FC, y = log_p_value)) +
  
  # Change colour of dots according to significance
  geom_point(aes(color = significance), size = 2, alpha = 0.8) +
  
  # Define custom colours
  scale_color_manual(values = c("Significant" = "black", "Not significant" = "grey")) +
  
  # Líneas de corte
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +  
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  
  
  # Añadir etiquetas para las proteínas más significativas
  geom_label_repel(
    data = subset(data_volcano, p_adj < 0.05 & abs(log2_FC) > 1.2),  # Significant protein only
    aes(label = Protein_ID),
    size = 5, 
    fill = "white",        # White background
    color = "black",       # Black text
    label.size = 0.5,      # Rectangle border thickness
    box.padding = 0.5,     # Spacing around text
    max.overlaps = 15      # Avoid excessive overlaps
  ) +
  
  # Customisation of the chart
  labs(title = "Volcano Plot - Differentially expressed proteins",
       x = "Log2 Fold Change", y = "-Log10(p-value)") +
  
  theme_minimal() +   # Apply theme minimal
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 14),  # X-axis title size
    axis.title.y = element_text(size = 14),  # Y-axis title size
    axis.text.x = element_text(size = 10),   # X-axis label size
    axis.text.y = element_text(size = 10)    # Y-axis label size
  )

dev.off()
