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
  filter(!is.na(log2_FC) & !is.na(p_adj)) %>% # Remove NA
  mutate(
    log_p_value = -log10(p_adj),
    significance = case_when(
      p_adj < 0.05 & abs(log2_FC) >= 1 ~ "Significant",
      TRUE ~ "Not significant")
  )

# Featured proteins
proteins_of_interest <- c("P12830", "O00622", "P56199", "P17301", "P26006", "P13612", "P08648", "P23229", "O95471", "Q13835")

# Create Volcano Plot
png("VolcanoPlot.png", width = 3000, height = 2000, res = 300)  # Export with good resolution
ggplot(data_volcano, aes(x = log2_FC, y = log_p_value)) +
  
  # Dots coloured according to condition
  geom_point(aes(color = ifelse(Protein_ID %in% proteins_of_interest, "Marked", significance)), size = 2, alpha = 0.8) +
  
  # Define custom colours
  scale_color_manual(values = c("Significant" = "black", "Not significant" = "grey", "Marked" = "red")) +
  
  # Cutting lines
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +  
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  
  
  # Featured protein labels
  geom_label_repel(data = subset(data_volcano, Protein_ID %in% proteins_of_interest), 
                   aes(label = Protein_ID), 
                   size = 5, fill = "white", color = "red", label.size = 0.5, 
                   box.padding = 0.5, max.overlaps = 10) +  
  
  # Customisation of the chart
  labs(title = "Volcano Plot - Differentially expressed proteins",
       x = "Log2 Fold Change", y = "-Log10(p-value)") +
  
  theme_minimal() +   # theme_minimal is applied correctly
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 18),  # X-axis title size
    axis.title.y = element_text(size = 18),  # Y-axis title size
    axis.text.x = element_text(size = 14),   # X-axis label size
    axis.text.y = element_text(size = 14)    # Y-axis label size
  )

dev.off()
