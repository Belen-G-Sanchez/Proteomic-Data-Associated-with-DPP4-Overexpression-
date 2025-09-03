# -----------------------------
# 1. Install and load the necessary libraries
# -----------------------------
install.packages(c("pheatmap", "dplyr", "tidyverse"))
library(pheatmap)
library(dplyr)
library(tidyverse)

# -----------------------------
# 2. Upload data
# -----------------------------
# Upload the protein intensities file
# Adjust the separator and file name accordingly.
intensities <- read.csv("report.pg_matrix_WTvsDPP4_HumanProtein_Heatmap.csv", 
                        sep = ";", 
                        header = TRUE, 
                        stringsAsFactors = FALSE) 

# Load the differentially expressed proteins file
# It is assumed that this file contains at least the column "Protein_ID"
diff_expr <- read.csv("DE_results_human_WTvsDPP4.csv", 
                      sep = ";", 
                      header = TRUE, 
                      stringsAsFactors = FALSE)

report_pg_filtered <- diff_expr %>% 
  filter(abs(log2_FC) > 1.2)

# -----------------------------
# 3. Review and prepare identifiers
# -----------------------------
# Print the column names to verify that you have the correct identifiers
print("Columns en intensities:")
print(colnames(intensities))
print("Columns en diff_expr:")
print(colnames(diff_expr))
# **Make sure that in both files the identifier is called ‘Protein_ID ’**.

# If necessary, we make the identifiers in the intensity file unique
intensities$Protein.Group <- make.unique(as.character(intensities$Protein.Group))

# -----------------------------
# 4. # Filter to retain only differentially expressed proteins. 
# -----------------------------
# Filter intensities to keep only the proteins that are in diff_expr**
intensities_DE <- intensities %>% filter(Protein.Group %in% report_pg_filtered$Protein_ID) 

cat("Proteins in intensities:", nrow(intensities), "\n")
cat("Proteins in report_pg_filtered:", nrow(report_pg_filtered), "\n")
cat("Proteins filtered (differentially expressed):", nrow(intensities_DE), "\n")

# -----------------------------
# 5. Convert to matrix
# -----------------------------
# **Convert column identifiers to row names** 
intensities_DE <- intensities_DE %>% column_to_rownames("Protein.Group")

# **Convert the data frame to a matrix** 
heatmap_matrix <- as.matrix(intensities_DE)

# Force all data to be numeric
heatmap_matrix <- apply(heatmap_matrix, 2, as.numeric)
rownames(heatmap_matrix) <- rownames(intensities_DE)

# -----------------------------
# 6. (Optional) Normalise data by rows
# -----------------------------

# Replace infinite values and delete rows with NA
heatmap_matrix[is.infinite(heatmap_matrix)] <- NA
heatmap_matrix <- heatmap_matrix[complete.cases(heatmap_matrix), ]

# **Normalise by ranks (Z-score), to compare the relative expression pattern**
heatmap_matrix_norm <- t(scale(t(heatmap_matrix)))

# -----------------------------
# 7. Generate the heatmap
# -----------------------------
# Define the colour palette
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)

# Adjust chart size
png("heatmap.png", width = 3000, height = 2000, res = 300)  # Ajusta dimensiones y resolución
pheatmap(heatmap_matrix_norm, 
         color = heatmap_colors,
         scale = "none",                        
         clustering_distance_rows = "euclidean",  
         clustering_distance_cols = "euclidean",  
         clustering_method = "complete",        
         fontsize_row = 3, # Adjust the size of the protein names (rows)**
         fontsize_col = 4, # Adjust the size of the column names*.
         fontsize = 5,     # Reduction of font size
         cellheight = 2.75, # Reduce the height of cells to improve visibility
         cellwidth = 10,
         angle_col = 45,
         treeheight_row = 5, # Reduce dendrogram height of rows
         treeheight_col = 5,     
         main = "Differentially Expressed Proteins") 
dev.off()  # Closes the graphical device
