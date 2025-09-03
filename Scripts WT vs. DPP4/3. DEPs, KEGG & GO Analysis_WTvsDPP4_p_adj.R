install.packages("dplyr") 
library(dplyr)

# Upload the .csv file
report.pg_matrix <- read.csv ("DE_results_human_WTvsDPP4.csv",
                              sep = ";", 
                              header = TRUE, 
                              dec = ".", 
                              check.names = FALSE) 
# Filtering differentially expressed proteins
filtered_report.pg <- report.pg_matrix %>%
  filter(abs(log2_FC) > 1 & p_adj < 0.05)
# Save the filtered file as a new CSV file
write.csv(filtered_report.pg, "report.pg_matrix_filtering.csv", row.names = FALSE)

# Getting protein IDs for further enrichment of pathways
protein_ids <- filtered_report.pg$Protein_ID

# Install "clusterProfiler" package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler", force = TRUE)
library(clusterProfiler)

# Install org.Hs.eg.db package
BiocManager::install("org.Hs.eg.db", force = TRUE)
library(org.Hs.eg.db) # Database for humans

# Convert genes to Entrez ID
entrez_ids <- bitr(protein_ids, fromType = "UNIPROT", 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)

# Remove duplicates, keeping only the first mapping
entrez_ids_unique <- entrez_ids[!duplicated(entrez_ids$UNIPROT), ]

# Enrichment analysis in the KEGG database
enrich_kegg <- enrichKEGG(gene = entrez_ids$ENTREZID, 
                          organism = "hsa", 
                          pvalueCutoff = 0.05)

# View results
head(enrich_kegg)
dotplot(enrich_kegg, showCategory = 15) # Mostrar 10 principales categorías

# Enrichment of pathways with GO

# Install necessary packages
install.packages("BiocManager")  # If you do not have Bioconductor installed
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "ggplot2", force = TRUE))
library(clusterProfiler)  # For enrichment analysis
library(org.Hs.eg.db)     # Human gene annotation database
library(enrichplot)       # For display
library(ggplot2)          # Additional graphics


# Gene list in ENTREZ format
genes_entrez <- entrez_ids$ENTREZID  

# Perform GO (Biological Process) enrichment analysis
ego_BP <- enrichGO(gene = genes_entrez, 
                   OrgDb = org.Hs.eg.db, 
                   keyType = "ENTREZID",
                   ont = "BP",  # "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
                   pAdjustMethod = "BH", 
                   pvalueCutoff = 0.05, 
                   qvalueCutoff = 0.05, 
                   readable = TRUE)  # Converts ENTREZ to readable gene names

#Dot plot
dotplot(ego_BP, showCategory = 15) + 
ggtitle("Top 15 GO Biological Process terms")

#Bar plot 
barplot(ego_BP, showCategory = 15, title = "GO Enrichment Analysis")

#Network diagram (to see relationships between enriched terms)
cnetplot(ego_BP, categorySize = "pvalue", foldChange = NULL)

# Perform GO (Biological Process) enrichment analysis
ego_MF <- enrichGO(gene = genes_entrez, 
                   OrgDb = org.Hs.eg.db, 
                   keyType = "ENTREZID",
                   ont = "MF",  # "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
                   pAdjustMethod = "BH", 
                   pvalueCutoff = 0.05, 
                   qvalueCutoff = 0.05, 
                   readable = TRUE)  # Converts ENTREZ to readable gene names

#Dot plot
dotplot(ego_MF, showCategory = 15) + 
  ggtitle("Top 15 GO Molecular Function terms")

# Perform GO (Cellular Compoment) enrichment analysis
ego_CC <- enrichGO(gene = genes_entrez, 
                   OrgDb = org.Hs.eg.db, 
                   keyType = "ENTREZID",
                   ont = "CC",  # "BP" for Biological Process, "MF" for Molecular Function, "CC" for Cellular Component
                   pAdjustMethod = "BH", 
                   pvalueCutoff = 0.05, 
                   qvalueCutoff = 0.05, 
                   readable = TRUE)  # Converts ENTREZ to readable gene names

#Dot plot
dotplot(ego_CC, showCategory = 15) + 
  ggtitle("Top 15 GO Cellular Component terms")
