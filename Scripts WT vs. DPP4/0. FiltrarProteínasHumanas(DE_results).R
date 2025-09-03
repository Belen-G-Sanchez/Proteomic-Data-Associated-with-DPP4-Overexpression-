#Choose directory where to find the files

#Packages needed to run the R: To read excel files install readxl and to work with FASTA files install Biostrings.
install.packages("readxl")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")
install.packages("writexl")
install.packages("stringi")
install.packages("stringr")
install.packages("dplyr")
library(readxl)
library(Biostrings)
library(openxlsx)
library(dplyr)
library(writexl)

#Read excel and select the column we want to focus on ‘Protein.Group’.
archivo_excel <- "DE_results_WTvsDPP4.xls"
datos_excel <- read_excel(archivo_excel)
ids_excel <- datos_excel$Protein.Group

#Read FASTA 
archivo_fasta <- "HumanProteomeUP000005640_9606.fasta"
#For proteins
secuencias <- readAAStringSet(archivo_fasta)
#To extract the sequence IDs from FASTA sequences
ids_fasta <- names(secuencias)
#To extract the gene name, the gene is found before _HUMAN
ids_fasta_clean <- sapply(strsplit(ids_fasta, "_HUMAN"), `[`,1)
ids_fasta_protein_clean <- sapply(strsplit(ids_fasta, "_HUMAN "), `[`,2)
ids_fasta_protein_cleanest <- sapply(strsplit(ids_fasta_protein_clean, " OS=Homo"), `[`,1)

#To extract the protein name, the protein is found after _HUMAN
protein_names <- sapply(strsplit(ids_fasta, "_HUMAN"), `[`,1)
protein_names_clean <- str_remove(protein_names, "sp|")
protein_names_cleaned <- str_remove(protein_names_clean, "tr")
protein_names_cleanest <- gsub("\\|", " ", protein_names_cleaned)
protein_names_cleanesttwo <- sapply(strsplit(protein_names_cleanest, " "), `[`, 3) #to remove the gene name (before the 3rd space)

#To remove what we are not interested in from the name 
#First, install str function
gene_name <- str_remove(ids_fasta_clean, "sp|") # to remove ‘sp’ before the name
gene_name_clean <- gsub("\\|", " ", gene_name) #to remove ‘|’ before the name 
gene_name_cleaned <- str_remove(gene_name_clean, "tr") #to remove ‘’tr’ before the name
gene_name_cleanest <- sapply(str_split(gene_name_cleaned, " "), `[`, 2) #to remove the gene name (after the 2nd space)

#Link the Excel column with the Fasta file to find the names of common proteins between our data and the human species.
ids_comunes <- intersect(ids_excel, gene_name_cleanest)

#Save the cleaned fasta file in Excel
df2 <- data.frame(Protein.Group = gene_name_cleanest, stringsAsFactors = FALSE) #name of the column of interest from FASTA
write_xlsx(df2, "gene_name_cleanest.xlsx")
#Join uniprot ID with gene name
protein_gene <- data.frame(Protein.Group = gene_name_cleanest, Protein.Names = ids_fasta_protein_cleanest)
write_xlsx(protein_gene, "proteinname_geneID.xlsx")
# Join protein name with ID of uniprot
protein_names_gene <- data.frame(Protein.Names = protein_names_cleanesttwo, Protein.Group = gene_name_cleanest)
write_xlsx(protein_names_gene, "proteinname_geneID_proteinname.xlsx")

#Read first excel file of results as data frame 
df1 <- read_excel(archivo_excel)

#Merge files keeping only rows with IDs in common
df_comun <- inner_join(df1, df2, by = "Protein.Group")
write_xlsx(df_comun, "datos_filtrados.xlsx")

gene_protein_comun <- df1 %>%
  left_join(protein_gene, by = "Protein.Group")  # Join for the column "Protein.Group"
write_xlsx(gene_protein_comun, "Excel_proteinname_ID.xlsx")

protein_protein_comun <- df1 %>%
  left_join(protein_names_gene, by = "Protein.Group")  # Join for the column "Protein.Group"
write_xlsx(protein_protein_comun, "Excel_proteinname_ID_ahora.xlsx")

results_without_contaminating <- inner_join(gene_protein_comun, df2, by = "Protein.Group")
write_xlsx(results_without_contaminating, "Excel_proteinname_ID_nocontaminating.xlsx")

results_without_contaminating_completo <- inner_join(protein_protein_comun, protein_gene, by = "Protein.Group")
write_xlsx(results_without_contaminating_completo, "DE_results_human.xlsx")
