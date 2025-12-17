# =============================================================================
# COVID-19 Patient Transcriptomics Data Preparation Script
# =============================================================================
# This script processes patient data, differential gene expression (DGE) results,
# and maps genes to KEGG pathways for downstream analysis.

# Load Required Packages
library(readxl)           # For reading Excel files (.xlsx)
library(clusterProfiler)  # For gene ontology and pathway enrichment analysis
library(org.Hs.eg.db)     # Human genome annotation database (Entrez Gene IDs)

# Set Working Directory
# Note: Update this path to match your project directory
setwd("~/Project2")

####################
### Patient Data ###
####################

# Read patient demographic and clinical data from supplementary material
# Source: mmc3.xlsx (Supplementary Data 3 from publication)
patient_data <- read_xlsx("DATA/mmc3.xlsx")

# The Excel file has headers in the second row, so we use that row as column names
colnames(patient_data) <- patient_data[2,]

# Remove the first two rows (header rows) and keep only the first 4 columns
# This extracts: Patient ID, Age, Sex, and Disease Severity
patient_data <- patient_data[-c(1:2), 1:4]

# Standardize column names for easier reference
colnames(patient_data) <- c("Patient", "Age", "Sex", "Severity")

# Standardize disease severity labels to lowercase for consistency
# Replaces variations like "Mild", "MILD", etc. with "mild"
patient_data$Severity[grepl("Mild", patient_data$Severity)] <- "mild"
patient_data$Severity[grepl("Severe", patient_data$Severity)] <- "severe"
patient_data$Severity[grepl("Moderate", patient_data$Severity)] <- "moderate"


####################
### DGE Analysis ###
####################

# Read differential gene expression (DGE) analysis results
# Source: mmc2.xlsx, sheet "Control_vs_COVID" (Supplementary Data 2)
# Contains genes with significant expression differences between control and COVID-19 patients
dge_results <- read_xlsx("DATA/mmc2.xlsx", sheet = "Control_vs_COVID")

# Remove any empty rows where the Gene column is NA
# This cleans the data for downstream processing
dge_results <- subset(dge_results, is.na(dge_results$Gene) == FALSE)

####################
### Human Genes Database ###
####################

# Load a sample gene list from the DOSE package to initialize the database
# This is used to access the gene annotation system
data(geneList, package = "DOSE")

# Create a Gene Ontology enrichment object to map Entrez Gene IDs to gene symbols
# This provides a mapping between gene identifiers (Entrez IDs) and gene symbols
# readable = TRUE converts IDs to human-readable gene symbols
gene_db <- enrichGO(names(geneList), OrgDb = org.Hs.eg.db, readable = TRUE)

# Convert the gene mapping to a data frame
# Creates two columns: Gene (symbol) and geneID (Entrez ID)
gene_db <- data.frame(Gene = gene_db@gene2Symbol, geneID = names(gene_db@gene2Symbol))

# Filter to only include genes that are present in our DGE results
# This ensures we only work with genes from our analysis
gene_db <- subset(gene_db, Gene %in% dge_results$Gene)


####################
### KEGG Pathways ###
####################

# Download KEGG pathway database for human (hsa = Homo sapiens)
# This retrieves all human pathways and their associated genes
pathway_db <- download_KEGG("hsa", keggType = "KEGG")

# Extract the mapping between pathway IDs and Entrez gene IDs
# KEGGPATHID2EXTID contains: pathway ID -> gene ID relationships
pathgenes_db <- pathway_db$KEGGPATHID2EXTID
colnames(pathgenes_db) <- c("hsaID", "geneID")

# Extract pathway names corresponding to pathway IDs
# KEGGPATHID2NAME contains: pathway ID -> pathway name relationships
pathway_db <- pathway_db$KEGGPATHID2NAME
colnames(pathway_db) <- c("hsaID", "Pathway")

# Merge pathway-gene mapping with our gene database
# This adds gene symbols to pathway-gene relationships
# Result: hsaID, geneID, Gene (symbol)
pathways <- merge(pathgenes_db, gene_db, by = "geneID")

# Add pathway names to complete the mapping
# Final result: hsaID, geneID, Gene (symbol), Pathway (name)
# This gives us a complete mapping: gene -> pathway for all genes in our dataset
pathways <- merge(pathways, pathway_db, by = "hsaID")


####################
### Save Results ###
####################

# Save processed DGE results to CSV for downstream analysis
# row.names = FALSE: Don't write row numbers as first column
# na = "": Write empty string for missing values
# quote = FALSE: Don't quote strings (reduces file size)
write.csv(dge_results, "R/DGE_results.csv", row.names = FALSE, na = "", quote = FALSE)

# Save gene-pathway mappings to CSV
# Contains: hsaID (KEGG pathway ID), geneID (Entrez ID), Gene (symbol), Pathway (name)
write.csv(pathways, "R/DGE_pathways.csv", row.names = FALSE, na = "", quote = FALSE)
