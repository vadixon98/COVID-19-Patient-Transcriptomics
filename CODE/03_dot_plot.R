# =============================================================================
# Pathway Enrichment Dot Plot Visualization Script
# =============================================================================
# This script creates a dot plot to visualize KEGG pathway enrichment results.
# The plot shows up-regulated and down-regulated pathways, with dot size and
# color representing the number of genes and statistical significance, respectively.

# Load Required Packages
library(ggplot2)  # For creating plots and graphics
library(ggpubr)   # For arranging multiple plots together (ggarrange)

# Set Working Directory
# Note: Update this path to match your project directory
setwd("~/Project2")


####################
### Data Preparation ###
####################

# Load Data from Previous Scripts
# DGE results: differential gene expression analysis results
dge_results <- read.csv("R/DGE_results.csv")
# Pathways: gene-to-pathway mappings from KEGG database
pathways <- read.csv("R/DGE_pathways.csv")

# Define Pathways of Interest to Display
# Up-regulated pathways: typically immune response, inflammation, and signaling pathways
# These are pathways enriched in COVID-19 patients compared to controls
up_paths <- c("Cytokine-cytokine receptor interaction", "JAK-STAT signaling pathway",
              "Complement and coagulation cascades", "Hematopoietic cell lineage",
              "Chemokine signaling pathway", "Inflammatory bowel disease",
              "Toll-like receptor signaling pathway", "IL-17 signaling pathway",
              "TGF-beta signaling pathway", "Th1 and Th2 cell differentiation")

# Down-regulated pathways: typically metabolic and cellular function pathways
# These are pathways that show reduced activity in COVID-19 patients
down_paths <- c("Ribosome", "Oxidative phosphorylation", "Viral myocarditis",
                "Protein processing in endoplasmic reticulum", "Oxytocin signaling pathway",
                "Type I diabetes mellitus", "Phagosome", "Amyotrophic lateral sclerosis",
                "Ferroptosis", "Allograft rejection")

# Filter Pathways to Only Include Those of Interest
# This reduces the dataset to pathways we want to visualize
pathways <- subset(pathways, Pathway %in% c(up_paths, down_paths))

# Merge Pathway Data with DGE Results
# This adds gene expression statistics (fold change, p-values) to pathway information
# Links each gene in a pathway with its differential expression results
pathways <- merge(pathways, dge_results, by = "Gene")


####################
### Calculate Pathway Statistics ###
####################

# Create Empty Columns for Pathway-Level Summaries
# NumGenes: total number of genes in each pathway that are in our dataset
pathways$NumGenes <- NA

# pMean: mean adjusted p-value (padj) for all genes in each pathway
# Lower values indicate more significant enrichment
pathways$pMean <- NA

# Calculate Statistics for Each Pathway
# Loop through each pathway and compute summary statistics
for (path in c(up_paths, down_paths)) {                  # Iterate through pathway names
  df <- subset(pathways, Pathway == path)                # Extract data for this pathway
  n <- length(df$Gene)                                   # Count unique genes in pathway
  p <- mean(df$padj)                                     # Calculate mean adjusted p-value
  pathways$NumGenes[pathways$Pathway == path] <- n       # Store gene count
  pathways$pMean[pathways$Pathway == path] <- p          # Store mean p-value
}

# Create Summary Data Frame
# Extract unique pathway-level summaries (one row per pathway)
pathways <- unique(pathways[c("Pathway", "NumGenes", "pMean")])

# Split Pathways into Up- and Down-Regulated Groups
# This allows us to plot them separately with different styling
up_pathways <- subset(pathways, Pathway %in% up_paths)
down_pathways <- subset(pathways, Pathway %in% down_paths)


####################
### Basic Up-Regulated Pathways Plot ###
####################

# Create Basic Dot Plot for Up-Regulated Pathways
# X-axis: Number of genes in the pathway (indicates pathway size/coverage)
# Y-axis: Pathway names (one row per pathway)
# Color: Mean adjusted p-value (darker/lighter indicates statistical significance)
# Size: Number of genes (larger dots = more genes in pathway)
ggplot(up_pathways,
       aes(x = NumGenes,             # X-axis: number of genes per pathway
           y = Pathway,              # Y-axis: pathway names (categorical)
           color = pMean,            # Color scale: mean adjusted p-value (significance)
           size = NumGenes)) +       # Point size: number of genes (redundant but useful)
  geom_point()                       # Create dot plot


####################
### Order Pathways by Gene Count ###
####################

# Create Ordered List of Pathways
# Order pathways by number of genes (ascending: smallest to largest)
# This will be used to arrange pathways from bottom to top in the plot
up_order <- up_pathways$Pathway[order(up_pathways$NumGenes)]

### BASE PLOT LAYER
ggplot(up_pathways, aes(x = NumGenes, y = Pathway)) +
  
  ### PLOT ELEMENTS
  geom_point(aes(color = pMean, size = NumGenes)) +     # Colored and sized dots
  scale_y_discrete(limits = up_order)                    # Apply custom ordering to y-axis
                                                        # Pathways ordered by gene count



####################
### Final Up-Regulated Pathways Plot ###
####################

# Order Pathways by Number of Genes (Ascending Order)
# Smallest pathways at bottom, largest at top
up_order <- up_pathways$Pathway[order(up_pathways$NumGenes)]

# Create Finalized Up-Regulated Pathways Plot
up.plot <-
  ggplot(up_pathways, aes(x = NumGenes, y = Pathway)) +
  
  # Plot Elements
  geom_point(aes(color = pMean, size = NumGenes)) +           # Colored and sized dots
  
  # Customize Y-Axis Ordering
  scale_y_discrete(limits = up_order) +                        # Order pathways by gene count
  
  # Color Scale: Blue (low p-value, significant) to Red (high p-value, less significant)
  # Use full range from all pathways for consistent color scale across plots
  scale_color_gradient(low = "blue", high = "red",
                       limits = c(min(pathways$pMean, na.rm = TRUE),
                                  max(pathways$pMean, na.rm = TRUE))) +
  
  # Size Scale: Consistent size range across both plots
  # Ensures dot sizes are comparable between up and down plots
  scale_size_continuous(limits = c(min(pathways$NumGenes, na.rm = TRUE),
                                   max(pathways$NumGenes, na.rm = TRUE))) +
  
  # Theme and Labels
  theme_light() +                                              # Light theme with grid
  theme(axis.title = element_blank()) +                       # Remove axis titles
  labs(color = "Mean\nadjusted\np-value",                     # Legend for color
       size  = "Number\nof genes")                            # Legend for size

# Display the Plot
print(up.plot)



####################
### Down-Regulated Pathways Plot ###
####################

# Order Down-Regulated Pathways by Number of Genes (Ascending Order)
down_order <- down_pathways$Pathway[order(down_pathways$NumGenes)]

### BASE PLOT LAYER
down.plot <-
  ggplot(down_pathways, aes(x = NumGenes, y = Pathway)) +
  
  ### PLOT ELEMENTS
  geom_point(aes(color = pMean, size = NumGenes)) +           # Colored and sized dots
  
  # Customize Axes for Mirror Effect
  scale_y_discrete(limits = down_order, position = "right") +  # Y-axis labels on right side
  scale_x_reverse() +                                          # Reverse x-axis (mirror effect)
                                                               # Creates back-to-back appearance
  
  # Color Scale: Same as up-regulated plot for consistency
  scale_color_gradient(low = "blue", high = "red",
                       limits = c(min(pathways$pMean, na.rm = TRUE),
                                  max(pathways$pMean, na.rm = TRUE))) +
  
  # Size Scale: Same range as up-regulated plot
  scale_size_continuous(limits = c(min(pathways$NumGenes, na.rm = TRUE),
                                   max(pathways$NumGenes, na.rm = TRUE))) +
  
  # Theme and Labels
  theme_light() +                                              # Light theme with grid
  theme(axis.title = element_blank()) +                       # Remove axis titles
  labs(color = "Mean\nadjusted\np-value",                     # Legend for color
       size  = "Number\nof genes")                            # Legend for size

# Display the Plot
print(down.plot)



####################
### Combine Up and Down Plots ###
####################

# Ensure ggpubr Package is Loaded (for ggarrange function)
library(ggpubr)  # Package for arranging multiple plots

# Combine Both Plots into a Single Figure
# Creates a back-to-back dot plot showing both up- and down-regulated pathways
dot.plot <-
  ggarrange(up.plot, down.plot,
            ncol = 1,              # Stack plots vertically (up on top, down on bottom)
            common.legend = TRUE,  # Use a single shared legend (color and size scales)
            legend = "right",      # Place legend on the right side of the figure
            bgcolor = "white")     # White background color

# Display the Combined Figure
print(dot.plot)

# Export Combined Plot to PDF
# Save the final figure for publication or presentation
ggsave("PLOTS/dot.pdf", dot.plot,
       width = 6, height = 10, units = "in")
