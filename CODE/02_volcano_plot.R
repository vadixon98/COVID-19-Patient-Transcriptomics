# =============================================================================
# Volcano Plot Visualization Script
# =============================================================================
# This script creates a volcano plot to visualize differential gene expression
# results. A volcano plot displays statistical significance (-log10 p-value)
# versus fold change (log2FoldChange) for each gene.

# Load Required Packages
library(ggplot2)   # For creating plots and graphics
library(ggrepel)   # For adding non-overlapping text labels to plots

# Set Working Directory
# Note: Update this path to match your project directory
setwd("~/Project2")

# Load Data from Previous Script
# Reads the processed differential gene expression results
dge_results <- read.csv("R/DGE_results.csv")


####################
### Basic Volcano Plot ###
####################

# Create the basic volcano plot structure
# X-axis: log2FoldChange (positive = upregulated, negative = downregulated)
# Y-axis: -log10(pvalue) (higher values = more statistically significant)
# 
# A volcano plot helps identify genes that are both:
# - Statistically significant (high -log10 p-value, top of plot)
# - Biologically relevant (large absolute fold change, far from center)

### BASE PLOT LAYER
ggplot(dge_results,                 # Data frame containing DGE results
       aes(x = log2FoldChange,      # X-axis: log2 fold change (COVID vs Control)
           y = -log10(pvalue))) +   # Y-axis: negative log10 of p-value (significance)
  
  ### PLOT ELEMENTS
  geom_point()                      # Add points for each gene


####################
### Add Colors to Distinguish Up/Down Regulation ###
####################

# Create a new column to categorize genes as up- or down-regulated
# This will be used to color-code points in the volcano plot
dge_results$Expression <- NA

# Categorize genes based on fold change direction:
# ... Up-regulated: log2FoldChange > 0 (genes increased in COVID-19 patients)
dge_results$Expression[dge_results$log2FoldChange > 0] <- "up"

# ... Down-regulated: log2FoldChange < 0 (genes decreased in COVID-19 patients)
dge_results$Expression[dge_results$log2FoldChange < 0] <- "down"

### BASE PLOT LAYER
ggplot(dge_results, 
       aes(x = log2FoldChange,     # X-axis: fold change
           y = -log10(pvalue))) +  # Y-axis: statistical significance
  
  ### PLOT ELEMENTS
  geom_point(aes(color = Expression))  # Color points based on regulation direction


####################
### Add Gene Labels ###
####################

# Define genes of interest that were cited in the publication
# Up-regulated genes: chemokines and cytokines involved in immune response
up_genes <- c("CXCL5", "CXCL12", "CCL2", "CCL4", "CXCL10", 
              "IFIH1", "IFI44", "IFIT1", "IL6", "IL10")

# Down-regulated genes: typically involved in basic cellular functions
down_genes <- c("RPL41", "RPL17", "SLC25A6", "CALM1", "TUBA1A")

# Create a subset of the data containing only the genes we want to label
# This will be used to add text labels to specific points on the plot
volcano_labels <- subset(dge_results, Gene %in% c(up_genes, down_genes))

### BASE PLOT LAYER
ggplot(dge_results, 
       aes(x = log2FoldChange,        # X-axis: fold change
           y = -log10(pvalue))) +     # Y-axis: statistical significance
  
  ### PLOT ELEMENTS
  geom_point(aes(color = Expression)) +  # Colored points for all genes
  geom_text_repel(data = volcano_labels, # Use subset data for labels
                  color = "black",       # Black text color for labels
                  aes(x = log2FoldChange,
                      y = -log10(pvalue),
                      label = Gene))     # Map Gene column to label text
                                        # ggrepel automatically avoids overlaps


####################
### Add Hollow Circles to Highlight Labeled Genes ###
####################

# Add a second layer of points (hollow circles) around labeled genes
# This helps visually highlight the genes of interest on the plot

### BASE PLOT LAYER
ggplot(dge_results,
       aes(x = log2FoldChange,
           y = -log10(pvalue))) +
  
  ### PLOT ELEMENTS
  geom_point(aes(color = Expression)) +       # Colored points for all genes
  geom_text_repel(data = volcano_labels,
                  color = "black",            # Black text for gene labels
                  aes(x = log2FoldChange,
                      y = -log10(pvalue),
                      label = Gene)) +        # Add gene name labels
  geom_point(data = volcano_labels,           # Overlay points for labeled genes
             shape = 1,                       # Hollow circles (shape 1)
             color = "black",                 # Black border
             aes(x = log2FoldChange,          # Same x coordinate
                 y = -log10(pvalue)))         # Same y coordinate
                                             # Creates a ring effect around labeled points


####################
### Add Reference Lines ###
####################

# Add dashed reference lines to help interpret the plot:
# Vertical lines mark fold change thresholds:
#   - log2FoldChange = 1  (2-fold increase, doubled expression)
#   - log2FoldChange = -1 (2-fold decrease, halved expression)
# Horizontal line marks the baseline (y = 0)
# These lines help identify genes with meaningful fold changes

### BASE PLOT LAYER
ggplot(dge_results,
       aes(x = log2FoldChange,
           y = -log10(pvalue))) +
  
  ### PLOT ELEMENTS
  geom_point(aes(color = Expression)) +                    # Colored points
  geom_text_repel(data = volcano_labels,
                  color = "black",
                  aes(x = log2FoldChange,
                      y = -log10(pvalue),
                      label = Gene)) +                     # Gene labels
  geom_point(data = volcano_labels,
             shape = 1, color = "black",
             aes(x = log2FoldChange,
                 y = -log10(pvalue))) +                    # Highlight circles
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +  # Vertical lines at ±2-fold change
  geom_hline(yintercept = 0, linetype = "dashed")           # Horizontal line at baseline


####################
### Final Customization and Export ###
####################

# Create the final polished volcano plot with custom styling
# Save it as an object so it can be viewed and exported

### BASE PLOT LAYER
volcano.plot <-
  ggplot(dge_results,
         aes(x = log2FoldChange,
             y = -log10(pvalue))) +
  
  ### PLOT ELEMENTS
  geom_point(aes(color = Expression)) +                  # Colored points
  geom_text_repel(data = volcano_labels,
                  color = "black",
                  aes(x = log2FoldChange,
                      y = -log10(pvalue),
                      label = Gene)) +                   # Gene labels
  geom_point(data = volcano_labels,
             shape = 1, color = "black",
             aes(x = log2FoldChange,
                 y = -log10(pvalue))) +                  # Highlight circles
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +  # Fold change thresholds
  geom_hline(yintercept = 0, linetype = "dashed") +     # Baseline
  
  ### CUSTOMIZATION
  scale_color_manual(values = c("blue", "red")) +       # Blue for down, red for up
  theme_classic() +                                     # Clean theme without grid
  theme(panel.border = element_rect(color = "grey", fill = NA),  # Grey border around plot
        legend.position = "none") +                     # Remove legend (colors are intuitive)
  labs(title = "Volcano plot representing upregulated and downregulated genes",
       subtitle = "Victoria Dixon")                     # Add title and author name

### VIEW PLOT
# Display the plot in the R graphics device
volcano.plot

### SAVE PLOT
# Export the plot as a PDF file for publication or presentation
ggsave("PLOTS/volcano.pdf", volcano.plot, device = "pdf",
       width = 6, height = 8, units = "in")

