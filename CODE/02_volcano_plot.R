# load required packages
library(ggplot2)   # plotting
library(ggrepel)   # add text

# set working directory
setwd("~/Project2")

# load data from previous script
dge_results <- read.csv("R/DGE_results.csv")


####################
### basic plot ###
####################

### BASE
ggplot(dge_results,                 # plot data
       aes(x = log2FoldChange,      # x axis
           y = -log10(pvalue))) +   # y axis
  
  ### PLOT
  geom_point()                        # add points


####################
### add colors ###
####################

# create empty column for expression variable
dge_results$Expression <- NA

# in the expression column, add whether gene is ...
# ... up-regulated (log2FoldChange > 0)
dge_results$Expression[dge_results$log2FoldChange > 0] <- "up"

# ... or down-regulated (log2FoldChange < 0)
dge_results$Expression[dge_results$log2FoldChange < 0] <- "down"

### BASE
ggplot(dge_results, 
       aes(x = log2FoldChange,     # x axis
           y = -log10(pvalue))) +  # y axis
  
  ### PLOT
  geom_point(aes(color = Expression))  # points with colors


####################
### add labels ###
####################

# genes cited in article
up_genes <- c("CXCL5", "CXCL12", "CCL2", "CCL4", "CXCL10", 
              "IFIH1", "IFI44", "IFIT1", "IL6", "IL10")
down_genes <- c("RPL41", "RPL17", "SLC25A6", "CALM1", "TUBA1A")

# create its own data frame
volcano_labels <- subset(dge_results, Gene %in% c(up_genes, down_genes))

### BASE
ggplot(dge_results, 
       aes(x = log2FoldChange,        # x axis
           y = -log10(pvalue))) +     # y axis
  
  ### PLOT
  geom_point(aes(color = Expression)) +  # points with colors
  geom_text_repel(data = volcano_labels, # source data
                  color = "black",       # black text
                  aes(x = log2FoldChange,
                      y = -log10(pvalue),
                      label = Gene))     # label mapping


####################
### named points ###
####################

### BASE
ggplot(dge_results,
       aes(x = log2FoldChange,
           y = -log10(pvalue))) +
  
  ### PLOT
  geom_point(aes(color = Expression)) +
  geom_text_repel(data = volcano_labels,
                  color = "black",           # black text labels
                  aes(x = log2FoldChange,
                      y = -log10(pvalue),
                      label = Gene)) +        # label mapping
  geom_point(data = volcano_labels,           # source data
             shape = 1,                       # hollow black circles
             color = "black",
             aes(x = log2FoldChange,          # same x axis
                 y = -log10(pvalue)))         # same y axis


####################
### dashed lines ###
####################

# add vertical lines at 
# log2FoldChange = 1 (doubled) and 
# log2FoldChange = -1 (halved)
# add horizontal line at y = 0

### BASE
ggplot(dge_results,
       aes(x = log2FoldChange,
           y = -log10(pvalue))) +
  
  ### PLOT
  geom_point(aes(color = Expression)) +
  geom_text_repel(data = volcano_labels,
                  color = "black",
                  aes(x = log2FoldChange,
                      y = -log10(pvalue),
                      label = Gene)) +
  geom_point(data = volcano_labels,
             shape = 1, color = "black",
             aes(x = log2FoldChange,
                 y = -log10(pvalue))) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +  # vertical dashed lines
  geom_hline(yintercept = 0, linetype = "dashed")           # horizontal dashed line


####################
### customize ###
####################

### BASE
volcano.plot <-
  ggplot(dge_results,
         aes(x = log2FoldChange,
             y = -log10(pvalue))) +
  
  ### PLOT
  geom_point(aes(color = Expression)) +
  geom_text_repel(data = volcano_labels,
                  color = "black",
                  aes(x = log2FoldChange,
                      y = -log10(pvalue),
                      label = Gene)) +
  geom_point(data = volcano_labels,
             shape = 1, color = "black",
             aes(x = log2FoldChange,
                 y = -log10(pvalue))) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  ### MODS
  scale_color_manual(values = c("blue", "red")) +   # change to blue & red
  theme_classic() +                                # remove grid lines
  theme(panel.border = element_rect(color = "grey", fill = NA),  # add grey border
        legend.position = "none") +                # remove legend
  labs(title = "Volcano plot representing upregulated and downregulated genes",
       subtitle = "Victoria Dixon")                # add title and your name

### view
volcano.plot

### save
ggsave("PLOTS/volcano.pdf", volcano.plot, device = "pdf",
       width = 6, height = 8, units = "in")

