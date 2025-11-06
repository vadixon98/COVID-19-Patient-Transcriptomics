####################
### dot plot ###
####################

# load required packages
library(ggplot2)
library(ggpubr)

# set working directory
setwd("~/Project2")


####################
### data prep ###
####################

# load data from previous script
dge_results <- read.csv("R/DGE_results.csv")
pathways <- read.csv("R/DGE_pathways.csv")

# pathways shown in plot
up_paths <- c("Cytokine-cytokine receptor interaction", "JAK-STAT signaling pathway",
              "Complement and coagulation cascades", "Hematopoietic cell lineage",
              "Chemokine signaling pathway", "Inflammatory bowel disease",
              "Toll-like receptor signaling pathway", "IL-17 signaling pathway",
              "TGF-beta signaling pathway", "Th1 and Th2 cell differentiation")

down_paths <- c("Ribosome", "Oxidative phosphorylation", "Viral myocarditis",
                "Protein processing in endoplasmic reticulum", "Oxytocin signaling pathway",
                "Type I diabetes mellitus", "Phagosome", "Amyotrophic lateral sclerosis",
                "Ferroptosis", "Allograft rejection")

# subset by those pathways
pathways <- subset(pathways, Pathway %in% c(up_paths, down_paths))

# combine with DGE results
pathways <- merge(pathways, dge_results, by = "Gene")


####################
### average p values ###
####################

# create empty columns for number of genes per pathway ...
pathways$NumGenes <- NA

# ... and the mean adjusted p value per pathway
pathways$pMean <- NA

for (path in c(up_paths, down_paths)) {                  # loop through pathway names
  df <- subset(pathways, Pathway == path)                # subset pathway data by pathway name
  n <- length(df$Gene)                                   # get number of genes in that pathway
  p <- mean(df$padj)                                     # calculate mean adjusted p value
  pathways$NumGenes[pathways$Pathway == path] <- n       # add number of genes to column
  pathways$pMean[pathways$Pathway == path] <- p          # add mean p value to column
}

# subset to up/down regulated pathways
pathways <- unique(pathways[c("Pathway", "NumGenes", "pMean")])
up_pathways <- subset(pathways, Pathway %in% up_paths)
down_pathways <- subset(pathways, Pathway %in% down_paths)


####################
### up-regulated pathways plot ###
####################

# base plot for up-regulated pathways
ggplot(up_pathways,
       aes(x = NumGenes,             # x-axis: number of genes
           y = Pathway,              # y-axis: pathway name
           color = pMean,            # color mapped to mean adjusted p-value
           size = NumGenes)) +       # size mapped to number of genes
  geom_point()


####################
### order axis ###
####################

# y axis is ordered by number of genes
up_order <- up_pathways$Pathway[order(up_pathways$NumGenes)]

### BASE
ggplot(up_pathways, aes(x = NumGenes, y = Pathway)) +
  
  ### PLOT
  geom_point(aes(color = pMean, size = NumGenes)) +
  scale_y_discrete(limits = up_order)  # order y axis



####################
### up-regulated pathways ###
####################

# order y-axis by number of genes (largest at top)
up_order <- up_pathways$Pathway[order(up_pathways$NumGenes)]

up.plot <-
  ggplot(up_pathways, aes(x = NumGenes, y = Pathway)) +
  geom_point(aes(color = pMean, size = NumGenes)) +
  scale_y_discrete(limits = up_order) +                        # reorder by NumGenes
  scale_color_gradient(low = "blue", high = "red",
                       limits = c(min(pathways$pMean, na.rm = TRUE),
                                  max(pathways$pMean, na.rm = TRUE))) +
  scale_size_continuous(limits = c(min(pathways$NumGenes, na.rm = TRUE),
                                   max(pathways$NumGenes, na.rm = TRUE))) +
  theme_light() +
  theme(axis.title = element_blank()) +
  labs(color = "Mean\nadjusted\np-value",
       size  = "Number\nof genes")

# show the plot
print(up.plot)



####################
### down-regulated pathways ###
####################

# order y axis by number of genes
down_order <- down_pathways$Pathway[order(down_pathways$NumGenes)]

### BASE
down.plot <-
  ggplot(down_pathways, aes(x = NumGenes, y = Pathway)) +
  
  ### PLOT
  geom_point(aes(color = pMean, size = NumGenes)) +
  scale_y_discrete(limits = down_order, position = "right") +  # ordered y-axis on the right
  scale_x_reverse() +                                          # reverse x-axis
  scale_color_gradient(low = "blue", high = "red",
                       limits = c(min(pathways$pMean, na.rm = TRUE),
                                  max(pathways$pMean, na.rm = TRUE))) +
  scale_size_continuous(limits = c(min(pathways$NumGenes, na.rm = TRUE),
                                   max(pathways$NumGenes, na.rm = TRUE))) +
  theme_light() +
  theme(axis.title = element_blank()) +
  labs(color = "Mean\nadjusted\np-value",
       size  = "Number\nof genes")

# display the plot
print(down.plot)



####################
### combine ###
####################

library(ggpubr)  # ensure loaded for ggarrange()

dot.plot <-
  ggarrange(up.plot, down.plot,
            ncol = 1,              # stack vertically
            common.legend = TRUE,  # shared legend
            legend = "right",      # legend on right
            bgcolor = "white")     # white background

# display the combined figure
print(dot.plot)

# export to PDF
ggsave("PLOTS/dot.pdf", dot.plot,
       width = 6, height = 10, units = "in")
