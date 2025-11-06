install.packages("tidyverse")
install.packages("ggpubr")
install.packages("ggrepel")
install.packages("Rtools")
installed.packages("Rtools")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")

