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

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.21")

BiocManager::install("clusterProfiler", force = TRUE)
.libPaths()  # check paths
install.packages(c("boot", "Matrix"), lib = "C:/Users/<YourUserName>/R/win-library/4.5")

