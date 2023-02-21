setRepositories(graphics = getOption("menu.graphics"),
                ind = NULL, addURLs = character())

#you need to install rtools

#  see https://cran.r-project.org/bin/windows/Rtools/index.html
#  After installation is complete, you need to perform one more step
#  to be able to compile R packages: you need to put the location 
#  of the Rtools make utilities (bash, make, etc) on the PATH. 
#  The easiest way to do so is to run the following line from Rstudio

writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")

# Now restart R, and verify that make can be found, 
#  which should show the path to your Rtools installation.

#test if worked with the following line

Sys.which("make")

#---------------------------------------------
#  for mac users:
#---------------------------------------------
#  You need to install xcode

#  https://clanfear.github.io/CSSS508/docs/compiling.html

# Open a terminal window
#  Click Spotlight search in the top right of your screen,
#  then search for "Terminal"
#  Copy and paste the following into the terminal,
#  then press enter.

# xcode-select --install

#  You will probably need to provide your password to 
#  enable installing the software.
#
#  Follow any onscreen instructions and wait for it to finish.
#
#------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.11")

rm(list = ls()); #erases all data in environment.

install.packages("tidyr")
install.packages("tidyverse")
install.packages("stringr")
install.packages("purrr")
install.packages("ggplot2")
install.packages("gplots")
install.packages("RColorBrewer")
install.packages("Matrix")

BiocManager::install("edgeR")
BiocManager::install("limma")
BiocManager::install("DESeq2")
BiocManager::install("apeglm")
BiocManager::install("biomaRt")
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("EnsDb.Hsapiens.v86")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("GenomicRanges")
BiocManager::install("ComplexHeatmap")
BiocManager::install("circlize")
BiocManager::install("dendextend")
BiocManager::install("stringr")
BiocManager::install("topGO")
BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
BiocManager::install("Rsamtools")
BiocManager::install("DiffBind")
BiocManager::install("SingleCellExperiment")
BiocManager::install("scater")
BiocManager::install("DropletUtils")
BiocManager::install("scran")
BiocManager::install("tidyverse")
BiocManager::install("Seurat")
BiocManager::install("celldex")
BiocManager::install("SingleR")

