# About
# This script reads in raw data > preprocesses data > to form a data matrix
# Data matrix is used for further stats analysis

# ============== Clears Memory ======
# clears all objects includes hidden objects
rm(list = ls(all.names = TRUE)) 

# frees up memory and reports the memory usage.
gc() 

# ============== Loads Packages =======
library(readxl)
library(dplyr)
library(data.table)
library(stringr)
library(tidyverse)
library(ggrepel)

# ============== 1. Read & Selects from Excel File ======
feno_raw <- fread('FENO_MOUSE_PG.csv') %>%
  na_if(., 0) %>%
  # renames Gene Symbol column
  rename("Gene Symbol" = "PG.ProteinAccessions") %>%
  
  # removes NAs in Gene Symbol
  na.omit("Gene Symbol") %>%
  
  # adds number to the end of duplicate gene symbols (ie Sptbn1-2)
  group_by(`Gene Symbol`) %>%
  mutate(`GS_count` = 1:n()) %>%
  mutate(`Gene Symbol` = ifelse(`GS_count` == 1, 
                                `Gene Symbol`, 
                                paste0(`Gene Symbol`, "-", `GS_count`))) %>%
  select(-c(`GS_count`))
  mutate_all(function(x) as.numeric(as.character(x))) %>%
  na.omit()

  
library(IMIFA)

transform_data <- function(x) {
  x <- log10(x)
  rowmed <- apply(x,1,median)
  x <- sweep(x,1,rowmed,"-")
  x <- data.frame(pareto_scale(x, centering = TRUE))
} 
# applies function to transform data
feno_normalized <- transform_data(feno_raw[,-1])  

# ====== Scales data =====
library(IMIFA)
transform_data <- function(x) {
  x <- log10(x)
  rowmed <- apply(x,1,median)
  x <- sweep(x,1,rowmed,"-")
  x <- data.frame(pareto_scale(x, centering = TRUE))
} 
# applies function to transform data
feno_normalized <- transform_data(feno_raw[,-1])

rownames(feno_normalized) <- feno_raw$`Gene Symbol`

# ===== runs PCA =====
res.pca <- PCA(feno_normalized,  graph = FALSE)

fviz_pca_biplot(res.pca,
                repel = FALSE,
                geom = "point",
                show.clust.cent = TRUE,
                label = "all",
                ggtheme = theme_minimal(),
                alpha = 0)
