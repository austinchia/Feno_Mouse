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

# ============== 1. Read & Selects from Excel File ======

# reads raw data
feno_raw <- read_excel('Feno_Mouse_Dataset.xlsx', sheet = 'Feno_Mouse_Dataset', na = c("", "NA")) %>%
  # replaces "Filtered" with NAs
  na_if(., "Filtered") %>%
  
  # removes NAs
  na.omit() %>%
  
  # renames Gene Symbol column
  rename("Gene Symbol" = "PG.ProteinNames") %>%
  
  # adds number to the end of duplicate gene symbols (ie Sptbn1-2)
  group_by(`Gene Symbol`) %>%
  mutate(`GS_count` = 1:n()) %>%
  mutate(`Gene Symbol` = ifelse(`GS_count` == 1, 
                                `Gene Symbol`, 
                                paste0(`Gene Symbol`, "-", `GS_count`))) %>%
  column_to_rownames(., var = "Gene Symbol") %>%
  
  # remove useless columns
  select(., -c(`PG.ProteinAccessions`, `PG.ProteinDescriptions`)) %>%
  
  mutate_all(function(x) as.numeric(as.character(x)))

# plots pca

pca_res <- prcomp(feno_raw, scale=TRUE)
summary(pca_res)

var_explained <- pca_res$sdev^2/sum(pca_res$sdev^2)

pca_res$x %>% 
  as.data.frame %>%
  ggplot(aes(x=PC1,y=PC2)) + geom_point(size=1) +
  theme_bw(base_size=32) + 
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  theme(legend.position="top")
