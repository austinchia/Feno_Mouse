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

# reads raw data
feno_raw <- read_excel('Feno_Mouse_Dataset.xlsx', sheet = 'Feno_Mouse_Dataset', na = c("", "NA")) %>%
  # replaces "Filtered" with NAs
  na_if(., "Filtered") %>%
  
  # renames Gene Symbol column
  rename("Gene Symbol" = "PG.ProteinNames") %>%
  
  # removes NAs in Gene Symbol
  na.omit(`Gene Symbol`) %>%
  
  # adds number to the end of duplicate gene symbols (ie Sptbn1-2)
  group_by(`Gene Symbol`) %>%
  mutate(`GS_count` = 1:n()) %>%
  mutate(`Gene Symbol` = ifelse(`GS_count` == 1, 
                                `Gene Symbol`, 
                                paste0(`Gene Symbol`, "-", `GS_count`))) %>%
  # column_to_rownames(., var = "Gene Symbol") %>%
  
  # remove useless columns
  select(., -c(`PG.ProteinAccessions`, `PG.ProteinDescriptions`)) %>%
  
  mutate_all(function(x) as.numeric(as.character(x)))
library(IMIFA)
transform_data <- function(x) {
  x <- log10(x)
  rowmed <- apply(x,1,median)
  x <- sweep(x,1,rowmed,"-")
  x <- data.frame(pareto_scale(x, centering = TRUE))
} 
# applies function to transform data
feno_normalized <- transform_data(feno_raw[,-1])

#  === runs method 1 pca =====

pca_res <- prcomp(t(feno_normalized), scale=TRUE)
summary(pca_res)



# var_explained <- pca_res$sdev^2/sum(pca_res$sdev^2)
# 
# pca_res$x %>% 
#   as.data.frame %>%
#   ggplot(aes(x=PC1,y=PC2)) + geom_point(size=1) +
#   theme_bw(base_size=32) + 
#   labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
#        y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
#   theme(legend.position="top")

# ==== runs method 2 pca =====

pca.var <- pca_res$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

pca.data <- data.frame(Sample = rownames(pca_res$x),
                       X = pca_res$x[,1],
                       Y = pca_res$x[,2])

ggplot(data=pca.data, aes(x=X, y=Y, label=Sample)) +
  geom_text_repel(max.overlaps = Inf, point.size = NA) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("My PCA Graph")

# ==== factoextra method 3 ====
library(factoextra)
library("FactoMineR")

res.pca <- PCA(t(feno_normalized),  graph = FALSE)


fviz_pca_biplot(pca_res,
             repel = FALSE,
             geom = "point",
             show.clust.cent = TRUE,
             ellipse.type = "norm",
             ggtheme = theme_minimal(),
             addEllipses = FALSE,
             col.ind=df$Kingdom,
             alpha = 0)

# =========== Metaboanalyst ================

library(MetaboAnalystR)

  # initializing object
  mSet <- InitDataObjects("pktable", "stat", FALSE)
  # loading in data
  mSet <- Read.TextData(mSet, "FENO_MOUSE_PG_OG.csv", "colu", "disc");
  # data check
  mSet <- SanityCheckData(mSet)
  mSet <- ReplaceMin(mSet);
  
  # filtering features
  mSet <- FilterVariable(mSet, "none", "F", 25)
  
  # median normalization, log10 transformation, pareto scaling
  mSet <- PreparePrenormData(mSet)
  mSet <- Normalization(mSet, "MedianNorm", "LogNorm", "ParetoNorm", ratio=FALSE, ratioNum=20)
  mSet <- PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
  mSet <- PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)
  
  mSet<-PCA.Anal(mSet)
  mSet<-PlotPCAPairSummary(mSet, "pca_pair_0_", "png", 72, width=NA, 5)
  mSet<-PlotPCA2DScore(mSet, "pca_score2d_5_", "png", 72, width=NA, 1,2,0.95,0,0)



