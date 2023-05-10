## Digital Epidemiology - HW 01
## Authors: Cruoglio Antonella, Mascolo Davide, Napoli Mario

## This file is used only to get data from GDC Data Portal

# Import Utils ------------------------------------------------------------

library(dplyr)
library(psych)
library(ggplot2)
library(ggnet)
library(igraph)
library(GGally)
library(sna)
library(network)
library(gplots)
library(biomaRt)
library(stringr)
library(TCGAbiolinks)
library(SummarizedExperiment)

# Download Data -----------------------------------------------------------
## Create directory
proj <- "TCGA-THCA"
dir.create(file.path(proj))

## Query for cancer data
rna_query_C <- GDCquery(project = proj,
                        data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification",
                        workflow.type = "STAR - Counts",
                        sample.type = "Primary Tumor")
## Download data
GDCdownload(query = rna_query_C, directory = "GDCdata",
            method = "api")
rna_data_C      <- GDCprepare(rna_query_C)
norm_data_C     <- rna_data_C@assays@data@listData[["fpkm_unstrand"]]

## Get Expression data
rna_expr_data_C    <- assay(rna_data_C)
## Get genomics metadata
rna_genes_info_C   <- rowRanges(rna_data_C)
## Get Metadata
rna_samples_info_C <- colData(rna_data_C)

## Set row and col names
rownames(norm_data_C) <- rownames(rna_expr_data_C)
colnames(norm_data_C) <- colnames(rna_expr_data_C)

## Write csv
## Save normalize data and other info
write.csv(norm.data.C, "norm_rna_expr_data_C.csv", sep = ",",
          row.names = TRUE, col.names = TRUE, quote = FALSE)
## Non-normalized data
write.csv(rna_expr_data_C, "rna_expr_data_C.csv",
          sep = ",", row.names = TRUE, col.names = TRUE,
          quote = FALSE)
write.table(rna_samples_info_C@listData$patient,
            file = file.path(proj,paste(proj, "rna_patients_C.txt",
                                        sep = "")),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(rna_genes_info_C@ranges@NAMES,
            file = file.path(proj,paste(proj, "rna_genes_C.txt",
                                        sep = "")),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

## ****************************************************************** 

## Query for non-cancer data
rna_query_N <- GDCquery(project = proj,
                        data.category = "Transcriptome Profiling", 
                        data.type = "Gene Expression Quantification", 
                        workflow.type = "STAR - Counts", 
                        sample.type = "Solid Tissue Normal")
## Download data
GDCdownload(query = rna_query_N, directory = "GDCdata",
            method = "api")
rna_data_N      <- GDCprepare(rna_query_N)
rna_expr_data_N <- assay(rna_data_N)
norm_data_N     <- rna_data_N@assays@data@listData[["fpkm_unstrand"]]

## Set row and col names
rownames(norm_data_N) <- rownames(rna_expr_data_N)
colnames(norm_data_N) <- colnames(rna_expr_data_N)

## Write csv
## Save normalize data and other info
write.csv(norm_data_N, "norm_rna_expr_data_N.csv",sep = ",",
          row.names = TRUE, col.names = TRUE, quote = FALSE)
## Non-normalized data
write.csv(rna_expr_data_C, "rna_expr_data_C.csv",
          sep = ",", row.names = TRUE, col.names = TRUE,
          quote = FALSE)
write.table(rna_samples_info_C@listData$patient,
            file = file.path(proj,paste(proj, "rna_patients_C.txt",
                                        sep = "")),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(rna_genes_info_C@ranges@NAMES,
            file = file.path(proj,paste(proj, "rna_genes_C.txt",
                                        sep = "")),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

## clinical info
clinical_query <- GDCquery_clinic(project = proj, type = "clinical",
                                  save.csv = FALSE)
