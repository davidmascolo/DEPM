## Digital Epidemiology - HW 01
## Authors: Cruoglio Antonella, Mascolo Davide, Napoli Mario


# Import Utils ------------------------------------------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DT)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(psych)
library(ggnet)
library(igraph)
library(GGally)
library(sna)
library(network)
library(gplots)
library(biomaRt)
library(stringr)
library(DESeq2)
library(WGCNA)



# Load Data ---------------------------------------------------------------
## Import cancer patients
rna_data_C <- read.csv("rna_expr_data_C.csv", header = T)
## Import non-cancer patients
rna_data_N <- read.csv("rna_expr_data_N.csv", header = T)
## Check row names
x_rna_data_C <- rna_data_C$X
rna_data_C$X <- NULL

x_rna_data_N <- rna_data_N$X
rna_data_N$X <- NULL

## Set row names
rownames(rna_data_C) <- x_rna_data_C
rownames(rna_data_N) <- x_rna_data_N


## ************************************************************


# 1. Preprocessing -----------------------------------------------------------
## Match columns
ncol(rna_data_C) #504
ncol(rna_data_N) #59

## Row names of clinical query and column names of norm.expr.data
## do not match
## The patient is codified by the first 12 charachters 

## Normal patients
unique(substr(colnames(rna_data_N), 1,12))
length(unique(substr(colnames(rna_data_N), 1,12))) #59
## Cancer patients
unique(substr(colnames(rna_data_C), 1,12))
length(unique(substr(colnames(rna_data_C), 1,12))) #504
## For cancer patients we have more samples

## who has more samples?
patients_N <- substr(colnames(rna_data_N), 1,12)
table(patients_N)
table(patients_N)[order(table(patients_N))]
## We have the same number of samples for all the non-cancer
## patients

## Save their index in the column list
patients_N_idx <- match(patients_N,
                        substr(colnames(rna_data_C),
                               1,12))

## Now, we consider only the patients for which there is one sample
## Cancer patients
expr_C <- as.data.frame(rna_data_C)
expr_C <- expr_C[,patients_N_idx]
## Non-cancer patients
expr_N <- as.data.frame(rna_data_N)

## Renaming the patients in the expression matrices
## Cancer patients
colnames(expr_C) <- substr(colnames(expr_C), 1,12)
colnames(expr_C) 
length(colnames(expr_C)) ## 59

## Non-cancer patients
colnames(expr_N) <- substr(colnames(expr_N), 1,12)
colnames(expr_N)
length(colnames(expr_N)) ## 59

## Matching
intersect(ncol(expr_N), ncol(expr_C))
## Check
setdiff(colnames(expr_N), colnames(expr_C))

## By this point, we consider only the common patients
expr_C <- expr_C %>%
  dplyr::select(intersect(colnames(expr_C), colnames(expr_N))) 

## Check Na values
table(is.na(expr_C))
table(is.na(expr_N))
## There aren't Na

## How many genes have no zeros in the data frame?
## Cancer patients
sum(rowSums(expr_C == 0) == 0) ## 18904
## Extract their names
no_zeros_genes_C <- rownames(expr_C)[rowSums(expr_C == 0) == 0]

## Non-cancer patients
sum(rowSums(expr_N == 0) == 0) ## 20295
no_zeros_genes_N <- rownames(expr_N)[rowSums(expr_N == 0) == 0] #these are their names

length(intersect(no_zeros_genes_C, no_zeros_genes_N)) ## 18322

## Let's consider only these genes 
filtr_expr_C <- expr_C[intersect(no_zeros_genes_C,
                                 no_zeros_genes_N),]
filtr_expr_N <- expr_N[intersect(no_zeros_genes_N,
                                 no_zeros_genes_C),]
## Match
all(rownames(filtr_expr_C) == rownames(filtr_expr_N))

## Final step
final_expr_C <- filtr_expr_C
final_expr_N <- filtr_expr_N


## ************************************************************


# 2. Differentially Expressed Genes (DEGs) --------------------------------

## Add informations to columns name
## Cancer patients
for (i in 1:ncol(final_expr_C)){
  colnames(final_expr_C)[i] <- paste(colnames(final_expr_C)[i],
                                     "_cancer", sep = "")
}
## Normal patients
for (i in 1:ncol(final_expr_N)){
  colnames(final_expr_N)[i] <- paste(colnames(final_expr_N)[i],
                                     "_normal", sep = "")
}
## Check
colnames(final_expr_C)
colnames(final_expr_N)

## Merge data
dat_final <- cbind(final_expr_C, final_expr_N)

## Check
dim(dat_final)

## Create a new data with information about the status of the patients
patients_condition <- data.frame(Patients = colnames(dat_final),
                                 Condition = c(rep("cancer",
                                                   ncol(final_expr_C)),
                                               rep("normal",
                                                   ncol(final_expr_N))))
patients_condition$Condition    <- factor(patients_condition$Condition)
rownames(patients_condition)    <- patients_condition$Patients
patients_condition$Patients <- NULL
## Creating DESeq2 object
DESeq_dt <- DESeqDataSetFromMatrix(countData = dat_final,
                                   colData = patients_condition,
                                   design = ~ Condition)
dim(DESeq_dt)
## 18323 rows and 118 columns

## Pre-filtering
## We apply a pre-filtering operation, selecting only the genes
## that have at least 10 reads.
keep <- rowSums(counts(DESeq_dt)) >= 10
DESeq_dt_filtering <- DESeq_dt[keep, ]
DESeq_dt_filtering$Condition <- relevel(DESeq_dt_filtering$Condition,
                                        ref = "normal")

dim(DESeq_dt_filtering)
## As we can see, we reduce the number of observations (genes)

## Identify DEGs specifying the thresholds setting using p-value and
## Fold Change FC >= |1.2|
## We apply Benjamini-Hochberg correction for multiple testing (FDR)
DESeq_output <- DESeq(DESeq_dt_filtering)
res          <- results(DESeq_output, alpha = 0.05,
                        lfcThreshold = 1.2, pAdjustMethod = "BH")
## Check
summary(res)
## We obtain a set of 669 genes
## LFC > 1.20 ---> Up-regulated genes   (520, 2.80%)
## LFC <-1.20 ---> Down regulated genes (149, 0.81%)

## MA Plot 
ggpubr::ggmaplot(res, main = expression("Group 1" %->% "Group 2"),
                 fdr = 0.05, fc = 1.2, size = 0.7,
                 palette = c("#B31B21", "#1465AC", "grey"),
                 legend = "top", top = 20,
                 font.label = c("bold", 11),
                 font.legend = "bold",
                 font.main = "bold",
                 ggtheme = theme_minimal())

## Volcano Plot
volcano_dat <- as.data.frame(res)
volcano_dat <- volcano_dat %>% 
  mutate(Expression = case_when(log2FoldChange >= 1.2 & padj <= 0.05 ~ "Up-regulated",
                                log2FoldChange <= -1.2 & padj <= 0.05 ~ "Down-regulated",
                                TRUE ~ "Unchanged"))
## Plot
ggplot(data = volcano_dat, aes(x = log2FoldChange,
                               y = -log10(padj),
                               col = Expression)) +
  geom_point() + 
  theme_minimal() +
  scale_color_manual(values = c("green", "grey", "red")) +
  geom_vline(xintercept = c(-1.2, 1.2), col = "black") +
  geom_hline(yintercept = -log10(0.05), col = "black") +
  ggtitle("Volcano Plot")


## ************************************************************


# 3. Co-expression networks -----------------------------------------------
## Computation
## In this step, we use the normalize data
## Import cancer patients
norm_rna_data_C <- read.csv("norm_rna_expr_data_C.csv", header = T)
## Import non-cancer patients
norm_rna_data_N <- read.csv("norm_rna_expr_data_N.csv", header = T)

# 3.1 Pre-processing step -------------------------------------------------
## Check row names
x_rna_data_C <- norm_rna_data_C$X
norm_rna_data_C$X <- NULL
x_rna_data_N <- norm_rna_data_N$X
norm_rna_data_N$X <- NULL

## Set row names
rownames(norm_rna_data_C) <- x_rna_data_C
rownames(norm_rna_data_N) <- x_rna_data_N

## Select index of DEGs genes
filter_norm_data <- volcano_dat %>%
  filter(Expression !=  "Unchanged")
## Check
dim(filter_norm_data) ## 669

## Normal patients
unique(substr(colnames(norm_rna_data_N), 1,12))
length(unique(substr(colnames(norm_rna_data_N), 1,12))) #59
## Cancer patients
unique(substr(colnames(norm_rna_data_C), 1,12))
length(unique(substr(colnames(norm_rna_data_C), 1,12))) #504
## For cancer patients we have more samples

## who has more samples?
patients_N <- substr(colnames(norm_rna_data_N), 1,12)
table(patients_N)
table(patients_N)[order(table(patients_N))]
## We have the same number of samples for all the non-cancer
## patients

## Save their index in the column list
patients_N_idx <- match(patients_N,
                        substr(colnames(rna_data_C),
                               1,12))

## Now, we consider only the patients for which there is one sample
## Cancer patients
expr_C <- as.data.frame(norm_rna_data_C)
expr_C <- expr_C[,patients_N_idx]
## Non-cancer patients
expr_N <- as.data.frame(norm_rna_data_N)

## Renaming the patients in the expression matrices
## Cancer patients
colnames(expr_C) <- substr(colnames(expr_C), 1,12)
colnames(expr_C) 
length(colnames(expr_C)) ## 59

## Non-cancer patients
colnames(expr_N) <- substr(colnames(expr_N), 1,12)
colnames(expr_N)
length(colnames(expr_N)) ## 59

## Matching
intersect(ncol(expr_N), ncol(expr_C))
## Check
setdiff(colnames(expr_N), colnames(expr_C))

## By this point, we consider only the common patients
expr_C <- expr_C %>%
  dplyr::select(intersect(colnames(expr_C), colnames(expr_N))) 

## Check Na values
table(is.na(expr_C))
table(is.na(expr_N))
## There aren't Na

## How many genes have no zeros in the data frame?
## Cancer patients
sum(rowSums(expr_C == 0) == 0) ## 18904
## Extract their names
no_zeros_genes_C <- rownames(expr_C)[rowSums(expr_C == 0) == 0]

## Non-cancer patients
sum(rowSums(expr_N == 0) == 0) ## 20295
no_zeros_genes_N <- rownames(expr_N)[rowSums(expr_N == 0) == 0] #these are their names

length(intersect(no_zeros_genes_C, no_zeros_genes_N)) ## 18322

## Let's consider only these genes 
filtr_expr_C <- expr_C[intersect(no_zeros_genes_C,
                                 no_zeros_genes_N),]
filtr_expr_N <- expr_N[intersect(no_zeros_genes_N,
                                 no_zeros_genes_C),]
## Match
all(rownames(filtr_expr_C) == rownames(filtr_expr_N))

## Final step
norm_final_expr_C <- filtr_expr_C
norm_final_expr_N <- filtr_expr_N

## Subset DEGs reads using normalize data
## Cancer patients
norm_coexp_C <- norm_final_expr_C %>% 
  filter(row.names(norm_final_expr_C) %in%
           row.names(filter_norm_data))
## Non-cancer patients
norm_coexp_N <- norm_final_expr_N %>% 
  filter(row.names(norm_final_expr_N) %in%
           row.names(filter_norm_data))
## Check
dim(norm_coexp_C)
dim(norm_coexp_N)


## ************************************************************


# 3.2 Computation ---------------------------------------------------------

## Implement a function that computes the hubs of the network
compute_adjacency <- function(data, cor_type = NULL){
  ## Input:  Expression Data
  ## Output: Adjacency matrix
  
  ## Correlation matrix 
  cor_mat       <- cor(t(data), method = cor_type)
  diag(cor_mat) <- 0
  ## Correlation matrix (p-value)
  cor_padj <- corr.p(cor_mat, nrow(cor_mat),
                     adjust = "fdr", ci = FALSE)$p
  cor_padj[lower.tri(cor_padj)] <- t(cor_padj)[lower.tri(cor_padj)]
  ## Build adjacency matrix
  adj_mat_1 <- ifelse(cor_mat >= 0.7, 1,
                      ifelse(cor_mat <= -0.7,-1, 0))
  adj_mat_2 <- ifelse(abs(cor_padj) > 0.05, 0, 1) 
  adj_mat   <- adj_mat_1 * adj_mat_2
  
  ## Return
  return(adj_mat)
}

## Implement a function that returns the list of the hubs
mart_names <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
compute_hubs <- function(adj_matrix, mart = mart_names){
  ## Input: Adjacency Matrix, mart
  ## Output: List of hubs, Degree, Level of quantile
  
  ## Compute degree
  degree <- sort(rowSums(abs(adj_matrix)), decreasing = TRUE)
  ## Compute quantile
  q_95 <- quantile(degree[degree > 0], 0.95)
  ## Find the hubs (5% of the nodes with highest degree values)
  hubs <- degree[degree >= q_95]
  
  ## Let's get the hubs' names
  hubs_clean <- (str_split(names(hubs),
                           "[.]", simplify = TRUE))[,1]
  hubs_names <- getBM(filters = "ensembl_gene_id",
                      attributes = c("ensembl_gene_id", "hgnc_symbol"),
                      values = hubs_clean, mart = mart)
  
  ## Let's order them by degree
  hubs_ord <- (str_split(names(sort(hubs, decreasing = TRUE)),
                         "[.]",simplify = TRUE))[,1]
  hubs_names <- hubs_names[match(hubs_ord,
                                 hubs_names$ensembl_gene_id),]
  ## Check names
  final_hubs <- hubs_names[,2]
  
  ## Return
  return(list("code_hubs" = hubs,
              "final_hubs" = final_hubs,
              "degree" = degree,
              "q_95" = q_95))
}

## Implement a function that plot the network with the hubs
plot_graph <- function(net, hubs, title, subtitle = NULL){
  ## Input:  Network, Hubs, Title, Subtitle
  ## Output: Network Plot
  
  ## Get Names
  net %v% "type"  <- ifelse(network.vertex.names(net) %in%
                              names(hubs),
                            "hub", "non-hub")
  net %v% "color" <- ifelse(net %v% "type" == "hub",
                            "red", "blue")
  set.edge.attribute(net, "edgecolor",
                     ifelse(net %e% "weights" > 0, "green", "blue"))
  
  ## Get vertex
  coord <- gplot.layout.fruchtermanreingold(net, NULL)
  net %v% "x" = coord[, 1]
  net %v% "y" = coord[, 2]
  
  ## Return Plot
  ggnet2(net, color = "color", alpha = 0.7, size = 2,
         mode = c("x","y"),
         edge.color = "edgecolor", edge.alpha = 1,
         edge.size = 0.15) + guides(size = "none") +
    ggtitle(title,
            subtitle = subtitle)
}

## Cancer
## Adjacency Matrix Cancer (Pearson)
adj_mat_c <- compute_adjacency(norm_coexp_C, cor_type = "pearson")

# 3.3 - Analysis ----------------------------------------------------------

## Cancer hubs
hubs_c <- compute_hubs(adj_mat_c)

## Plot degree distribution in order to check if the network
## is scale free
df1           <- data.frame(cbind(hubs_c$degree))
colnames(df1) <- "Degree"
(hist_cancer <- ggplot(df1, aes(x = Degree)) +
  geom_histogram(fill = "blue", alpha = 0.7, bins = 20) +
  ggtitle("Degree Distribution (Cancer Co-Expression Network)") +
  xlab("Degree") +
  ylab("Frequency") +
  theme_minimal())
## Scale-free

## Compute Cancer Network
net_c <- network(adj_mat_c, matrix.type = "adjacency",
                 ignore.eval = FALSE, names.eval = "weights",
                 directed = FALSE)
## Compute density
network.density(net_c)
## Giant component
nrow(component.largest(net_c, result = "graph")) ## 483

## How many positive/negative correlations?
sum(adj_mat_c == 1)  ## 4796
sum(adj_mat_c == -1) ## 32

## Plot Cancer hubs
plot_graph(net_c, hubs_c$code_hubs,
           title = "Cancer Co-Expression Network Hubs")

## ----------------------------------------------------

# 3.3 - Analysis ----------------------------------------------------------
## Normal
## Adjacency Matrix Normal (Pearson)
adj_mat_n <- compute_adjacency(norm_coexp_N, cor_type = "pearson")

## Normal hubs
hubs_n <- compute_hubs(adj_mat_n)

## Plot degree distribution in order to check if the network
## is scale free
df2           <- data.frame(cbind(hubs_n$degree))
colnames(df2) <- "Degree"
(hist_normal <- ggplot(df2, aes(x = Degree)) +
  geom_histogram(fill = "blue", alpha = 0.7, bins = 20) +
  ggtitle("Degree Distribution (Normal Co-Expression Network)") +
  xlab("Degree") +
  ylab("Frequency") +
  theme_minimal())
## Scale-free

## Compute Cancer Network
net_n <- network(adj_mat_n, matrix.type = "adjacency",
                 ignore.eval = FALSE, names.eval = "weights",
                 directed = FALSE)
## Compute density
network.density(net_n)
## Giant component
nrow(component.largest(net_n, result = "graph")) ## 296

## How many positive/negative correlations?
sum(adj_mat_n == 1)  ## 4912
sum(adj_mat_n == -1) ## 62

## Plot Cancer hubs
plot_graph(net_n, hubs_n$code_hubs,
           title = "Normal Co-Expression Network Hubs")

## Comparison between cancer and non-cancer hubs
intersect(hubs_c$final_hubs, hubs_n$final_hubs)

## Identify the hubs characterizing only cancer tissue
hubs_c$final_hubs


## At this point, we make the same operations using a different
## correlation measure ---> Spearman

## Cancer (Spearman)
## Adjacency
adj_mat_c_sp <- compute_adjacency(norm_coexp_C, cor_type = "spearman")

## Hubs
hubs_c_sp <- compute_hubs(adj_mat_c_sp)

## Plot degree distribution in order to check if the network
## is scale free
df3           <- data.frame(cbind(hubs_c_sp$degree))
colnames(df3) <- "Degree"
ggplot(df3, aes(x = Degree)) +
  geom_histogram(fill = "blue", alpha = 0.7, bins = 20) +
  ggtitle("Degree Distribution (Cancer Co-Expression Network)",
          subtitle = "using Spearman Correlation") +
  xlab("Degree") +
  ylab("Frequency") +
  theme_minimal()
## Scale-free

## Compute Cancer Network
net_c_sp <- network(adj_mat_c_sp, matrix.type = "adjacency",
                    ignore.eval = FALSE, names.eval = "weights",
                    directed = FALSE)
## Compute density
network.density(net_c_sp)
## Giant component
nrow(component.largest(net_c_sp, result = "graph")) ## 454

## How many positive/negative correlations?
sum(adj_mat_c_sp == 1)  ## 7978
sum(adj_mat_c_sp == -1) ## 1074

## Plot Cancer hubs
plot_graph(net_c_sp, hubs_c_sp$code_hubs,
           title = "Cancer Co-Expression Network Hubs",
           subtitle = "using Spearman Correlation")

## -----------------------------------------------------

## Normal (Spearman)
## Adjacency
adj_mat_n_sp <- compute_adjacency(norm_coexp_N, cor_type = "spearman")

## Hubs
hubs_n_sp <- compute_hubs(adj_mat_n_sp)

## Plot degree distribution in order to check if the network
## is scale free
df4           <- data.frame(cbind(hubs_n_sp$degree))
colnames(df4) <- "Degree"
ggplot(df4, aes(x = Degree)) +
  geom_histogram(fill = "blue", alpha = 0.7, bins = 20) +
  ggtitle("Degree Distribution (Normal Co-Expression Network)",
          subtitle = "using Spearman Correlation") +
  xlab("Degree") +
  ylab("Frequency") +
  theme_minimal()
## Scale-free

## Compute Normal Network
net_n_sp <- network(adj_mat_n_sp, matrix.type = "adjacency",
                    ignore.eval = FALSE, names.eval = "weights",
                    directed = FALSE)
## Compute density
network.density(net_n_sp)
## Giant component
nrow(component.largest(net_n_sp, result = "graph")) ## 178

## How many positive/negative correlations?
sum(adj_mat_n_sp == 1)  ## 1898
sum(adj_mat_n_sp == -1) ## 74

## Plot Cancer hubs
plot_graph(net_n_sp, hubs_n_sp$code_hubs,
           title = "Normal Co-Expression Network Hubs",
           subtitle = "using Spearman Correlation")

## Compare the hubs
## Cancer - Cancer (Pearson Vs. Spearman)
intersect(hubs_c$final_hubs, hubs_c_sp$final_hubs)
## Normal - Normal (Pearson Vs. Spearman)
intersect(hubs_n$final_hubs, hubs_n_sp$final_hubs)
## Cancer - Normal (Spearman)
intersect(hubs_c_sp$final_hubs, hubs_n_sp$final_hubs)
## Cancer - Normal (Pearson)
intersect(hubs_c$final_hubs, hubs_n$final_hubs)

## Not common
subset(hubs_c_sp$final_hubs, !hubs_c_sp$final_hubs %in% hubs_n_sp$final_hubs)

## At this point, we want to try with a different approach,
## using soft-t.
## Cancer network

# 3.4 Bonus Point ---------------------------------------------------------

## Implement the functions used to compute the soft thresholding
soft_thresholding <- function(data, h_cutoff = 0.90){
  ## Input:  Transpose data, cut off for R^2
  ## Output: Plots[Scale-free topology, Mean connectivity]
  
  ## Choose a set of soft-thresholding powers
  powers <- c(1:20)
  ## Call the network topology analysis function
  sft <- pickSoftThreshold(t(data), powerVector = powers,
                           verbose = 5)
  
  ## Plot the results
  par(mfrow = c(1,2));
  cex1 = 0.9;
  plot(sft$fitIndices[,1],
       -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab = "Soft Threshold (power)",
       ylab = "Scale Free Topology Model Fit,signed R^2",
       type = "n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1],
       -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels = powers, cex = cex1, col = "red");
  
  ## This line corresponds to using an R^2 cut-off of h
  abline(h = h_cutoff, col = "red")
  
  ## Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab = "Soft Threshold (power)",
       ylab = "Mean Connectivity", type = "n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5],
       labels = powers, cex = cex1, col = "red")
}

soft_adjacency <- function(data, beta){
 ## Input:  Data, Beta
 ## Output: Plots[Degree distribution, Scale Free]
  
  ## Correlation matrix
  cor_mat       <- cor(t(data), method = "pearson")
  diag(cor_mat) <- 0
  
  ## Soft adjacency
  soft_adj <- abs(cor_mat)^beta
  
  ## Compute degree
  connectivity <- apply(soft_adj, 2, sum)
  
  hist(connectivity, main = "Distribution of Connectivity",
       xlab = "Connectiity", col = "blue")

  
  ## Scale Free Plot
  scaleFreePlot(connectivity)
}

## Cancer
## 1.
soft_thresholding(norm_coexp_C)
## Beta = 11, R^2 = 0.92, slope = -1.890

## 2.
soft_adjacency(norm_coexp_C, 11)

## Normal
## 1.
soft_thresholding(norm_coexp_N)
## Beta = 8, R^2 = 0.91, slope = -1.30
## 2.
soft_adjacency(norm_coexp_N, 8)

# 3.5 - Bonus Point ----------------------------------------------------------

## Cancer Network
## At this point, we compute the hubs of the network using different
## Centrality Index (CI) measures and we make a comparison with the
## genes obtained in the previous steps.

## Compute

## Implement a function to compute different CI measures
## (Betwness, Closeness, Eigenvector)
CI_measures <- function(net, CI_type, mart = mart_names){
  ## Input:  Network, CI measure, Mart
  ## Output: Top 5% nodes with highest CI

  CI         <- CI_type(net, gmode = "graph")
  CI_idx     <- (net %v% "vertex.names")[order(CI, decreasing = TRUE)]
  CI_names_5 <- head(CI_idx, floor(0.05 * length(CI_idx)))
  ## Let's get the hubs' names
  hubs_clean_CI <- (str_split(CI_names_5, "[.]", simplify = TRUE))[,1]
  hubs_names    <- getBM(filters = "ensembl_gene_id",
                                   attributes = c("ensembl_gene_id",
                                                  "hgnc_symbol"),
                                   values = hubs_clean_CI,
                                   mart = mart)
  ## Return
  return(hubs_names[,2])
}

## Betweenness
betweenness_top_5_c <- CI_measures(net_c, betweenness)
## Closeness
closeness_top_5_c   <- CI_measures(net_c, closeness)
## Eigenvector
eigenvector_top_5_c <- CI_measures(net_c, evcent)


## Intersection

## Betwness
cat(paste(intersect(hubs_c$final_hubs, betweenness_top_5_c),
          collapse = "\n"))
## Closeness
cat(paste(intersect(hubs_c$final_hubs, closeness_top_5_c),
          collapse = "\n"))
## Eigenvector
cat(paste(intersect(hubs_c$final_hubs, eigenvector_top_5_c),
          collapse = "\n"))
## --------------------------------------------------------

## Normal Network
## Betweenness
betweenness_top_5_n <- CI_measures(net_n, betweenness)
## Closeness
closeness_top_5_n   <- CI_measures(net_n, closeness)
## Eigenvector
eigenvector_top_5_n <- CI_measures(net_n, evcent)


## Intersection

## Betwness
cat(paste(intersect(hubs_n$final_hubs, betweenness_top_5_n),
          collapse = "\n"))
## Closeness
cat(paste(intersect(hubs_n$final_hubs, closeness_top_5_n),
          collapse = "\n"))
## Eigenvector
cat(paste(intersect(hubs_n$final_hubs, eigenvector_top_5_n),
          collapse = "\n"))



# 4. Differential Co-expressed Network ------------------------------------
## Computation

## Using normalize data
## Check
dim(norm_coexp_C)
dim(norm_coexp_N)

## Compute correlations
## Cancer data
co_net_corr_C <- cor(t(norm_coexp_C), method = "pearson")
diag(co_net_corr_C) <- 0

## Normal data
co_net_corr_N <- cor(t(norm_coexp_N), method = "pearson")
diag(co_net_corr_N) <- 0

## Calculation of differential correlations
## Apply Fisher z-transformation
z_c <- 0.5 * log((1 + co_net_corr_C) / (1 - co_net_corr_C))
z_n <- 0.5 * log((1 + co_net_corr_N) / (1 - co_net_corr_N))

## Sample size for each of the condition
n_c <- ncol(norm_coexp_C)
n_n <- ncol(norm_coexp_N)

## Z-score to evaluate the correlation
Z <- (z_c - z_n) / sqrt(1/(n_c - 3) + (1/(n_n - 3)))

## Threshold --> how do we find the optimal one?
t <- 3

## Adjacency Matrix a_ij = 0, if |Z| < 3.
adj_differential <- ifelse(abs(Z) < t, 0, 1)

## Generate network
net_diff_coex <- graph_from_adjacency_matrix(adj_differential,
                                             mode = "undirected",
                                             diag = FALSE)

## Plot network
ggnet2(net_diff_coex, color = "blue", alpha = 0.7, size = 2,
       edge.color = "grey", edge.alpha = 1, edge.size = 0.15) + 
  guides(size = "none") +
  ggtitle("Differential Co-Expression Network") +
  theme_minimal() +
  xlab("") + ylab("")

## Analysis
## Compute the degree index
degree_diff_coex <- sort(rowSums(adj_differential),
                         decreasing = TRUE)

## Who is the most connected hub?
(hub_most_connected <- degree_diff_coex[1])
## Extract the neighbours
idx   <- which(colnames(adj_differential) == names(hub_most_connected))
neigh <- names(which(adj_differential[,idx] == 1))
neigh <- c(neigh, names(hub_most_connected))

## Plot most connected hub with neighbours
subgraph_hubs <- subgraph(net_diff_coex, neigh)
V(subgraph_hubs)$color <- ifelse(V(subgraph_hubs)$name == names(hub_most_connected),
                                 "red", "blue")

## Plot network 
plot(subgraph_hubs, color = V(subgraph_hubs)$color,
     vertex.label = NA,
     vertex.size = 5)

## Plot degree distribution to check if the graph is a scale-free
df5           <- data.frame(cbind(degree_diff_coex))
colnames(df5) <- "Degree"
ggplot(df5, aes(x = Degree)) +
  geom_histogram(fill = "blue", alpha = 0.7, bins = 20) +
  ggtitle("Degree Distribution (Differential Co-Expression Network)",
          subtitle = "using t = 3") +
  xlab("Degree") +
  ylab("Frequency") +
  theme_minimal()

## Compute hubs
## how big is the degree of the most connected nodes?
(q_diff_coex <- quantile(degree_diff_coex[degree_diff_coex > 0],
                         0.95))
## The 5% of the most connected nodes have a degree greater than 28.65


## Find the hubs (5% of the nodes with highest degree values)
hubs_diff_coex <- degree_diff_coex[degree_diff_coex >= q_diff_coex]
## Genes
names(hubs_diff_coex)
## How many?
length(hubs_diff_coex)

## Get names
hubs_coexp_clean  <- (str_split(names(hubs_diff_coex),
                                "[.]", simplify = TRUE))[,1]
hubs_coexp_names  <- getBM(filters = "ensembl_gene_id",
                           attributes = c("ensembl_gene_id",
                                          "hgnc_symbol"),
                           values = hubs_coexp_clean,
                           mart = mart_names)

## Common with cancer
cat(intersect(hubs_coexp_names[,2], hubs_c$final_hubs),
    collapse = '\n')
## Common with Normal
cat(intersect(hubs_coexp_names[,2], hubs_n$final_hubs),
    collapse = '\n')


## *******************************************************

# 5. Patient Similarity Network (PSN) -------------------------------------

## Cancer Network
## Compute cosine similarity
dist_eucl <- as.matrix(1 / (1+dist(t(norm_coexp_C),
                                   method = "euclidean")))
dev.off()
x11()
## Heatmap
heatmap.2(dist_eucl, main = "Cancer Patients (Euclidean Distance)")
diag(dist_eucl) <- 0

## We implement the network
psn_network <- graph.adjacency(dist_eucl,
                               mode = "undirected",
                               weighted = TRUE)

## Louvian Algorithm for Community Detection
lc <- cluster_louvain(psn_network)
## We can see the communities
communities(lc)

## Create new attributes
V(psn_network)$community <- lc$membership
rain <- c("red", "blue", "green")
V(psn_network)$color <- rain[V(psn_network)$community]
## Set edges colors
E(psn_network)$color <- apply(as.data.frame(get.edgelist(psn_network)), 1, 
                              function(x) ifelse((V(psn_network)$community[which(lc$names == x[1])] ==
                                                    V(psn_network)$community[which(lc$names == x[2])]),
                                                 rain[V(psn_network)$community[which(lc$names == x[1])]],
                                                 "grey"))
## Plot
plot(psn_network, vertex.size = 4,
     vertex.label = NA, edge.color = E(psn_network)$color,
     main = "")

## Normal Network
## Compute cosine similarity
dist_eucl_n <- as.matrix(1 / (1+dist(t(norm_coexp_N),
                                     method = "euclidean")))
dev.off()
x11()
## Heatmap
heatmap.2(dist_eucl_n, main = "Normal Patients (Euclidean Distance)")
diag(dist_eucl_n) <- 0

## We implement the network
psn_network_n <- graph.adjacency(dist_eucl_n,
                                 mode = "undirected",
                                 weighted = TRUE)

## Louvian Algorithm for Community Detection
lc_n <- cluster_louvain(psn_network_n)
## We can see the communities
communities(lc_n)

## Create new attributes
V(psn_network_n)$community <- lc_n$membership
V(psn_network_n)$color <- rain[V(psn_network_n)$community]
## Set edges colors
E(psn_network_n)$color <- apply(as.data.frame(get.edgelist(psn_network_n)), 1, 
                                function(x) ifelse((V(psn_network_n)$community[which(lc$names == x[1])] ==
                                                      V(psn_network_n)$community[which(lc$names == x[2])]),
                                                   rain[V(psn_network_n)$community[which(lc$names == x[1])]], "grey"))
## Plot
plot(psn_network_n, vertex.size = 4,
     vertex.label = NA, edge.color = E(psn_network_n)$color,
     main = "")

## END
