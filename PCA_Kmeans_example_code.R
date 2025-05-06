library(ggplot2)
library(readxl)


# --- read gct data
read_gct <- function(file) {
  lines <- readLines(file)
  
  # First line is version
  version <- lines[1]
  
  # Second line: number of rows and columns
  dims <- strsplit(lines[2], "\t")[[1]]
  n_rows <- as.integer(dims[1])
  n_cols <- as.integer(dims[2])
  
  # Read the actual data
  df <- read.delim(file, skip = 2, header = TRUE, stringsAsFactors = FALSE)
  
  # Extract expression matrix
  expr <- as.matrix(df[, -(1:2)])  # Remove Name and Description
  rownames(expr) <- df$Name
  
  return(list(
    version = version,
    data = df,
    expression = expr
  ))
}

setwd("~/Downloads")


gtex_data <- read_gct("GTEx_Analysis_2022-06-06_v10_RNASeQCv2.4.2_gene_reads.gct")
gtex_attr <- read_xlsx("GTEx_Analysis_v10_Annotations_SampleAttributesDD.xlsx")
gtex_metadata <- read.delim("GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt", sep = "\t", header = TRUE)
gtex_age <- read.delim("GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt", sep = "\t", header = TRUE)

gtex_age$SUBJID <- gsub("-", ".", gtex_age$SUBJID)

gtex_tissues <- unique(gtex_metadata$SMTS)

tissue_IDs<- function(tissue) {
  IDs <- gsub("-", ".", gtex_metadata$SAMPID[gtex_metadata$SMTS==tissue])
  tissue <- gsub(" ", "_", tissue)
  var_name <- paste(tissue, "IDs", sep="_")
  assign(var_name, IDs, envir=.GlobalEnv)
  return(var_name)
}

tissue_list <- NULL
df_name <- NULL

for (i in gtex_tissues) {
  tissue_type <- tissue_IDs(i)
  tissue_list <- c(tissue_list, tissue_type)
  iterable_tissue <- get(tissue_type)
  
  assign(tissue_type, iterable_tissue, envir=.GlobalEnv)
}

all_tissues <- c()
for (i in tissue_list) {
  all_tissues <- c(all_tissues, get(i))
}

tissue_df <- data.frame()
for (i in tissue_list) {
  tissue_type <- gsub("_IDs", "", i)
  tissue_type <- gsub("_", " ", tissue_type)
  tissue_df <- rbind(tissue_df, data.frame(tissue=tissue_type, sampleID = get(i)))
}

tissue_df <- tissue_df |> mutate(SUBJID = sub("^([^.]+\\.[^.]+)\\..*", "\\1", sampleID))

ages <- gtex_age |> select(SUBJID, AGE)

age_tissues <- merge(tissue_df, ages, by = "SUBJID")

age_limiter <- function(age) {
  as.list(age_tissues[age_tissues$AGE == age, ]$sampleID)
}

age_seventy <- age_limiter("70-79")
age_sixty <- age_limiter("60-69")
age_fifty <- age_limiter("50-59")
age_forty <- age_limiter("40-49")
age_thirty <-age_limiter("30-39")
age_twenty <- age_limiter("20-29")


gene_gatherer <- function(age_list) {
  gtex_data$expression[,colnames(gtex_data$expression) %in% age_list]
}

seventy_expr <- gene_gatherer(age_seventy)
sixty_expr <- gene_gatherer(age_sixty)
fifty_expr <- gene_gatherer(age_fifty)
forty_expr <-gene_gatherer(age_forty)
thirty_expr <- gene_gatherer(age_thirty)
twenty_expr <- gene_gatherer(age_twenty)

popular_genes <- function(gene_age_expr) {
  gene_age_expr[rowSums(gene_age_expr > 20) >= 20, ]
}

seventy_expr_filtered <- popular_genes(seventy_expr)
sixty_expr_filtered <- popular_genes(sixty_expr)
fifty_expr_filtered <- popular_genes(fifty_expr)
forty_expr_filtered <- popular_genes(forty_expr)
thirty_expr_filtered <- popular_genes(thirty_expr)
twenty_expr_filtered <- popular_genes(twenty_expr)

PCA_graph <- function(filtered_expr, total_age_tissues) {
  gene_vars <- apply(filtered_expr, 1, var)
  top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:500])  # top 500 most variable genes
  filtered_expr <- filtered_expr[top_genes, ]
  expr <- filtered_expr
  expr_t <- t(expr)
  pca_result <- prcomp(expr_t, center = TRUE, scale. = TRUE)
  pca_df <- as.data.frame(pca_result$x[, 1:2])
  pca_df$Sample <- rownames(pca_df)
  pca_tissues <- total_age_tissues |> select(sampleID, tissue)
  merger <- left_join(pca_df, pca_tissues, by = c("Sample" = "sampleID"))
  pca_df$Tissue <- merger$tissue
  
  return_list <- list(expr_t, pca_df)
  
  return(return_list)
}

seventy_PCA <- PCA_graph(seventy_expr_filtered, age_tissues)
sixty_PCA <- PCA_graph(sixty_expr_filtered, age_tissues)
fifty_PCA <- PCA_graph(fifty_expr_filtered, age_tissues)
forty_PCA <- PCA_graph(forty_expr_filtered, age_tissues)
thirty_PCA <- PCA_graph(thirty_expr_filtered, age_tissues)
twenty_PCA <- PCA_graph(twenty_expr_filtered, age_tissues)


ggplot(seventy_PCA[[2]], aes(x = PC1, y = PC2, color=Tissue)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA: First Two Principal Components", x = "PC1", y = "PC2")

ggplot(sixty_PCA[[2]], aes(x = PC1, y = PC2, color=Tissue)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA: First Two Principal Components", x = "PC1", y = "PC2")

ggplot(fifty_PCA[[2]], aes(x = PC1, y = PC2, color=Tissue)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA: First Two Principal Components", x = "PC1", y = "PC2")

ggplot(forty_PCA[[2]], aes(x = PC1, y = PC2, color=Tissue)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA: First Two Principal Components", x = "PC1", y = "PC2")

ggplot(thirty_PCA[[2]], aes(x = PC1, y = PC2, color=Tissue)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA: First Two Principal Components", x = "PC1", y = "PC2")

ggplot(twenty_PCA[[2]], aes(x = PC1, y = PC2, color=Tissue)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA: First Two Principal Components", x = "PC1", y = "PC2")


K_means <- function(expr_t, pca_df) {
  expr_scaled <- scale(expr_t)
  kmeans_result <- kmeans(expr_scaled, centers = 2)
  cluster_df <- data.frame(
    Sample = rownames(expr_scaled),
    Cluster = as.factor(kmeans_result$cluster)
  )
  
  pca_df$Cluster_Kmean <- cluster_df$Cluster
  
  return(pca_df)
}

seventy_K_means <- K_means(seventy_PCA[[1]], seventy_PCA[[2]])
sixty_K_means <- K_means(sixty_PCA[[1]], sixty_PCA[[2]])
fifty_K_means <- K_means(fifty_PCA[[1]], fifty_PCA[[2]])
forty_K_means <- K_means(forty_PCA[[1]], forty_PCA[[2]])
thirty_K_means <- K_means(thirty_PCA[[1]], thirty_PCA[[2]])
twenty_K_means <- K_means(twenty_PCA[[1]], twenty_PCA[[2]])

ggplot(seventy_K_means, aes(x = PC1, y = PC2, color = Cluster_Kmean)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "K-means Clustering (k=2) on Expression Data")

ggplot(sixty_K_means, aes(x = PC1, y = PC2, color = Cluster_Kmean)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "K-means Clustering (k=2) on Expression Data")

ggplot(fifty_K_means, aes(x = PC1, y = PC2, color = Cluster_Kmean)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "K-means Clustering (k=2) on Expression Data")

ggplot(forty_K_means, aes(x = PC1, y = PC2, color = Cluster_Kmean)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "K-means Clustering (k=2) on Expression Data")

ggplot(thirty_K_means, aes(x = PC1, y = PC2, color = Cluster_Kmean)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "K-means Clustering (k=2) on Expression Data")

ggplot(twenty_K_means, aes(x = PC1, y = PC2, color = Cluster_Kmean)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "K-means Clustering (k=2) on Expression Data")



