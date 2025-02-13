library(Seurat)
library(spdep)
library(clusterProfiler)
library(org.Hs.eg.db) 
source("functions.R")


args <- commandArgs(trailingOnly = TRUE)

# Check if arguments are provided
if (length(args) < 3) {
    stop("Please provide 3 arguments: dataset_path, gene_panel_path and score_path")
}

# Parse the arguments
dataset_path <- args[1]
gene_panel_path <- args[2]
score_path<-args[3]

data <- readRDS(dataset_path)
patterns <- c("DeprecatedCodeword", "NegcontrolCodeword", "UnassignedCodeword","NegControlProbe","NegControlCodeword")
genes_to_remove <- grepl(paste(patterns, collapse = "|"), rownames(data))
data <- data[!genes_to_remove, ]


gene_panel <- readLines(gene_panel_path)

panel_score <- calculate_panel_score(data=data, celltype_column="celltype", panel=gene_panel)
scores_df<-as.data.frame(panel_score)

write.csv(scores_df, score_path, row.names = FALSE)

# Rscript scores.R demo_data/seurat_main_cohort.rds demo_data/AML.txt scores.csv

