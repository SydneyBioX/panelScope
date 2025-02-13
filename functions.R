binary_ct_marker_scores <- function(data, celltype_column, panel){
  
  celltypes <- unname(unlist(data[[celltype_column]]))
  unique_celltypes <- unique(unname(unlist(data[[celltype_column]])))
  n_celltypes <- length(unique_celltypes)
  counts_mtx <- GetAssayData(object = data, assay = "RNA", layer = "data")
  #counts_mtx <- GetAssay(data, "RNA")@data
  prop_ratio_matrices <- list()
  
  tmp_binary_mtx <- as.matrix(counts_mtx[rownames(counts_mtx) %in% panel,])
  tmp_binary_mtx <- ifelse(tmp_binary_mtx > 0, 1, 0)
  n_genes <- nrow(tmp_binary_mtx)
    
  prop_ratio_matrix <- matrix(NA, nrow=n_genes, ncol=n_celltypes)
  rownames(prop_ratio_matrix) <- rownames(tmp_binary_mtx)
  colnames(prop_ratio_matrix) <- unique_celltypes
    
  for (j in 1:ncol(prop_ratio_matrix)) {
      
    celltype_j <- which(celltypes == unique_celltypes[j])
    celltype_not_j <- which(celltypes != unique_celltypes[j])
      
    celltype_j_prop <- rowMeans(tmp_binary_mtx[, celltype_j]) 
    celltype_not_j_prop <- rowMeans(tmp_binary_mtx[, celltype_not_j]) 
      
    prop_ratio_matrix[,j] <- log((celltype_j_prop+1) / (celltype_not_j_prop+1))
  }

  return(prop_ratio_matrix)
}

calculate_morans <- function(ct_marker_scores){
  n_rows <- nrow(ct_marker_scores)
  n_cols <- ncol(ct_marker_scores)
    
  coords <- expand.grid(1:n_rows, 1:n_cols)
  neighbors <- dnearneigh(coords, d1=0, d2= 5, row.names = NULL)
  weights <- nb2listw(neighbors, style="S")
    
  data_vector <- as.vector(as.matrix(ct_marker_scores))
  morans_I <- moran.test(data_vector, listw = weights)
  prop_zeroes <- sum((data_vector > -0.01 & data_vector < 0.01)+.Machine$double.xmin)/length(data_vector)
  score <- 1-unname(morans_I$estimate[1]*prop_zeroes)
  
  return(score)
}

calculate_correlation_ratio <- function(data, panel){
  tmp_data <- data[rownames(data) %in% panel,]
  #cor_matrix <- cor(as.matrix(t(GetAssay(tmp_data,"RNA")@data)), method="spearman")
  cor_matrix <- cor(as.matrix(t(LayerData(tmp_data,layer = "data"))), method="spearman")
  cor_matrix <- cor_matrix[row(cor_matrix) < col(cor_matrix)]
  cor_matrix <- cor_matrix[!is.na(cor_matrix)]
  cor_ratio <- log2((sum(cor_matrix<=0)/sum(cor_matrix>0)))
    
  return(cor_ratio)
}




calculate_pathway_score <- function(data, panel, species="org.Hs.eg.db"){
  
  results_df <- data.frame(`# genes`=NA,
                           `% genes`=NA,
                           `# pathways`=NA,
                           `max qvalue`=NA,
                           `diversity score`=NA)
  gene_pathway_df <- list()
  
  ego <- enrichGO(gene=panel,
                  keyType="SYMBOL",
                  OrgDb=species,
                  ont="BP",
                  pAdjustMethod="BH",
                  pvalueCutoff=0.05,
                  qvalueCutoff=0.05,
                  maxGSSize=500,
                  readable=TRUE)
  
  ego_df <- data.frame(ego@result[ego@result$p.adjust<0.05,])
  genes <- unique(unlist(strsplit(ego_df$geneID,"/")))
  n_genes <- length(genes)
  
  ego_simplified <- clusterProfiler::simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
  ego_simplified_df <- ego_simplified@result
  
  diversity_score <- 1-((nrow(ego_df) - nrow(ego_simplified_df))/nrow(ego_df))
  
  tmp_df <- data.frame(`# genes`=n_genes,
                       `% genes`=n_genes/length(panel),
                       `# pathways`=nrow(ego_df),
                       `max qvalue`=ifelse(max(ego_df$qvalue) < 0, NA, max(ego_df$qvalue)),
                       `diversity score`=ifelse(is.na(diversity_score), NA, diversity_score))
  results_df <- rbind(results_df, tmp_df)
  
  
  gene_pathway_df <- do.call(rbind, gene_pathway_df)
  results_df[,2:ncol(results_df)] <- round(results_df[,2:ncol(results_df)], digits=3)
  colnames(results_df) <- c("# genes", "% genes", "# pathways", "max qvalue", "diversity score")
  results_df <- results_df[-1,]
  score <- results_df$`% genes`
  return(score)
}

# x_col: x coordinates of cells
# y_col: y coordinates of cells
# The higher the better
spatial_moransI <- function(data, panel, x_col, y_col){
  exp_mtx <- as.matrix(GetAssayData(object = data, assay = "RNA", slot = "data"))
  #exp_mtx <- as.matrix(GetAssay(data, "RNA")$data)
  coords <- data.frame(x=data@meta.data[[x_col]], y=data@meta.data[[y_col]])
  neighbors <- dnearneigh(as.matrix(coords), d1=0, d2= 10, row.names = NULL)
  weights <- nb2listw(neighbors, style="W", zero.policy = TRUE)
  
  panel <- panel[panel %in% rownames(exp_mtx)]
  tmp_morans_I_vec <- numeric(length(panel))
  for(j in 1:length(panel)){
    
    data_vector <- exp_mtx[panel[j],]
    m <- moran.test(data_vector, listw = weights)
    tmp_morans_I_vec[[j]] <- unname(m$estimate[1])
    
  }
  
  return(median(tmp_morans_I_vec))
}

# Run this one and pass it to the function to get the transcriptional variability score.
reference_clustering <- function(data){
  
  orig_data <- NormalizeData(data)
  orig_data <- FindVariableFeatures(orig_data, selection.method = "vst", nfeatures = 2000)
  orig_data <- ScaleData(orig_data, features = rownames(orig_data))
  orig_data <- RunPCA(orig_data, features = VariableFeatures(object = orig_data))
  orig_data <- FindNeighbors(orig_data, dims = 1:10)
  orig_data <- FindClusters(orig_data, resolution = 0.5)
  orig_clusters <- orig_data$seurat_clusters
  
  return(orig_clusters)
}

# NMI: the higher the better
transcriptional_variability <- function(data, panel, orig_clusters){
  
  tmp_data <- data[rownames(data) %in% panel,]
  tmp_data <- NormalizeData(tmp_data)
  tmp_data <- FindVariableFeatures(tmp_data, selection.method = "vst", nfeatures = nrow(tmp_data))
  tmp_data <- ScaleData(tmp_data, features = rownames(tmp_data))
  tmp_data <- RunPCA(tmp_data, features = VariableFeatures(object = tmp_data))
  tmp_data <- FindNeighbors(tmp_data, dims = 1:2)
  tmp_data <- FindClusters(tmp_data, resolution = 0.5)
  tmp_clusters <- tmp_data$seurat_clusters
  nmi <- aricode::NMI(orig_clusters, tmp_clusters)
    
  return(nmi)
}







calculate_panel_score <- function(data, celltype_column, panel){
  binary_mtx <- binary_ct_marker_scores(data=data, celltype_column=celltype_column, panel=panel)
  cts <- calculate_morans(binary_mtx)
  cat('1')
  cor_ratio <- calculate_correlation_ratio(data=data, panel=panel)
  cat('2')
  pathway <- calculate_pathway_score(data=data, panel=panel, species="org.Hs.eg.db")
  cat('3')
  spatial <- spatial_moransI(data=data, panel=panel, x_col='x_centroid', y_col='y_centroid')
  cat('4')
  ref_clu <- reference_clustering(data=data)
  tv <- transcriptional_variability(data=data, panel=panel, orig_clusters=ref_clu)
  return(list(Celltype_specificity_score = cts,  Correlation_ratio = cor_ratio, Pathway_score=pathway, Spatial_score=spatial,Transcriptional_variability=tv) )
  
}
