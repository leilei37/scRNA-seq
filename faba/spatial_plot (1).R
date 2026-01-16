library(ggplot2)
library(cowplot)
library(Seurat)
library(RColorBrewer)

load('L5_intissue_seurat.RData')
plot_clusters <- function(seurat_obj, clusters, reduction = "umap") {
  meta_df <- seurat_obj@meta.data
  
  # Figure 1: Full-color spatial
  p1 <- ggplot(meta_df, aes(x = imagecol, y = imagerow, color = seurat_clusters)) +
    geom_point(size = 0.5) +
    coord_fixed() +
    theme_minimal() +theme(legend.position = "none") +
    labs(title = "Spatial Distribution (All Clusters)",
         x = "X coordinate", y = "Y coordinate")
  
#Figure 2: Highlight the spatial of the specified cluster
  p2 <- ggplot(meta_df, aes(x = imagecol, y = imagerow)) +
    geom_point(color = "lightgrey", size = 0.5) +
    geom_point(data = subset(meta_df, seurat_clusters %in% clusters),
               aes(color = seurat_clusters), size = 0.6) +
    scale_color_manual(values = setNames(brewer.pal(max(3, length(clusters)), "Set1"),
                                         clusters)) +
    coord_fixed() +
    theme_minimal() +
    labs(title = paste0("Spatial Distribution (Highlighted: ",
                        paste(clusters, collapse = ", "), ")"),
         x = "X coordinate", y = "Y coordinate")
  
  # Figure 3: Full-color DimPlot
  p3 <- DimPlot(seurat_obj, reduction = reduction, group.by = "seurat_clusters",label = T) +
    ggtitle(paste("DimPlot (All Clusters, ", toupper(reduction), ")", sep = ""))
  
# Figure 4: Highlight the DimPlot of the specified cluster
  p4 <- DimPlot(seurat_obj, reduction = reduction, group.by = "seurat_clusters",
                cells.highlight = WhichCells(seurat_obj, idents = clusters)) +
    scale_color_manual(values = c("grey", "red")) +
    ggtitle(paste("DimPlot (Highlighted: ",
                  paste(clusters, collapse = ", "), ")", sep = ""))
  
  # 2×2 
  plot_grid(p1, p2, p3, p4, ncol = 2)
}

# example
plot_clusters(seurat_obj, c("0"), reduction = "umap")
plot_clusters(seu, c("5",'7'), reduction = "umap")
seu@meta.data[["imagecol"]]<-seu@meta.data[["x"]]
seu@meta.data[["imagerow"]]<-seu@meta.data[["y"]]

library(ggplot2)
library(cowplot)
library(Seurat)

plot_gene <- function(seurat_obj, gene, reduction = "umap") {
  # 提取基因表达
  expr <- FetchData(seurat_obj, vars = gene)
  meta_df <- cbind(seurat_obj@meta.data, expr)
  colnames(meta_df)[ncol(meta_df)] <- "expr"
  
  p1 <- ggplot(meta_df, aes(x = imagecol, y = imagerow, color = seurat_clusters)) +
    geom_point(size = 0.5) +
    coord_fixed() +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(title = "Spatial Distribution (All Clusters)",
         x = "X coordinate", y = "Y coordinate")
  
  p2 <- ggplot(meta_df, aes(x = imagecol, y = imagerow, color = expr)) +
    geom_point(size = 0.5) +
    scale_color_gradient(low = "lightgrey", high = "red") +
    coord_fixed() +
    theme_minimal() +
    labs(title = paste0("Spatial Expression (", gene, ")"),
         x = "X coordinate", y = "Y coordinate", color = "Expr")
  
  p3 <- DimPlot(seurat_obj, reduction = reduction, group.by = "seurat_clusters",label = T) +
    ggtitle(paste("DimPlot (All Clusters, ", toupper(reduction), ")", sep = ""))
  
  p4 <- FeaturePlot(seurat_obj, features = gene, reduction = reduction) +
    ggtitle(paste("DimPlot Expression (", gene, ")"))
  
  plot_grid(p1, p2, p3, p4, ncol = 2)
  #p4
}

plot_gene(seurat_obj, "Vfaba.Hedin2.R2.1g010869", reduction = "umap")
plot_gene(seurat_obj, "Vfaba.Hedin2.R2.1g010655", reduction = "umap")
plot_gene(seu, "Vfaba.Hedin2.R2.1g010869", reduction = "umap")


