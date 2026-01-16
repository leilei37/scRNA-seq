# # single-cell analysis package
# library(Seurat)
# 
# # plotting and data science packages
# library(tidyverse)
# library(cowplot)
# library(patchwork)
# 
# # co-expression network analysis packages:
# library(WGCNA)
# library(hdWGCNA)
# 
# # using the cowplot theme for ggplot
# theme_set(theme_cowplot())
# 
# load('~/Downloads/sc/sp/L5_intissue_seurat.RData')
# #get metacell matrix
# DefaultAssay(seurat_obj) <- "RNA"
# 
# seurat_obj <- SetupForWGCNA(
#   seurat_obj,
#   gene_select = "fraction", # the gene selection approach
#   fraction = 0.0005, # fraction of cells that a gene needs to be expressed in order to be included
#   wgcna_name = "faba" # the name of the hdWGCNA experiment
# )
# seurat_obj <- MetacellsByGroups(
#   seurat_obj = seurat_obj,
#   group.by = c('seurat_clusters'), # specify the columns in seurat_obj@meta.data to group by
#   reduction = 'umap', # select the dimensionality reduction to perform KNN on
#   k = 25, # nearest-neighbors parameter
#   max_shared = 10, # maximum number of shared cells between two metacells
#   ident.group = 'seurat_clusters' # set the Idents of the metacell seurat object
# )
# #seurat_obj@misc[["faba"]][["wgcna_genes"]]
# 
# seurat_obj <- NormalizeMetacells(seurat_obj)
# seurat_obj <- SetDatExpr(
#   seurat_obj,
#   group_name = seurat_obj@misc[["faba"]][["wgcna_metacell_obj"]]@meta.data[["seurat_clusters"]],
#   group.by = "seurat_clusters"
# )
# exp.sc.all <- GetDatExpr(seurat_obj)
# 
# save(exp.sc.all, file = "~/Downloads/sc/sp/sp_metacell_expression.RData")

load('sp_metacell_expression.RData')
library(data.table)
#calculate top 500 correlations
calculate_top500_correlations <- function(expression_data) {
  gene_list <- colnames(expression_data)
  
  # **1. Calculate the entire correlation matrix**
  cor_matrix <- cor(expression_data, method = "pearson")
  
  # **2. Convert to long format **
  cor_df <- as.data.table(melt(cor_matrix, varnames = c("gene1", "gene2"), value.name = "correlation"))
  
  # **3. Keep only combinations of different genes (remove the line gene1 == gene2)
  cor_df <- cor_df[gene1 != gene2]
  
  # **4. Calculate P-values (optional, speed up the calculation can not save P-values)
  cor_df[, p.value := corPvalueStudent(correlation, nrow(expression_data))]
  
  # **5. Calculate FDR**
  cor_df[, fdr := p.adjust(p.value, method = "BH")]
  
  cor_df[, correlation := round(correlation, 2)]  # The correlation retains 2 decimal places
  cor_df[, p.value := signif(p.value, 3)]  # The p-value retains 3 significant digits
  cor_df[, fdr := signif(fdr, 3)] 
  
  # **6. Get the Top 500 gene pairs for each gene **
  cor_df <- cor_df[order(gene1, -abs(correlation))]  # Sorted by gene and in descending order of relevance
  cor_df <- cor_df[, head(.SD, 500), by = gene1]  # Take the first 500 related genes of each gene
  
  return(cor_df)
}

top500_correlations <- calculate_top500_correlations(exp.sc.all)
