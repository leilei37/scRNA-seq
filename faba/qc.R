library(stringr)
library(dplyr)
suppressPackageStartupMessages(library(Seurat))
library(Matrix)
library(patchwork)
library(ggplot2)
library(cowplot)
library(purrr)
library(tibble)
#devtools::install_github('satijalab/seurat-wrappers')
suppressMessages(require(SeuratWrappers))
library(metap)
library(multtest)
library(plotly)
library(scater)
library(DESeq2)
library(RColorBrewer)
library(viridis)
library(VennDiagram)
library(UpSetR)
library(ComplexUpset)
packageVersion('Seurat')
anno<-read.table('~/Downloads/genome/HedinV2.mercator.functional.txt',sep = '\t',header = T)

p7s2 <- Read10X("~/Downloads/legume/faba/faba_iseq_filter_cr/") 
p7s2 <- Read10X("~/Downloads/legume/faba/5sampleiseq/dap20-2filter_cr/") 
p7s2 <- Read10X("~/Downloads/legume/faba/5sampleiseq/dap25-1filter_cr/") 
p7s2 <- Read10X("~/Downloads/legume/faba/5sampleiseq/dap30-1filter_cr/") 
p7s2 <- Read10X("~/Downloads/legume/faba/5sampleiseq/dpi7filter_cr/") 
p7s2 <- Read10X("~/Downloads/legume/faba/5sampleiseq/root7filter_cr/") 

p7s2 <- Read10X("~/Downloads/legume/faba/round2_nucleu/36hair1_filter_cr/") 
p7s2 <- Read10X("~/Downloads/legume/faba/round2_nucleu/36hair2_filter_cr/") 
p7s2 <- Read10X("~/Downloads/legume/faba/round2_nucleu/d15r1_filter_cr/") 
p7s2 <- Read10X("~/Downloads/legume/faba/round2_nucleu/") 
p7s2 <- Read10X("~/Downloads/legume/faba/round2_nucleu/") 
p7s2 <- Read10X("~/Downloads/legume/faba/round2_nucleu/") 


p7s2 <- CreateSeuratObject(counts = p7s2, project = "d15_r1", min.cells = 3, min.features = 200)
p7s2 <- AddMetaData(object = p7s2, metadata = "d15_r1", col.name = "Condition")

p7s2[["percent.mt"]] <- PercentageFeatureSet(p7s2, pattern = "LD")

#gene.cp<-read.table('~/Downloads/legume/faba/chloroplast.name.txt')[,1]
#gene.cp<-intersect(gene.cp,rownames(p7s2))
#p7s2[["percent.chloroplast"]] <- PercentageFeatureSet(p7s2, features = gene.cp)
p7s2[["percent.chloroplast"]] <- PercentageFeatureSet(p7s2, pattern = "nbis")

#gene.ssp<-gsub('\\.\\d+$','',anno$IDENTIFIER[grep('seed storage protein',anno$mercator4v6.0)])
#gene.ssp<-intersect(gene.ssp,rownames(p7s2))
#p7s2[["percent.ssp"]] <- PercentageFeatureSet(p7s2, features = gene.ssp)

library(DropletUtils)
counts_matrix <- GetAssayData(p7s2, slot = "counts")
if (!inherits(counts_matrix, "dgCMatrix")) {
  counts_matrix <- as(counts_matrix, "dgCMatrix")
}
bcrank <- barcodeRanks(counts_matrix)
# Only showing unique points for plotting speed.
uniq <- !duplicated(bcrank$rank)
plot(bcrank$rank[uniq], bcrank$total[uniq], log="xy",
     xlab="Rank", ylab="Total UMI count", cex.lab=1.2)
abline(h=metadata(bcrank)$inflection, col="darkgreen", lty=2)
abline(h=metadata(bcrank)$knee, col="dodgerblue", lty=2)
hist(p7s2$nCount_RNA,labels = T,breaks = 1000,xlim = c(0,5000))
hist(p7s2$nFeature_RNA,labels = T,breaks = 500,xlim = c(0,2000))
median(p7s2$nCount_RNA)
median(p7s2$nFeature_RNA)

VlnPlot(p7s2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.chloroplast"), ncol = 5)
plot1 <- FeatureScatter(p7s2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(p7s2, feature1 = "nCount_RNA", feature2 = "percent.chloroplast")
plot1+plot2
FeatureScatter(p7s2,feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#detected doublet cells
library(scDblFinder)
panc.sce <- as.SingleCellExperiment(p7s2)
sce <- scDblFinder(panc.sce, samples="Condition",clusters=T)
p7s2@meta.data[["scDblFinder.class"]]<-sce$scDblFinder.class
table(sce$scDblFinder.class)
#table(p7s2_filt$scDblFinder.class,p7s2_filt$not_ambient)

#filter doublet cells
p7s2_sin <- subset(p7s2, subset = scDblFinder.class  == "singlet")
#p7s2_noa <- subset(p7s2_sin, subset = not_ambient  == "TRUE")

VlnPlot(p7s2_sin, features = c("nFeature_RNA", "nCount_RNA","percent.mt", "percent.chloroplast"), ncol = 4)
plot1 <- FeatureScatter(p7s2_sin, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(p7s2_sin, feature1 = "nCount_RNA", feature2 = "percent.chloroplast")
plot1+plot2
FeatureScatter(p7s2_sin,feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#filter cells based on percent.mt & percent.chloroplast
#p7s2_filt <- subset(p7s2_sin, subset = nCount_RNA > 200 & nFeature_RNA > 200 &
p7s2_filt <- subset(p7s2_sin, subset = 
                      nFeature_RNA < 15000 & nCount_RNA<20000 &
                     percent.mt < 0.3 & percent.chloroplast < 0.1)
                     # percent.mt < 20 & percent.chloroplast < 10)

VlnPlot(p7s2_filt, features = c("nFeature_RNA", "nCount_RNA","percent.mt", "percent.chloroplast"), ncol = 4)
plot1 <- FeatureScatter(p7s2_filt, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(p7s2_filt, feature1 = "nCount_RNA", feature2 = "percent.chloroplast")
plot1+plot2
FeatureScatter(p7s2_filt,feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
saveRDS(p7s2_filt,file = '~/Downloads/legume/faba/round2_nucleu/faba_36hai_r2_filt.rds')
readRDS('~/Downloads/legume/faba/round2_nucleu/faba_36hai_r1_filt.rds')->p7s1_filt
readRDS('~/Downloads/legume/faba/round2_nucleu/faba_36hai_r2_filt.rds')->p7s2_filt

#Integrating using SCtransform and reference
faba.list <- list(p7s1_filt, p7s2_filt)
#faba.list <- list(p7s1, p7s2)

faba.list <- lapply(X = faba.list, FUN = function(x) {
#  x <- SCTransform(x, vars.to.regress = c('nCount_RNA','nFeature_RNA'), verbose = TRUE)
   x <- SCTransform(x, vars.to.regress = c("percent.mt",  "percent.chloroplast",'nCount_RNA'), verbose = TRUE)
#   x <- SCTransform(x, vars.to.regress = 'percent.ssp',verbose = TRUE)
 })
faba.features <- SelectIntegrationFeatures(object.list = faba.list, nfeatures = 3000)
faba.list <- PrepSCTIntegration(object.list = faba.list, anchor.features = faba.features)
faba.anchors <- FindIntegrationAnchors(object.list = faba.list, normalization.method = "SCT",
                                       anchor.features = faba.features, reference = c(1,2))
faba.integrated <- IntegrateData(anchorset = faba.anchors, normalization.method = "SCT")

#library(harmony)
#faba.integrated <- RunHarmony(
#  faba.integrated,
#  group.by.vars = "Condition",
#  theta = 2, 
#  lambda = 0.5,
#  assay.use = "SCT"
#)
faba.integrated <- RunUMAP(faba.integrated, assay.use = "SCT",  reduction = "harmony", dims = 1:50)


faba.integrated<-PrepSCTFindMarkers(faba.integrated)
markers <- FindAllMarkers(faba.integrated,min.diff.pct = 0.2, logfc.threshold = 0.5,only.pos = TRUE)
#markertop10<-markers
markertop10<-markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
#pdf('Downloads/label_medicago96r2_top5.pdf')
DotPlot(faba.integrated,features=unique(c(markertop10$gene)))+coord_flip()
FeaturePlot(faba.integrated,'Vfaba.Hedin2.R2.1g000768',label = T)
saveRDS(faba.integrated,file = '~/Downloads/legume/faba/round2_nucleu/faba_36hai_r2_filt.rds')
