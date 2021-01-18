#Description: this R script provides the single cell data analysis for maunscript:Molecular characterization of selectively vulnerable neurons in Alzheimer's Disease
#Note: The manuscript provides data with 2 brain regions: entorhinal cortex (EC) and superior frontal gyrus (SFG), this script is specifically for EC, and for SFG, the structure of the code is the same, for minor parameter differences, please see Methods part in our maunscript.
#Dataset: GEO -- GSE147528


# dir_work <- ""  ##provide a dir_work 

setwd(dir_work)

dir_data = "../data/processed_data"

library(DropletUtils)
library(scran)
library(BiocSingular)
library(scater)
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(grid)
library(cowplot)

ec_astrocyte <- readRDS(file = "../data/processed_data/EC_astrocytes.rds")

clusters <- quickCluster(ec_astrocyte, BSPARAM=IrlbaParam())
table(clusters)

ec_astrocyte <- computeSumFactors(ec_astrocyte, min.mean=0.1, cluster=clusters)
summary(sizeFactors(ec_astrocyte))

ec_astrocyte <- logNormCounts(ec_astrocyte)

ec_astrocyte_seurat = as.Seurat(ec_astrocyte)

ec_astrocyte <- modelGeneVar(ec_astrocyte)
hvgtop1000.ec_astrocyte <- getTopHVGs(ec_astrocyte, n=1000)

stage0_samples = ec_astrocyte_seurat[, ec_astrocyte_seurat$BraakStage == 0]
stage2_samples = ec_astrocyte_seurat[, ec_astrocyte_seurat$BraakStage == 2]
stage6_samples = ec_astrocyte_seurat[, ec_astrocyte_seurat$BraakStage == 6]



s0_samples = colnames(stage0_samples)
s6_samples = colnames(stage6_samples)

#==================================CODE WITH BATCH EFFECT CORRECTION===========================================#
batch_a_idx = which(ec_astrocyte_seurat@meta.data$SampleBatch == 'A')
batch_a_obj = ec_astrocyte_seurat[, batch_a_idx]

batch_b_idx = which(ec_astrocyte_seurat@meta.data$SampleBatch == 'B')
batch_b_obj = ec_astrocyte_seurat[, batch_b_idx]


batch_c_idx = which(ec_astrocyte_seurat@meta.data$SampleBatch == 'C')
batch_c_obj = ec_astrocyte_seurat[, batch_c_idx]

batch_d_idx = which(ec_astrocyte_seurat@meta.data$SampleBatch == 'D')
batch_d_obj = ec_astrocyte_seurat[, batch_d_idx]


sample_list = list()
sample_list[['a']] <- batch_a_obj
sample_list[['b']] <- batch_b_obj
sample_list[['c']] <- batch_c_obj
sample_list[['d']] <- batch_d_obj


sample_anchors <- FindIntegrationAnchors(object.list = sample_list, dims = 1:30, reduction = "cca", l2.norm = TRUE)

sample_integrated <- IntegrateData(anchorset = sample_anchors, dims = 1:30)

DefaultAssay(sample_integrated) <- "integrated"

sample_integrated <- ScaleData(sample_integrated,  verbose = FALSE, features = hvgtop1000.ec_astrocyte)

sample_integrated <- RunPCA(sample_integrated, features = hvgtop1000.ec_astrocyte, npcs = 30)

sample_integrated <- FindNeighbors(sample_integrated, dims = 1:12)

sample_integrated = FindClusters(sample_integrated, resolution = 0.2)

# sample_integrated <- RunTSNE(object = sample_integrated, dims = 1:12)
# TSNEPlot(object = sample_integrated) + scale_color_manual(values=c("#4682B4", "#3CB371", "#DC143C", "#FFA500")) 

sample_integrated <- RunUMAP(object = sample_integrated, dims = 1:12)
UMAPPlot(object = sample_integrated) + scale_color_manual(values=c("#4682B4", "#3CB371", "#DC143C", "#FFA500")) 


dir_plot = '../plot/ec_umap.pdf'
pdf(dir_plot, width = 9, height = 10)
DimPlot(sample_integrated, reduction = "umap", label = FALSE, pt.size = 0.8) + scale_color_manual(values=c("#4682B4", "#3CB371", "#FF0000", "#FFA500")) + NoLegend() + ggtitle("Entorhinal Cortex (EC) Astrocytes") +
  theme(plot.title = element_text(size=25, hjust = 0.5), axis.title.x = element_text(size=22), axis.title.y = element_text(size=22), axis.text.x = element_text(size=22), axis.text.y = element_text(size=22), legend.text=element_text(size=22))
dev.off()



DefaultAssay(sample_integrated) <- "RNA"
all.genes = rownames(sample_integrated@assays$RNA@counts)




marker_features = c("GFAP", "CD44", "HSPB1", "TNC", "SLC1A2", "SLC1A3", "GLUL", "NRXN1", "CADM2", "PTN", "GPC5")
dir_plot = "../plot/ec_astrocyte_violin_plot.png"
png(dir_plot, width = 10, height = 18, units = 'in', res = 300)
p1 <- FeaturePlot(sample_integrated, features = "GFAP", slot = "data", min.cutoff = "q9") +
  theme(plot.title = element_text(size = 25, face="bold.italic"), axis.title.x = element_text(size=22), axis.title.y = element_text(size=22), axis.text.x = element_text(size=22), axis.text.y = element_text(size=22), legend.text=element_text(size=20)) 
p1 <- p1 + scale_color_gradientn(colors = c("grey", "white", "red"),limits=c(0, 4), breaks=c(0,2,4))
p2 <- FeaturePlot(sample_integrated, features = "CD44", slot = "data", min.cutoff = "q9") +
  theme(plot.title = element_text(size = 25, face="bold.italic"), axis.title.x = element_text(size=22), axis.title.y = element_text(size=22), axis.text.x = element_text(size=22), axis.text.y = element_text(size=22), legend.text=element_text(size=20))
p2 <- p2 + scale_color_gradientn(colors = c("grey", "white", "red"),limits=c(0, 3), breaks=c(0,1,2,3))
p3 <- FeaturePlot(sample_integrated, features = "HSPB1", slot = "data", min.cutoff = "q9") +
  theme(plot.title = element_text(size = 25, face="bold.italic"), axis.title.x = element_text(size=22), axis.title.y = element_text(size=22), axis.text.x = element_text(size=22), axis.text.y = element_text(size=22), legend.text=element_text(size=20))
p3 <- p3 + scale_color_gradientn(colors = c("grey", "white", "red"),limits=c(0, 4), breaks=c(0,2,4))
p4 <- FeaturePlot(sample_integrated, features = "TNC", slot = "data", min.cutoff = "q9") +
  theme(plot.title = element_text(size = 25, face="bold.italic"), axis.title.x = element_text(size=22), axis.title.y = element_text(size=22), axis.text.x = element_text(size=22), axis.text.y = element_text(size=22), legend.text=element_text(size=20))
p4 <- p4 + scale_color_gradientn(colors = c("grey", "white", "red"),limits=c(0, 3), breaks=c(0,1,2,3))
p5 <- FeaturePlot(sample_integrated, features = "SLC1A2", slot = "data", min.cutoff = "q9") +
  theme(plot.title = element_text(size = 25, face="bold.italic"), axis.title.x = element_text(size=22), axis.title.y = element_text(size=22), axis.text.x = element_text(size=22), axis.text.y = element_text(size=22), legend.text=element_text(size=20))
p5 <- p5 + scale_color_gradientn(colors = c("grey", "white", "red"),limits=c(1, 5), breaks=c(1,3,5))
p6 <- FeaturePlot(sample_integrated, features = "SLC1A3", slot = "data", min.cutoff = "q9") +
  theme(plot.title = element_text(size = 25, face="bold.italic"), axis.title.x = element_text(size=22), axis.title.y = element_text(size=22), axis.text.x = element_text(size=22), axis.text.y = element_text(size=22), legend.text=element_text(size=20))
p6 <- p6 + scale_color_gradientn(colors = c("grey", "white", "red"),limits=c(0, 4), breaks=c(0,2,4))
p7 <- FeaturePlot(sample_integrated, features = "GLUL", slot = "data", min.cutoff = "q9") +
  theme(plot.title = element_text(size = 25, face="bold.italic"), axis.title.x = element_text(size=22), axis.title.y = element_text(size=22), axis.text.x = element_text(size=22), axis.text.y = element_text(size=22), legend.text=element_text(size=20))
p7 <- p7 + scale_color_gradientn(colors = c("grey", "white", "red"),limits=c(0, 4), breaks=c(0,2,4))
p8 <- FeaturePlot(sample_integrated, features = "NRXN1", slot = "data", min.cutoff = "q9") +
  theme(plot.title = element_text(size = 25, face="bold.italic"), axis.title.x = element_text(size=22), axis.title.y = element_text(size=22), axis.text.x = element_text(size=22), axis.text.y = element_text(size=22), legend.text=element_text(size=20))
p8 <- p8 + scale_color_gradientn(colors = c("grey", "white", "red"),limits=c(0, 4), breaks=c(0,2,4))
p9 <- FeaturePlot(sample_integrated, features = "CADM2", slot = "data", min.cutoff = "q9") +
  theme(plot.title = element_text(size = 25, face="bold.italic"), axis.title.x = element_text(size=22), axis.title.y = element_text(size=22), axis.text.x = element_text(size=22), axis.text.y = element_text(size=22), legend.text=element_text(size=20))
p9 <- p9 + scale_color_gradientn(colors = c("grey", "white", "red"),limits=c(0, 4), breaks=c(0,2,4))
p10 <- FeaturePlot(sample_integrated, features = "PTN", slot = "data", min.cutoff = "q9") +
  theme(plot.title = element_text(size = 25, face="bold.italic"), axis.title.x = element_text(size=22), axis.title.y = element_text(size=22), axis.text.x = element_text(size=22), axis.text.y = element_text(size=22), legend.text=element_text(size=20))
p10 <- p10 + scale_color_gradientn(colors = c("grey", "white", "red"),limits=c(0, 3), breaks=c(0,1,2,3))
p11 <- FeaturePlot(sample_integrated, features = "GPC5", slot = "data", min.cutoff = "q9") +
  theme(plot.title = element_text(size = 25, face="bold.italic"), axis.title.x = element_text(size=22), axis.title.y = element_text(size=22), axis.text.x = element_text(size=22), axis.text.y = element_text(size=22), legend.text=element_text(size=20))
p11 <- p11 + scale_color_gradientn(colors = c("grey", "white", "red"),limits=c(0, 3), breaks=c(0,2,4))


pi <- list(p1, p2, p3, p4,p5,p6,p7,p8,p9,p10,p11)
Reduce( `+`, pi ) + plot_layout(ncol = 2) + plot_annotation(
  title = 'Heatmap of DAA Marker Genes (Human EC)',
  theme = theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5))
)
# FeaturePlot(sample_integrated, features = marker_features , ncol = 2) + theme(plot.title = element_text(size = 25, face="bold.italic"), axis.title.x = element_text(size=22), axis.title.y = element_text(size=22), axis.text.x = element_text(size=22), axis.text.y = element_text(size=22), legend.text=element_text(size=20)) 
dev.off()
