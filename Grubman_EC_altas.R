#Description: this R script provides the single cell data analysis for maunscript: A single-cell altas of entorhinal cortex from individuals with Alzheimer's disease revelas cell-type-specific gene expression regulation
#Dataset: GEO -- GSE138852


# dir_work = ""  --provide a dir_work

setwd(dir_work)

dir_data = "../data"

library(Seurat)
library(dplyr)
library(Matrix)
library(BRETIGEA)
library(rray)
library(ggplot2)
library(patchwork)
library(grid)
library(cowplot)

dir_counts = paste(dir_data, "GSE138852_counts.csv", sep = '/')
raw_counts<-read.table(dir_counts,sep=',', header = T, row.names = 1)

ad_samples <- grep(pattern = "AD", x = colnames(raw_counts), value = TRUE)

ad_raw_counts = raw_counts[, ad_samples]

ct_samples <- grep(pattern = "Ct", x = colnames(raw_counts), value = TRUE)
ct_raw_counts = raw_counts[, ct_samples]

ad_obj <- CreateSeuratObject(counts = ad_raw_counts,project = "AD")
ad_obj$genotype <- "AD"

ct_obj <- CreateSeuratObject(counts = ct_raw_counts,project = "CT")
ct_obj$genotype <- "CONTROL"

data <- merge(x = ad_obj, y = ct_obj,  add.cell.ids = c("AD", "CT"), project = "JF_PROJ")

rm(raw_counts)
rm(ad_raw_counts)
rm(ad_obj)
rm(ct_raw_counts)
rm(ct_obj)

data <- NormalizeData(object = data, normalization.method = "LogNormalize",
                      scale.factor = 10000)

data <- FindVariableFeatures(data, selection.method = "mean.var.plot", mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
data <- ScaleData(data)


data <- RunPCA(data, features = VariableFeatures(object = data))
data <- FindNeighbors(data, dims = 1:19)
data <- FindClusters(data, resolution = 0.8)


data <- RunUMAP(object = data, dims = 1:19)
UMAPPlot(object = data, label = T)


# dir_plot = c('../plots/heatmap.png')
# png(dir_plot, width = 15, height = 6, units = 'in', res = 300)
# feature_genes = c("GRIN1", "SYT1", "RBFOX3", "SNAP25", "GRIN2A", "SLC17A7", "SATB2", "GAD1", "GAD2", "SST", "NPY", "PLP1", "MBP", "CLDN11", "MOG", "CAMK2A", "CCK", "ATP1A1", "PDE10A", "TAC1", "PENK", "SDK2", "COL5A1", "NPTX1", "PDE1A", "ETL4","SLC1A2","GJA1","AQP4","HEXB","CSF1R","C1QA","P2RY12","PDGFRA","VCAN","CSPG4","OLIG1","VTN","FLT1","CLDN5")
# DoHeatmap(data, features = feature_genes, size = 3, slot = "data")
# dev.off()


all_nuclei = colnames(data)

cluster3_nuclei = all_nuclei[which(data$seurat_clusters == 3)]
cluster5_nuclei = all_nuclei[which(data$seurat_clusters == 5)]

astrocyte_nuclei = c(cluster3_nuclei, cluster5_nuclei)
astrocyte_data = data[, astrocyte_nuclei]



#====subcluster astrocyte===#

astrocyte_data <- FindVariableFeatures(astrocyte_data, selection.method = "vst", nfeatures = 2000)

astrocyte_data <- ScaleData(astrocyte_data)

astrocyte_data <- RunPCA(astrocyte_data)

# astrocyte_data <- JackStraw(astrocyte_data, num.replicate = 50, dims = 40)
# 
# astrocyte_data <- ScoreJackStraw(astrocyte_data, dims = 1:40)
# 
# JackStrawPlot(astrocyte_data, dims = 1:40)


# 
astrocyte_data <- FindNeighbors(astrocyte_data, dims = 1:10)
astrocyte_data <- FindClusters(astrocyte_data, resolution = 0.1)
astrocyte_data <- RunTSNE(object = astrocyte_data, dims = 1:10)
TSNEPlot(object = astrocyte_data, label = TRUE)


dir_plot = c('../plots/astrocyte_subcluster.pdf')
pdf(dir_plot, width = 9, height = 9)
DimPlot(astrocyte_data, reduction = "tsne", label = FALSE, pt.size = 1.5) + scale_color_manual(values=c("#4682B4", "#3CB371", "#DC143C", "#FFA500", "#BA55D3")) + NoLegend() +
  theme(axis.title.x = element_text(size=22), axis.title.y = element_text(size=22), axis.text.x = element_text(size=22), axis.text.y = element_text(size=22), legend.text=element_text(size=22))
dev.off()

all.genes = rownames(astrocyte_data)

cluster_markers = FindMarkers(astrocyte_data, ident.1 = 1,  test.use = "MAST", features = all.genes, logfc.threshold = 0.0, min.pct = 0.000000001)
dir_micro_markers = c('../result/astrocyte_subcluster1_markers.txt')
write.table(cluster_markers, dir_micro_markers, sep = '\t', row.names = T, col.names = T)


#=====================stacked violin plot==========================#
# https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj,
                          feature,
                          pt.size = 0.1,
                          plot.margin = unit(c(-0.05, 0, -0.05, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... ) + scale_fill_manual(values=c("#4682B4", "#3CB371", "#DC143C", "#FFA500", "#BA55D3")) +
    xlab("") + ylab("") + ggtitle(feature) +
    theme(legend.position = "none",
          plot.title = element_text(size = 22),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = 20, angle = 0),
          axis.text.y = element_text(size = 20),
          plot.margin = plot.margin )
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0.2,
                          plot.margin = unit(c(-0.05, 0, -0.05, 0), "cm"),
                          ...) {

  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))

  # Add back x-axis title to bottom plot. patchwork is going to support this?
  # plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
  #   theme(axis.text.x=element_text(size = 20), axis.ticks.x = element_line())
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_blank(), axis.ticks.x = element_line())

  # change the y-axis tick to only max value
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x +
                            scale_y_continuous(limits = c(0, y)) +
                            expand_limits(y = y))

  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 2)

  return(p)
}

# features = c("GLUL", "NRXN1", "CADM2", "PTN", "GPC5")
features = c("GFAP", "CD44", "HSPB1", "TNC", "SLC1A2", "SLC1A3", "GLUL", "NRXN1", "CADM2", "PTN", "GPC5")
dir_plot = '../plots/stacked_vlnplot.png'
png(dir_plot, width = 12, height = 14, units = 'in', res = 300)
StackedVlnPlot(obj = astrocyte_data, features = features)
dev.off()
