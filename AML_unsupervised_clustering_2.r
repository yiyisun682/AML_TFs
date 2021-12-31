## R version 4.0.3 ; Seurat version 2.3.2
rm(list = ls())
library(Seurat)
library(ggplot2)
library(ggsignif)
path <- "/share/home/sunruixia/AML/707B/"

## load data
load(paste0(path, "AML707B.Rdata"))

# Normalizing the data
adata <- NormalizeData(adata, normalization.method = "LogNormalize",
    scale.factor = 10000) # 9115 297

# Identification of highly variable features (feature selection)
adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 300)
plot1 <- VariableFeaturePlot(adata)

# scale
all.genes <- rownames(adata)
adata <- ScaleData(adata, features = all.genes)
adata <- RunPCA(adata, features = VariableFeatures(object = adata))

# Determine the ‘dimensionality’ of the dataset
for (i in 12) {
    adata <- FindNeighbors(adata, dims = seq(1, i))
    for (j in 0.3) {
        adata <- FindClusters(adata, resolution = j)
# Run non-linear dimensional reduction (UMAP/tSNE)
        adata <- RunTSNE(adata, dims = seq(1, i))
# Two-dimensional tSNE display 
        n <- DimPlot(adata, reduction = "tsne", label = T)
        ggsave(n, file = paste0(path, "seurat_707B_pc", i, "_re_",
            j, "_tsne_cluster.pdf"), width = 5.2, height = 4)
    }
}
# ERF two-dimensional tSNE display
n <- FeaturePlot(adata, features = c("ERF"))
ggsave(n, file = paste0(path, "ERF.pdf"), height = 3.8, width = 4.1)

# find highly activated TFs
Idents(adata) <- "RNA_snn_res.0.3"
markers <- FindMarkers(adata, ident.1 = "2", ident.2 = "3",
    only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# annotation
colnames(adata@meta.data)[8] <- "ct"
adata@meta.data[which(adata@meta.data$ct == 0), "ct"] <- "C0-T/CTL/NK"
adata@meta.data[which(adata@meta.data$ct == 1), "ct"] <- "C1-ProMono/Mono"
adata@meta.data[which(adata@meta.data$ct == 2), "ct"] <- "C2-ProMono-like/Mono-like"
adata@meta.data[which(adata@meta.data$ct == 3), "ct"] <- "C3-Prog/GMP"
adata@meta.data[which(adata@meta.data$ct == 4), "ct"] <- "C4-earlyEry/lateEry"
adata@meta.data[which(adata@meta.data$ct == 5), "ct"] <- "C5-B"
adata@meta.data[which(adata@meta.data$ct == 6), "ct"] <- "C6-ProB"
adata@meta.data[which(adata@meta.data$ct == 7), "ct"] <- "C7-Plasma"

# tSNE plot of AML707B
n <- DimPlot(adata, group.by = "ct", reduction = "tsne", 
    label = TRUE, cols = c("#32CD32", "#EEEE00", "#FFFF32",
    "#EE7FEE", "#00FFFF", "#1E90FF", "#9BCF33"))
ggsave(n, file = paste0("seurat_707B_pc", i, "_re_0.3_label.pdf"),
    width = 4.2, height = 3.8)

# TFs used for enrichment analysis
up_reg <- markers[which(markers$p_val_adj <= 0.05 
    & markers$avg_log2FC >= log(1.5)), ]
down_reg <- markers[which(markers$p_val_adj <= 0.05 
    & markers$avg_log2FC <= -log(1.5)), ]

# pheatmap of C2, C3
library(pheatmap)
df_707B <- adata@assays$RNA@data
anno_707B <- adata@meta.data[, "ct"]
df_707B_sub <- df_707B[which(df_707B$ct) %in% c("C2", "C3")),
    which(rownames(df_707B) %in% c(up_reg,down_reg))]
anno_707B_sub <- anno_707B[anno_707B$ct %in% c("C2", "C3")]
n <- pheatmap(df_707B_sub,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_colnames = FALSE,
    show_rownames = TRUE,
    gaps_col = c(1959),
    gaps_row = c(10),
    scale = "row",
    annotation_col = anno_707B_sub,
    color = colorRampPalette(
    c("navyblue", "white", "firebrick3"))(100))
ggsave(n, file = pste0(path, "AML_707B_hmap.pdf"), height = 4, width = 5)

# gene sets scoring
load(paste0(path, "AML707B_RNA.Rdata"))
# cell cycle scoring
obj <- NormalizeData(obj, normalization.method = "LogNormalize", 
    scale.factor = 10000)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, 
    set.ident = FALSE)
data_select$Phase <- factor(data_select$Phase, levels = c("G1", "S", "G2M"))

#cell cycle barplot
n <- ggplot(data_select, aes(x = ct, fill = Phase)) +
    geom_bar(stat = "count", position = "fill") +
    scale_fill_manual(values = c("#619CFF", "#00BA38", "#F8776D")) +
    theme(panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank()) +
    xlab("phase")
ggsave(n, file = paste0(path, "phase.pdf"), width = 3, height = 4)

# AmiGO 2
cell_diff <- read.csv(paste0(path, "gene_set_cell_diff.csv"))
cycle_arrest <- read.csv(paste0(path, "gene_set_cell cycle_arrest.csv"))
obj <- AddModuleScore(object = obj, features = list(cell_diff$gene), 
    name = "cell_diff")
obj <- AddModuleScore(object = obj, features = list(cycle_arrest$gene),
    name = "cycle_arrest")
data_select <- obj@meta.data[which(obj@meta.data$ct %in% c("2", "3")),
    c("ct", "cell_diff", "cycle_arrest")]

# barplot
compaired <- list(c("2", "3"))
for (i in c(2, 3)) {
    name <- colnames(data_select)[i]
    p <- ggplot(data_select, aes(x = ct, y = get(name), 
            fill = ct)) +
        geom_boxplot() +
        theme(panel.grid = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.background = element_blank(), 
            legend.position = "none") +
        ylab("score") +
        xlab(colnames(data)[i]) +
        geom_signif(comparisons = compaired, 
            map_signif_level = T, 
            textsize = 3,
            test = "wilcox.test", 
            step_increase = 0.05, 
            tip_length = 0)
    ggsave(p, file = paste0(path, colnames(data_select)[i], ".pdf"),
        width = 4, height = 4)
}
