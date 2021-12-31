## R version 4.0.3 ; Seurat version 2.3.2
rm(list = ls())
library(Seurat)
library(monocle)
library(ggplot2)
library(loomR)
library(ggplot2)
library(ggsignif)

# read the Loom file  
path <- "/share/home/sunruixia/AML/AML556/monocle/"
align <- connect(filename = paste0(path, "AML556_adata.loom"), 
    mode = "r +", skip.validate = TRUE)
attr.df <- align$get.attribute.df(MARGIN = 2, col.names = "CellID")
sample <- data.frame(V1 = attr.df$CellID, ct = attr.df[, c("ct")])
sample$ct[which(sample$ct == 1)] <- "ProMono/Mono"
sample$ct[which(sample$ct == 2)] <- "ProMono-like/Mono-like"
sample$ct[which(sample$ct == 3)] <- "Prog/GMP"
sample$ct <- factor(sample$ct, levels <- c("Prog/GMP", "ProMono/Mono",
    "ProMono-like/Mono-like"))
gene <- data.frame(gene_short_name = align$row.attrs$Gene[],
    stringsAsFactors = FALSE)
data <- t(align[["matrix"]][, ])
align$close_all()

# marker genes of C1, C2, C3
load(paste0(path, "AML556_RNA.Rdata"))
obj_sub <- subset(obj, subset = ct %in% c(1, 2, 3))
obj_sub <- NormalizeData(obj, normalization.method = "LogNormalize", 
    scale.factor = 10000)
markers <- FindAllMarkers(obj_sub, only.pos = TRUE, 
    min.pct = 0.25, logfc.threshold = 0.25)

# top100 markers of C1, C2, C3
for (i in c(1, 2, 3)) {
    k <- markers[which(markers$cluster == i), ]
    if (nrow(k) >= 100) {
        diff_gene <- rbind(diff_gene, k[1:100, ])
    } 
    else {
        diff_gene <- rbind(diff_gene, k)
    }
}

# trajectory analysis by monocle 2
for (i in 100) {
    diff_gene <- read.csv(paste0(path, "diff_gene_", i, ".csv"),
    header = TRUE, stringsAsFactors = FALSE)
    pd <- new("AnnotatedDataFrame", data = sample)
    fd <- new("AnnotatedDataFrame", data = gene)
    cds <- newCellDataSet(as(as.matrix(data), "sparseMatrix"), 
        phenoData = pd, featureData = fd, lowerDetectionLimit = 0.5, 
        expressionFamily = VGAM::negbinomial.size())
    cds <- estimateSizeFactors(cds)
    cds <- estimateDispersions(cds)
    cds <- detectGenes(cds, min_expr = 0.5)
#choose genes that define a cell"s progress
    fData(cds)$number <- rownames(fData(cds))
    sub <- fData(cds)[fData(cds)$gene_short_name %in% diff_gene$gene, ]
    ordering_genes <- sub$number
#Sets the features (e.g. genes) to be used for ordering cells in pseudotime.
    cds <- setOrderingFilter(cds, ordering_genes)
#reduce data dimensionality
    cds <- reduceDimension(cds, max_components = 2,
    reduction_method = "DDRTree", residualModelFormulaStr = NULL, 
    verbose = FALSE)
#order cells along the trajectory
    cds <- orderCells(cds, root_state = 2)
    n <- plot_cell_trajectory(cds, color_by = "ct", show_branch_points = FALSE)
    ggsave(n, file = paste0(path, "cds_", i, ".png"))
    n <- plot_cell_trajectory(cds, color_by = "Pseudotime",
    show_branch_points = FALSE)
    ggsave(n, file = paste0(path, "cds_", i, "_Pseudotime.png"))
}

# boxplot of Pseudotime analysis
df <- pData(cds)[, c(1, 2, 5)]
compaired <- list(c("ProMono/Mono", "Prog/GMP"),
    c("ProMono/Mono", "ProMono-like/Mono-like"),
    c("Prog/GMP", "ProMono-like/Mono-like"))
n <- ggplot(df, aes(x = ct, y = Pseudotime, fill = ct))  +
    geom_boxplot() +
    theme(panel.grid = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(), 
        legend.position = "none",
        axis.text.x = element_text(size = 15,
        angle = 30, hjust = 1, vjust = 1), 
        axis.title = element_text(size = 15)) +
    ylab("Pseudotime") +
    geom_signif(comparisons = compaired, 
        map_signif_level = T, textsize = 5,
        test = "wilcox.test", step_increase = 0.05, 
        tip_length = 0)
ggsave(n, file = paste0(path, "Pseudotime_score.png"), 
    width = 4, height = 7)

## Analyzing Branches in Single-Cell Trajectories by BEAM
BEAM_res <- BEAM(cds, branch_point = 1, cores = 10)
BEAM_res <- BEAM_res[order(BEAM_res$qval), ]
BEAM_res <- BEAM_res[, c("gene_short_name", "pval", "qval")]
pdf(file = paste0(path, "BEAM_res_qval_10_6.pdf"))
plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,
     qval < 1e-6)), ], branch_point = 1, num_clusters = 4,
    cores = 20, use_gene_short_name = T, 
    show_rownames = F, return_heatmap = TRUE)
dev.off()

# generate the gene order of BEAM
row_cluster <- cutree(data$ph$tree_row, k = 4)
newOrder <- data$heatmap_matrix[data$ph$tree_row$order, ]
newOrder <- data.frame(newOrder)
newOrder$Cluster <- row_cluster[match(rownames(newOrder),
    rownames(data$heatmap_matrix))]
write.csv(newOrder, paste0(path, "cluster_list_out_pval_10_6.csv"))
gene_order <- data.frame(gene = rownames(newOrder), 
    cluster = newOrder[, c(201)])
write.csv(gene_order, file = paste0(path, "cds_cluster_list_out_pval_10_6.csv"))

# plot genes enriched in immune-related GO terms
gene <- c("GRN", "NEAT1", "BTG1", "JUN", "LYST", "CFD")
sub <- BEAM_res[which(BEAM_res$gene_short_name %in% m), 
    c("gene_short_name", "qval")]
sub <- sub[which(sub$qval <= 1e-10), ]$gene_short_name

# plot_genes_in_pseudotime shows two kinetic trends, 
# one for each lineage, instead of one
for (i in 1:5) {
    for (j in (i + 1):6) {
        n <- c(sub[i], sub[j])
        print(n)
        pic <- c()
        cds_genes <- row.names(subset(fData(cds), gene_short_name %in% n))
        try(pic <- plot_genes_branched_pseudotime(cds[cds_genes, ],
                            branch_point = 1,
                            color_by = "ct",
                            ncol = 1,
                            branch_labels = c("AML", "Normal")))
        ggsave(pic, file = paste0(path, "BEAM_res_", sub[i], 
        "_", sub[j], ".pdf"))
    }
}

## cellphoneDB
load(paste0(path, "AML556_RNA.Rdata"))
obj <- subset(obj, subset = ct %in% c(1, 2, 0))
obj <- NormalizeData(obj, normalization.method = "LogNormalize", 
    scale.factor = 10000)
write.table(as.matrix(obj@assays$RNA@data), "cellphonedb_count_T.txt", 
    sep = "\t", quote = F)
meta_data <- cbind(cell = rownames(obj@meta.data), 
    obj@meta.data[, "ct", drop = F])
write.table(meta_data, paste0(path, "cellphonedb_meta_T.txt"), 
    sep = "\t", quote = F, row.names = F)

# cellphonedb
cellphonedb method statistical-analysis cellphonedb_meta_T.txt 
    cellphonedb_count_T.txt --counts-data=gene_name --project-name=out_T

# Plotting statistical method results
cellphonedb plot dot-plot
cellphonedb plot heatmap-plot cellphonedb_meta.txt

# plot
pvalues <- read.table(paste0(path, "pvalues.txt"),
    header = T, sep = "\t", stringsAsFactors = F)
pvalues <- pvalues[, c(2, 12:20)]
pvalues_melt <- melt(pvalues)
pvalues_melt <- pvalues_melt[which(pvalues_melt$value <= 0.05),]
pvalues_melt <- pvalues_melt[pvalues_melt$variable %in% c("X2.0","X0.2","X1.0","X0.1"),]
pvalues_melt$variable <- as.character(pvalues_melt$variable)
pvalues_melt$variable[pvalues_melt$variable == "X1.0"] <- "ProMono_Mono|T_CTL_NK"
pvalues_melt$variable[pvalues_melt$variable == "X2.0"] <- "ProMono-like_Mono-like|T_CTL_NK"
pvalues_melt$variable[pvalues_melt$variable == "X0.2"] <- "T_CTL_NK|ProMono-like_Mono-like"
pvalues_melt$variable[pvalues_melt$variable == "X0.1"] <- "T_CTL_NK|ProMono_Mono"

means <- read.table(paste0(path, "means.txt"),
    header = T, sep = "\t", stringsAsFactors = F)
means <- means[, c(2, 12:20)]
means_melt <- melt(means)
means_melt <- means_melt[means_melt$variable %in% c("X2.0","X0.2","X1.0","X0.1"),]

pvalues_melt$joinlab <- paste0(pvalues_melt$interacting_pair, "_", pvalues_melt$variable)
means_melt$joinlab <- paste0(means_melt$interacting_pair, "_", means_melt$variable)
means_melt <- means_melt[which(means_melt$joinlab %in% pvalues_melt$joinlab),]
df <- merge(pvalues_melt, means_melt, by = "joinlab")
n <- ggplot(df, aes(variable.x, interacting_pair.x)) + 
    geom_point(aes(color = exps, size = -log10(pvalues + 0.0001)) ) +
    scale_size_continuous(range = c(1.5, 4)) +
	scale_color_gradient2(high="red", mid  ="darkblue") +
	theme_bw() + 
    theme(axis.text.x = element_text(angle = 45,hjust = 0.92,vjust = 1)) +
	xlab("inter_cluster") +
    ylab("inter_pair")
ggsave(n, file=paste0(path, "cellphonedb_T.png"),
    width=5, height=7)
