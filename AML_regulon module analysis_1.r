## R version 4.0.3
## Correlation analysis of transcription factors (TFs) based on
## regulon activity scores (RAS) of total AML datasets
rm(list = ls())
library(ggplot2)
library(ggsignif)

#### import RAS data of total AML datasets
path <- "/share/home/sunruixia/AML/total/"
data <- read.table(paste0(path,  "scenic_matrix2.txt"), # 35389 196
    header = T,  sep = "\t",  row.names = 1)
meta <- read.table(paste0(path,  "scenic_meta.txt"), # 35389 7
    header = T,  sep = "\t",  row.names = 1)
AML <- c("HSC-like", "Prog-like", "GMP-like",
    "ProMono-like", "Mono-like", "cDC-like")
BM <- c("HSC", "Prog", "GMP", "ProMono", "Mono", "cDC")

#### correlation analysis
cor_matrix_data <- cor(data,  method = c("pearson"))

# unsupervised analysis by pheatmap
library(pheatmap)
p <- pheatmap(cor_matrix_data, 
    cluster_rows = TRUE, 
    cluster_cols = TRUE, 
    show_colnames = FALSE, 
    show_rownames = FALSE, 
    cutree_rows = 6, 
    cutree_cols = 6, 
    treeheight_col = 20, 
    treeheight_row = 20, 
    border_color = NA, 
)
ggsave(p,  file = paste0(path,  "cor_data_6_modules.pdf"))

# TFs order of modules generated from heatmap
data_order <- colnames(cor_matrix_data[p$tree_row$order,  p$tree_col$order])
data_cluster <- data.frame(module = cutree(p$tree_row, k = 6))
data_cluster$TF <- rownames(data_cluster)
data_cluster <- data_cluster[data_order, ]
write.csv(data_cluster,  file = paste0(path,  "cor_data_6_modules.csv"))

# Comparison of RAS between AML and normal cells of 6 modules
module <- data.frame()
for (i in 1:6) {
    cluster <- data_cluster[which(data_cluster$module == i), c("TF")]
    df_cluster <- data[,  which(colnames(data) %in% cluster)]
    module1 <- data.frame(mean_value = apply(df_cluster,  1,  mean), 
        Module = paste0("module",  i))
    module1$ct <- meta$ct
    module <- rbind(module,  module1)
}

my_comparisons <- list(c("AML", "Normal"))
num <- unique(module$Module)
celltype <- c(AML, BM)

for (n in num) {
    M6_mean <- module[which(module$Module == n & module$ct %in% celltype), ]
    M6_mean[which(M6_mean$ct %in% AML), 2] <- "AML"
    M6_mean[which(M6_mean$ct %in% BM), 2] <- "Normal"
    M6_mean$Module <- factor(M6_mean$Module, 
        levels = c("AML", "Normal"))
# boxplot of RAS between AML cells and normal cells
    p1 <- ggplot(M6_mean,  aes(x = Module, y = mean_value, fill = Module)) + 
        geom_boxplot(outlier.shape = NA) +
        geom_signif(Module = my_comparisons,
            map_signif_level = T, 
            textsize = 2, 
            test = "wilcox.test", 
            step_increase = 0, 
            tip_length = 0) +
        theme_bw() +
        theme(panel.grid = element_blank(), 
            panel.background = element_rect(fill = "transparent"), 
            axis.line = element_line(colour = "black"), 
            axis.text.x = element_text(size = 12, 
            angle = 30, hjust = 1, vjust = 1)) +
        ylab("value") + 
        xlab(n)  # set the X-axis and Y-axis titles
    ggsave(p1, file = paste0(path, n, ".pdf"),
        height = 3, width = 6)
}

#### Transcriptional regulatory networks built by regulons
## in AML cells and corresponding normal cells
regulon_data <- read.table(paste0(path, "TFtarget.tab"), header = T, sep = "\t")
regulon_data$log10_im <- log10(regulon_data$importance) # 3179662 4

# density plot
p2 <- ggplot(regulon_data, aes(log10_im)) + geom_density()

# filter markers of AML and corresponding normal clusters
marker <- read.table(paste0(path, "scenic_ct_markers.txt"),
    header = T, sep = "\t")
marker_filter <- marker[which(marker$avg_logFC >= log(1.5)
    & marker$p_val_adj <= 0.05 & marker$cluster %in% celltype), ]

#　top10 markers of AML and corresponding normal clusters
for (i in celltype) {
    k <- marker_filter[which(marker_filter$cluster == i), ]
    if (nrow(k) >= 10) {
        total_top10 <- rbind(total_top10, k[1:10, ])
    } 
    else {
        total_top10 <- rbind(total_top10, k)
    }
}

# choose the threhold for filtering the TFs
df <- data.frame()
for (i in seq(0, 2, 0.01)) {
    data_log10 <- regulon_data[which(regulon_data$log10_im >= i),]
    uni_TF <- data.frame(unique_TF = length(unique(data_log10$TF)))
    uni_target <- data.frame(unique_target = length(unique(data_log10$target)))
    top10 <- data.frame(total_top10 = length(intersect(total_top10$gene,
        data_log10$TF)))
    total <- cbind(uni_TF, uni_target, top10)
    total$threshold  <-  i
    df <- rbind(df, total)
}

# density plot
library(reshape2)
df_melt <- melt(df, id.var = "threshold")
p3 <- ggplot(df_melt, aes(x = threshold, y = value, color = variable)) +
    geom_line(aes(color = variable))
ggsave(p3, file = paste0(path, "unique_TF_target.png"))

# choose 1.0 as the threhold, and obtained 49 TFs 
# with 1195 target genes for Cytoscape
data_log10 <- data[which(data$log10_im >= 1),]
TF <- intersect(total_top10$gene, data_log10$TF) #49
data_select <- data_log10[which(data_log10$TF %in% TF),]
write.table(data_select, file = paste0(path, "data_select.tab"),
    row.names = FALSE)

