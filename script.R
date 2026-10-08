# Single-cell RNA-seq of human bone marrow (PanglaoDB SRA779509 / SRS3805259)
# QC, clustering, marker genes and manual cell type annotation with Seurat v5.
# Run from the repository root. Input: data/SRA779509_SRS3805259.sparse.RData

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

stopifnot(packageVersion("Seurat") >= "5.0.0")
set.seed(42)

data_file <- "SRA779509_SRS3805259.sparse.RData"
fig_dir <- "results/figures"
tab_dir <- "results/tables"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, name, w = 8, h = 6) {
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}

# 1. Load data ------------------------------------------------------------------
if (!file.exists(data_file)) stop("Input file not found: ", data_file)
env <- new.env()
load(data_file, envir = env)
sm <- env$sm  # sparse count matrix, genes x cells

# Row names are "symbol_ENSEMBLID" (some with a "_PAR_Y" suffix): keep the symbol only
symbols <- sub("_ENSG[0-9.]+(_PAR_Y)?$", "", rownames(sm))
symbols <- gsub("_", "-", symbols)
rownames(sm) <- make.unique(symbols)

bone <- CreateSeuratObject(counts = sm, project = "sc_bone", min.cells = 3, min.features = 200)
message("Cells x genes in the raw object: ", ncol(bone), " x ", nrow(bone))

# 2. Quality control -----------------------------------------------------------
bone[["percent.mt"]]  <- PercentageFeatureSet(bone, pattern = "^MT-")
bone[["percent.rbp"]] <- PercentageFeatureSet(bone, pattern = "^RP[SL][0-9]")
bone[["percent.hb"]]  <- PercentageFeatureSet(bone, pattern = "^HB[ABDEGMQZ][0-9]*$")

qc_feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rbp")
save_plot(VlnPlot(bone, features = qc_feats, ncol = 4, pt.size = 0), "01_qc_violin_before", w = 12, h = 5)
save_plot(FeatureScatter(bone, "nCount_RNA", "percent.mt") +
            FeatureScatter(bone, "nCount_RNA", "nFeature_RNA"), "02_qc_scatter_before", w = 12, h = 5)

# Thresholds chosen from the plots above: high mt fractions occur almost only in low-UMI cells;
# the lower gene cutoff is kept low because erythroid cells have few detected genes;
# the upper gene cutoff sits at the edge of the dense cloud
qc <- list(min_features = 250, max_features = 2500, max_mt = 10)

md <- bone@meta.data
p_mt <- ggplot(md, aes(nCount_RNA, percent.mt, colour = percent.hb)) +
  geom_point(size = 0.4) + scale_colour_viridis_c(name = "% Hb") +
  geom_hline(yintercept = qc$max_mt, linetype = "dashed") + theme_classic()
p_ft <- ggplot(md, aes(nCount_RNA, nFeature_RNA, colour = percent.hb)) +
  geom_point(size = 0.4) + scale_colour_viridis_c(name = "% Hb") +
  geom_hline(yintercept = c(qc$min_features, qc$max_features), linetype = "dashed") + theme_classic()
save_plot(p_mt + p_ft + plot_layout(guides = "collect"), "02b_qc_scatter_thresholds", w = 14, h = 5)

n_before <- ncol(bone)
bone <- subset(bone, subset = nFeature_RNA > qc$min_features &
                 nFeature_RNA < qc$max_features & percent.mt < qc$max_mt)
message("Cells after QC: ", ncol(bone), " of ", n_before)
save_plot(VlnPlot(bone, features = qc_feats, ncol = 4, pt.size = 0), "03_qc_violin_after", w = 12, h = 5)

# 3. Normalisation and scaling ----------------------------------------------------
bone <- NormalizeData(bone, normalization.method = "LogNormalize", scale.factor = 10000)
bone <- CellCycleScoring(bone,
                         s.features = cc.genes.updated.2019$s.genes,
                         g2m.features = cc.genes.updated.2019$g2m.genes,
                         set.ident = TRUE)
bone <- FindVariableFeatures(bone, selection.method = "vst", nfeatures = 2000)
save_plot(LabelPoints(plot = VariableFeaturePlot(bone), points = head(VariableFeatures(bone), 10),
                      repel = TRUE), "04_variable_features")
bone <- ScaleData(bone, features = rownames(bone))

# 4. PCA, clustering, UMAP ---------------------------------------------------------
bone <- RunPCA(bone, features = VariableFeatures(bone), verbose = FALSE)
save_plot(ElbowPlot(bone, ndims = 30), "05_elbow_plot", w = 6, h = 4)

n_pcs <- 20
bone <- FindNeighbors(bone, dims = 1:n_pcs)
bone <- FindClusters(bone, resolution = 0.5)
bone <- RunUMAP(bone, dims = 1:n_pcs, seed.use = 42)
save_plot(DimPlot(bone, reduction = "umap", label = TRUE), "06_umap_clusters")

cc_plot <- bone@meta.data %>%
  count(seurat_clusters, Phase) %>%
  group_by(seurat_clusters) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x = seurat_clusters, y = percent, fill = Phase)) +
  geom_col() +
  ggtitle("Cell cycle phase per cluster")
save_plot(cc_plot, "07_cell_cycle_by_cluster", w = 7, h = 5)

# 5. Marker genes ----------------------------------------------------------------------
Idents(bone) <- "seurat_clusters"
bone_markers <- FindAllMarkers(bone, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(bone_markers, file.path(tab_dir, "markers_all_clusters.csv"), row.names = FALSE)

top10 <- bone_markers %>% group_by(cluster) %>% slice_max(n = 10, order_by = avg_log2FC)
write.csv(top10, file.path(tab_dir, "top10_markers_per_cluster.csv"), row.names = FALSE)
save_plot(DoHeatmap(bone, features = top10$gene) + NoLegend(), "08_heatmap_top10_markers", w = 20, h = 15)

# Same, ranked by the difference in detection rate inside vs outside the cluster
top10_specific <- bone_markers %>%
  filter(p_val_adj < 0.05, pct.1 >= 0.4) %>%
  mutate(pct_diff = pct.1 - pct.2) %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = pct_diff)
write.csv(top10_specific, file.path(tab_dir, "top10_specific_markers_per_cluster.csv"), row.names = FALSE)

# 6. Annotation support ----------------------------------------------------------------
# Pairwise contrasts between related clusters (positive log2FC = higher in the first cluster).
# All cells come from one sample, so p-values are inflated: interpret log2FC and pct.1/pct.2
contrasts <- list(
  c1_vs_c2 = c("1", "2"),
  c1_vs_c3 = c("1", "3"),
  c3_vs_c4 = c("3", "4")
)
for (nm in names(contrasts)) {
  res <- FindMarkers(bone, ident.1 = contrasts[[nm]][1], ident.2 = contrasts[[nm]][2],
                     min.pct = 0.1, logfc.threshold = 0.25)
  res$gene <- rownames(res)
  write.csv(res[order(-res$avg_log2FC), ], file.path(tab_dir, paste0("contrast_", nm, ".csv")),
            row.names = FALSE)
}

canonical_markers <- c("CD3E", "IL7R", "CD8A", "TRDC", "TRGC1", "TRGC2",
                       "NKG7", "GNLY", "KLRD1", "KLRF1",
                       "MS4A1", "CD79A", "IGHD", "TCL1A", "JCHAIN", "MZB1", "SDC1",
                       "CD14", "LYZ", "VCAN", "S100A8", "FCGR3A", "LST1", "MS4A7", "CDKN1C",
                       "C1QA", "C1QB", "CD163",
                       "ELANE", "MPO", "AZU1", "PRTN3",
                       "CLEC10A", "CD1C", "FCER1A", "LILRA4",
                       "HBM", "GYPA", "ALAS2", "HBB",
                       "CD34", "KIT", "AVP", "HLF", "CRHBP", "SPINK2")
save_plot(DotPlot(bone, features = canonical_markers) + RotatedAxis(),
          "10_dotplot_canonical_markers", w = 16, h = 6)

t_nk_markers <- c("CD3E", "CD3D", "TRAC", "CD4", "CD8A", "CD8B", "TRDC", "TRGC1", "TRGC2",
                  "CCR7", "SELL", "TCF7", "LEF1", "IL7R", "S100A4", "KLRB1", "GZMK", "GZMH",
                  "NCAM1", "KLRF1", "CD69")
t_nk <- subset(bone, idents = c("1", "2", "3", "4"))
save_plot(DotPlot(t_nk, features = t_nk_markers) + RotatedAxis(),
          "10b_dotplot_T_NK_subsets", w = 11, h = 4)

stress_genes <- c("FOS", "JUN", "JUNB", "DUSP1", "DUSP2", "RGS1", "HSPA1A", "HSPA1B",
                  "DNAJB1", "ZFP36L2")
save_plot(DotPlot(bone, features = stress_genes) + RotatedAxis(),
          "10c_dotplot_stress_genes", w = 9, h = 5)

# 7. Cell type annotation --------------------------------------------------------------
cluster_labels <- c(
  "0"  = "Classical monocytes",
  "1"  = "Memory T cells (IL7R+, effector memory-like)",
  "2"  = "Naive/central memory T cells",
  "3"  = "Cytotoxic effector CD8 T cells (TEMRA-like)",
  "4"  = "NK cells",
  "5"  = "Naive B cells",
  "6"  = "HSPC / myeloid-primed progenitors",
  "7"  = "Non-classical monocytes",
  "8"  = "Erythroid precursors",
  "9"  = "cDC2",
  "10" = "pDC"
)
stopifnot(setequal(levels(bone), names(cluster_labels)))

bone <- RenameIdents(bone, cluster_labels)
bone$cell_type <- Idents(bone)
save_plot(DimPlot(bone, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend(),
          "11_umap_annotated")

cell_type_counts <- as.data.frame(table(bone$cell_type))
colnames(cell_type_counts) <- c("cell_type", "n_cells")
write.csv(cell_type_counts, file.path(tab_dir, "cell_type_counts.csv"), row.names = FALSE)

# 8. Save ----------------------------------------------------------------------------------
saveRDS(bone, "results/bone_annotated.rds")
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
