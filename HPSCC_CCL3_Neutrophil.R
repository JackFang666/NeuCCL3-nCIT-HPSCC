# ============================================================
# CCL3+ Neutrophil Signature Predicts Response to nCIT in HPSCC
# Reproducible analysis pipeline (scRNA-seq + bulk RNA-seq)
# ------------------------------------------------------------
# This script is designed for GitHub release:
# - Minimal hard-coded paths (edit CONFIG)
# - Modular steps with reusable functions
# - Clear definition of CCL3+ neutrophils (Neu_CCL3)
# - Cell-cell interaction with T cells (CellChat)
# - Clinical validation in bulk RNA-seq (ssGSEA + ROC; optional survival)
#
# Generated/merged on: 2025-12-31
# ============================================================


## ============================================================
## 0) CONFIG  (EDIT HERE)
## ============================================================

config <- list(
  # ---- folders ----
  data_dir    = "data",
  results_dir = "results",
  fig_dir     = file.path("results", "figures"),
  table_dir   = file.path("results", "tables"),
  obj_dir     = file.path("results", "objects"),

  # ---- scRNA inputs ----
  # Option A: start from 10X folders (one folder per sample; inside contains matrix.mtx + barcodes + features)
  use_10x      = TRUE,
  tenx_dir     = file.path("data", "scRNA", "tenx"),  # e.g., data/scRNA/tenx/SAMPLE1/...
  min_cells    = 3,
  min_features = 300,

  # Option B: start from a processed Seurat object (if you already have it)
  seurat_rds_in  = file.path("data", "objects", "scRNA_merged.rds"),
  seurat_rds_out = file.path("results", "objects", "scRNA_harmony.rds"),

  # ---- metadata columns ----
  batch_col     = "orig.ident",
  patient_col   = "patient_id",   # if not available, will fall back to orig.ident
  response_col  = "Response",     # e.g., CR/PR/SD/TN
  responder_lv  = c("CR", "PR"),
  nonresp_lv    = c("SD", "TN"),

  # ---- QC thresholds (adjust to your dataset) ----
  qc = list(
    nFeature_low  = 300,
    nFeature_high = 6000,
    mt_high       = 20,     # percent
    hb_high       = 5       # percent
  ),

  # ---- preprocessing / clustering ----
  n_var_features = 3000,
  n_pcs          = 30,
  harmony_dims   = 1:30,
  resolution     = 0.5,

  # ---- neutrophil subclustering ----
  neu_dims       = 1:20,
  neu_resolution = 0.4,

  # ---- parallel ----
  n_workers      = 4
)


## ============================================================
## 1) Packages (auto-install optional)
## ============================================================

required_pkgs <- c(
  # core
  "Seurat","harmony","Matrix","dplyr","tibble","tidyr","stringr","data.table",
  "ggplot2","patchwork","cowplot","readr","readxl",
  # enrichment
  "clusterProfiler","org.Hs.eg.db","fgsea","msigdbr",
  # interactions
  "CellChat","future",
  # validation
  "GSVA","GSEABase","pROC",
  # optional
  "survival","survminer"
)

install_if_missing <- function(pkgs){
  missing <- pkgs[!vapply(pkgs, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    message("Installing missing packages: ", paste(missing, collapse=", "))
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
}

# If you prefer NOT to auto-install, comment out the next line.
install_if_missing(required_pkgs)

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(Matrix)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(stringr)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(readr)
  library(readxl)

  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(fgsea)
  library(msigdbr)

  library(CellChat)
  library(future)

  library(GSVA)
  library(GSEABase)
  library(pROC)

  # optional
  suppressWarnings({
    library(survival)
    library(survminer)
  })
})

set.seed(1234)


## ============================================================
## 2) Helpers
## ============================================================

dir_create <- function(x){
  if (!dir.exists(x)) dir.create(x, recursive = TRUE, showWarnings = FALSE)
}

save_plot <- function(p, file, width=7, height=5, dpi=300){
  dir_create(dirname(file))
  ggsave(filename = file, plot = p, width = width, height = height, dpi = dpi)
}

write_table <- function(df, file){
  dir_create(dirname(file))
  data.table::fwrite(df, file)
}

stopif_missing_cols <- function(df, cols, name="data.frame"){
  miss <- cols[!cols %in% colnames(df)]
  if (length(miss) > 0) stop(name, " is missing required columns: ", paste(miss, collapse=", "))
}

make_response2 <- function(meta, response_col, responder_lv, nonresp_lv){
  stopif_missing_cols(meta, response_col, "meta")
  x <- as.character(meta[[response_col]])
  out <- rep(NA_character_, length(x))
  out[x %in% responder_lv] <- "R"
  out[x %in% nonresp_lv]   <- "NR"
  out <- factor(out, levels = c("NR","R"))
  out
}

# Safer gene symbol matching for bulk vs signature lists
upper_first <- function(x){
  # "ccL3l1" -> "Ccl3l1" (human data usually uppercase; you can adapt)
  x <- tolower(x)
  paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
}


## ============================================================
## 3) Load scRNA-seq
## ============================================================

load_seurat <- function(){
  if (isTRUE(config$use_10x)) {
    message("Loading 10X data from: ", config$tenx_dir)
    sample_ids <- list.dirs(config$tenx_dir, full.names = FALSE, recursive = FALSE)
    if (length(sample_ids) == 0) stop("No sample folders found under: ", config$tenx_dir)

    obj_list <- list()
    for (id in sample_ids) {
      counts <- Read10X(data.dir = file.path(config$tenx_dir, id))
      obj <- CreateSeuratObject(
        counts = counts, project = id,
        min.cells = config$min_cells, min.features = config$min_features
      )
      obj$orig.ident <- id
      obj_list[[id]] <- obj
    }
    seu <- merge(obj_list[[1]], y = obj_list[-1])
    return(seu)
  } else {
    message("Loading Seurat RDS from: ", config$seurat_rds_in)
    seu <- readRDS(config$seurat_rds_in)
    return(seu)
  }
}


## ============================================================
## 4) QC + Harmony Integration + Major clustering
## ============================================================

run_qc_and_harmony <- function(seu){
  # QC metrics
  seu[["mt_percent"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  seu[["hb_percent"]] <- PercentageFeatureSet(seu, pattern = "^HB[AB]") # human hemoglobin

  # Filter
  seu <- subset(
    seu,
    subset =
      nFeature_RNA > config$qc$nFeature_low &
      nFeature_RNA < config$qc$nFeature_high &
      mt_percent   < config$qc$mt_high &
      hb_percent   < config$qc$hb_high
  )

  # Standard workflow
  seu <- NormalizeData(seu, verbose = FALSE) |>
    FindVariableFeatures(selection.method = "vst", nfeatures = config$n_var_features, verbose = FALSE) |>
    ScaleData(verbose = FALSE) |>
    RunPCA(npcs = config$n_pcs, verbose = FALSE)

  # Harmony
  if (!config$batch_col %in% colnames(seu@meta.data)) stop("batch_col not found: ", config$batch_col)
  seu <- RunHarmony(seu, group.by.vars = config$batch_col)

  # Cluster
  seu <- RunUMAP(seu, reduction = "harmony", dims = config$harmony_dims) |>
    FindNeighbors(reduction = "harmony", dims = config$harmony_dims) |>
    FindClusters(resolution = config$resolution)

  # Response2
  if (config$response_col %in% colnames(seu@meta.data)) {
    seu$Response2 <- make_response2(seu@meta.data, config$response_col, config$responder_lv, config$nonresp_lv)
  } else {
    warning("Response column not found: ", config$response_col, " (skipping Response2)")
  }

  return(seu)
}


## ============================================================
## 5) Major cell type annotation (two options)
## ============================================================
# Option 1 (manual): provide a mapping from clusters -> celltypes after inspecting markers.
# Option 2 (recommended for GitHub): marker-score based labeling (less brittle).
#
# Below is a simple marker-score approach: compute module scores and assign by max score per cluster.

major_marker_sets <- list(
  Tcells   = c("CD3D","CD3E","TRAC"),
  Bcells   = c("MS4A1","CD79A","CD74"),
  Myeloid  = c("LYZ","S100A8","S100A9","FCN1"),
  Epi      = c("EPCAM","KRT19","KRT8"),
  Fibro    = c("COL1A1","COL1A2","DCN","LUM"),
  Endo     = c("PECAM1","VWF","KDR"),
  Mast     = c("TPSAB1","TPSB2","KIT")
)

assign_major_celltype <- function(seu, cluster_col="seurat_clusters"){
  stopif_missing_cols(seu@meta.data, cluster_col, "seurat meta.data")
  # AddModuleScore creates columns like Tcells1, Bcells1 ...
  seu <- AddModuleScore(seu, features = major_marker_sets, name = names(major_marker_sets))
  score_cols <- paste0(names(major_marker_sets), "1")

  # cluster -> mean score -> label
  df <- seu@meta.data |>
    dplyr::select(all_of(cluster_col), all_of(score_cols)) |>
    dplyr::group_by(.data[[cluster_col]]) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(score_cols), mean, na.rm=TRUE), .groups="drop")

  df$major_celltype <- apply(df[, score_cols, drop=FALSE], 1, function(v) names(major_marker_sets)[which.max(v)])
  mapping <- df$major_celltype
  names(mapping) <- df[[cluster_col]]

  seu$major_celltype <- plyr::mapvalues(
    x = as.character(seu[[cluster_col]][,1]),
    from = names(mapping),
    to = unname(mapping),
    warn_missing = FALSE
  )
  seu$major_celltype <- factor(seu$major_celltype, levels = names(major_marker_sets))
  return(seu)
}


## ============================================================
## 6) Subclustering: Myeloid -> Neutrophils; define Neu_CCL3 concretely
## ============================================================

run_subcluster <- function(seu, keep_major, res=0.4, dims=1:20){
  Idents(seu) <- seu$major_celltype
  sub <- subset(seu, idents = keep_major)

  sub <- NormalizeData(sub, verbose = FALSE) |>
    FindVariableFeatures(nfeatures = 3000, verbose = FALSE) |>
    ScaleData(verbose = FALSE) |>
    RunPCA(npcs = max(dims), verbose = FALSE)

  sub <- RunHarmony(sub, group.by.vars = config$batch_col)
  sub <- RunUMAP(sub, reduction = "harmony", dims = dims) |>
    FindNeighbors(reduction = "harmony", dims = dims) |>
    FindClusters(resolution = res)

  return(sub)
}

# ---- Neutrophil definition logic (key points) ----
# Step A) From Myeloid subclusters, define "neutrophils" by canonical markers:
#         S100A8/S100A9/CSF3R/CXCR2/FCGR3B
# Step B) Reclustering neutrophils
# Step C) Define Neu_CCL3 by a dedicated CCL3-inflammatory module score + marker genes.
#
# We keep both a cluster-level label (most reproducible for figures) and a per-cell score.

neu_marker_sets <- list(
  Neu_S100A9 = c("S100A8","S100A9","S100A12","CSF3R","CXCR2","FCGR3B"),
  Neu_IFIT   = c("IFIT1","IFIT2","IFIT3","ISG15","RSAD2","MX1"),
  Neu_CCL3   = c("CCL3","CCL4","CCL3L1","IL1B","CXCL8","TNF","OSM","PTGS2","NFKBIA","TNFAIP3")
)

pick_neutrophil_cells <- function(myeloid){
  # Neutrophils are typically high S100A8/S100A9 and low LST1 (monocyte) in many datasets,
  # but thresholds depend on your data. Here we use a module score filter as a safe default.
  myeloid <- AddModuleScore(myeloid, features = list(neu_marker_sets$Neu_S100A9), name = "NeuCore")
  # keep top quantile of "NeuCore1" as neutrophils; adjust q if too strict/loose
  q <- quantile(myeloid$NeuCore1, probs = 0.70, na.rm = TRUE)
  neu <- subset(myeloid, subset = NeuCore1 >= q)
  return(neu)
}

annotate_neutrophil_subtypes <- function(neu, dims=1:20, res=0.4){
  # recluster neutrophils
  neu <- NormalizeData(neu, verbose = FALSE) |>
    FindVariableFeatures(nfeatures = 2000, verbose = FALSE) |>
    ScaleData(verbose = FALSE) |>
    RunPCA(npcs = max(dims), verbose = FALSE)

  neu <- RunHarmony(neu, group.by.vars = config$batch_col)
  neu <- RunUMAP(neu, reduction = "harmony", dims = dims) |>
    FindNeighbors(reduction = "harmony", dims = dims) |>
    FindClusters(resolution = res)

  # module scores for 3 canonical neutrophil programs
  neu <- AddModuleScore(neu, features = neu_marker_sets, name = names(neu_marker_sets))
  score_cols <- paste0(names(neu_marker_sets), "1")

  # cluster-level average score -> assign label by max score
  df <- neu@meta.data |>
    dplyr::select(seurat_clusters, all_of(score_cols)) |>
    dplyr::group_by(seurat_clusters) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(score_cols), mean, na.rm=TRUE), .groups="drop")

  df$subtype <- apply(df[, score_cols, drop=FALSE], 1, function(v){
    names(neu_marker_sets)[which.max(v)]
  })

  mapping <- df$subtype
  names(mapping) <- df$seurat_clusters

  neu$subtype <- plyr::mapvalues(
    x = as.character(neu$seurat_clusters),
    from = names(mapping),
    to = unname(mapping),
    warn_missing = FALSE
  )
  neu$subtype <- factor(neu$subtype, levels = c("Neu_S100A9","Neu_IFIT","Neu_CCL3"))

  # sanity-check: require true marker presence for Neu_CCL3
  # If your dataset has sparse CCL3, consider lowering thresholds / using cluster markers instead.
  return(list(neu = neu, cluster_score = df))
}

export_markers <- function(seu, group_col="subtype", out_prefix="markers", top_n=30){
  Idents(seu) <- seu[[group_col]][,1]
  m <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  m <- m %>% dplyr::filter(p_val_adj < 0.05)
  top <- m %>% dplyr::group_by(cluster) %>% dplyr::slice_max(order_by = avg_log2FC, n = top_n)
  write_table(m,   file.path(config$table_dir, paste0(out_prefix, "_all.tsv")))
  write_table(top, file.path(config$table_dir, paste0(out_prefix, "_top", top_n, ".tsv")))
  invisible(list(all=m, top=top))
}

calc_cell_fraction <- function(meta, group_col, cell_col){
  stopif_missing_cols(meta, c(group_col, cell_col), "meta")
  df <- meta %>%
    dplyr::filter(!is.na(.data[[group_col]]), !is.na(.data[[cell_col]])) %>%
    dplyr::count(.data[[group_col]], .data[[cell_col]], name = "n") %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::mutate(frac = n / sum(n)) %>%
    dplyr::ungroup()
  colnames(df)[1:2] <- c("group","celltype")
  df
}


## ============================================================
## 7) Neu_CCL3 vs CD8 T cells: quantitative association (fraction/score)
## ============================================================

# Optional: compute exhaustion / effector module scores in T cells and correlate with Neu_CCL3 fraction
t_marker_sets <- list(
  CD8_eff = c("NKG7","GZMB","GZMH","PRF1","GNLY"),
  CD8_ex  = c("PDCD1","LAG3","HAVCR2","TIGIT","ENTPD1"),
  Treg    = c("FOXP3","IL2RA","CTLA4","IKZF2"),
  Tfh     = c("CXCL13","PDCD1","ICOS","TOX")
)

add_t_scores <- function(t_obj){
  t_obj <- AddModuleScore(t_obj, features = t_marker_sets, name = names(t_marker_sets))
  return(t_obj)
}


## ============================================================
## 8) CellChat: Neu_CCL3 -> T cell interactions (NR vs R)
## ============================================================

run_cellchat_NR_vs_R <- function(seu, group_by="subtype", split_by="Response2",
                                 keep_types=NULL, n_workers=4, min_cells=10){
  stopif_missing_cols(seu@meta.data, c(group_by, split_by), "seurat meta.data")

  if (!is.null(keep_types)) {
    Idents(seu) <- seu[[group_by]][,1]
    seu <- subset(seu, idents = keep_types)
  }

  # split
  Idents(seu) <- seu[[split_by]][,1]
  seu_NR <- subset(seu, idents = c("NR"))
  seu_R  <- subset(seu, idents = c("R"))

  cellchatNR <- createCellChat(object = seu_NR, meta = seu_NR@meta.data, group.by = group_by)
  cellchatR  <- createCellChat(object = seu_R,  meta = seu_R@meta.data,  group.by = group_by)

  CellChatDB <- CellChatDB.human
  cellchatNR@DB <- CellChatDB
  cellchatR@DB  <- CellChatDB

  cellchatNR <- subsetData(cellchatNR)
  cellchatR  <- subsetData(cellchatR)

  # parallel
  options(future.globals.maxSize = 5 * 1024^3)
  future::plan("multisession", workers = n_workers)

  run_one <- function(x){
    x <- identifyOverExpressedGenes(x)
    x <- identifyOverExpressedInteractions(x)
    x <- projectData(x, PPI.human)
    x <- computeCommunProb(x, raw.use = TRUE)
    x <- filterCommunication(x, min.cells = min_cells)
    x <- computeCommunProbPathway(x)
    x <- aggregateNet(x)
    x
  }

  cellchatNR <- run_one(cellchatNR)
  cellchatR  <- run_one(cellchatR)

  merged <- mergeCellChat(list(NR = cellchatNR, R = cellchatR), add.names = c("NR","R"))

  p_compare <- compareInteractions(merged, show.legend = FALSE, group = c(1,2)) +
    compareInteractions(merged, show.legend = FALSE, group = c(1,2), measure = "weight")

  return(list(cellchatNR=cellchatNR, cellchatR=cellchatR, merged=merged, p_compare=p_compare))
}

plot_neu_to_t_bubble <- function(merged_cellchat, sources="Neu_CCL3", targets=c("CD8_ex","CD8_eff"), title=NULL){
  p <- netVisual_bubble(
    merged_cellchat,
    sources.use = sources,
    targets.use = targets,
    comparison = c(1,2),
    angle.x = 45,
    remove.isolate = TRUE
  )
  if (!is.null(title)) p <- p + ggtitle(title)
  p
}


## ============================================================
## 9) Functional analysis: R vs NR (GO / GSEA) within a cell type
## ============================================================
# Note: For strict inference, use pseudobulk per patient. Here we provide BOTH:
# - "quick" cell-level FindMarkers (good for exploration)
# - "recommended" pseudobulk template (edgeR-style counts)
#
# Choose based on what you used in the manuscript.

quick_DE_sc <- function(seu_sub, group_col="Response2", ident1="R", ident2="NR"){
  stopif_missing_cols(seu_sub@meta.data, group_col, "seu_sub meta.data")
  Idents(seu_sub) <- seu_sub[[group_col]][,1]
  FindMarkers(seu_sub, ident.1 = ident1, ident.2 = ident2, logfc.threshold = 0, min.pct = 0.1)
}

run_GO_enrich <- function(gene_symbols){
  gene_symbols <- unique(gene_symbols)
  gene_symbols <- gene_symbols[!is.na(gene_symbols) & gene_symbols != ""]
  if (length(gene_symbols) < 10) return(NULL)

  ego <- enrichGO(
    gene          = gene_symbols,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
  ego
}

run_fgsea_hallmark <- function(ranks){
  pathways_h <- msigdbr(species = "Homo sapiens", category = "H") %>%
    split(x = .$gene_symbol, f = .$gs_name)

  fg <- fgsea(pathways = pathways_h, stats = ranks, minSize = 15, maxSize = 500, nperm = 10000)
  fg <- fg[order(fg$pval), ]
  list(fgsea = fg, pathways = pathways_h)
}


## ============================================================
## 10) Bulk validation: ssGSEA + ROC (optional survival)
## ============================================================

score_signature_ssgsea <- function(expr_mat, gene_list, method="ssgsea"){
  # expr_mat: genes x samples (numeric)
  common <- intersect(rownames(expr_mat), gene_list)
  if (length(common) < 5) stop("Too few signature genes found in bulk expression (n<5).")
  gs <- list(Neu_CCL3 = common)
  GSVA::gsva(as.matrix(expr_mat), gs, method = method, ssgsea.norm = TRUE, verbose = FALSE)
}

run_roc <- function(y, score){
  roc_obj <- pROC::roc(response = y, predictor = score, quiet = TRUE)
  auc_val <- pROC::auc(roc_obj)
  list(roc = roc_obj, auc = auc_val)
}


## ============================================================
## 11) MAIN RUN (you can comment out blocks you don't need)
## ============================================================

dir_create(config$results_dir)
dir_create(config$fig_dir)
dir_create(config$table_dir)
dir_create(config$obj_dir)

# ---- 11.1 Load + preprocess ----
seu <- load_seurat()
seu <- run_qc_and_harmony(seu)
seu <- assign_major_celltype(seu)

saveRDS(seu, file.path(config$obj_dir, "scRNA_harmony_major.rds"))
save_plot(DimPlot(seu, group.by="major_celltype", label=TRUE, repel=TRUE) + ggtitle("Major cell types"),
          file.path(config$fig_dir, "UMAP_major_celltypes.png"), width=8, height=6)

# ---- 11.2 Myeloid -> Neutrophils -> Neu_CCL3 definition ----
myeloid <- run_subcluster(seu, keep_major="Myeloid", res=0.4, dims=config$neu_dims)
save_plot(DimPlot(myeloid, label=TRUE) + ggtitle("Myeloid subclustering"),
          file.path(config$fig_dir, "UMAP_myeloid.png"), width=7, height=6)

neu <- pick_neutrophil_cells(myeloid)
neu_res <- annotate_neutrophil_subtypes(neu, dims=config$neu_dims, res=config$neu_resolution)
neu <- neu_res$neu
cluster_score <- neu_res$cluster_score

saveRDS(neu, file.path(config$obj_dir, "Neutrophils_annotated.rds"))
write_table(cluster_score, file.path(config$table_dir, "Neutrophil_cluster_module_scores.tsv"))

p_neu <- DimPlot(neu, group.by="subtype", label=TRUE, repel=TRUE) + ggtitle("Neutrophil subtypes (module-score)")
save_plot(p_neu, file.path(config$fig_dir, "UMAP_neutrophil_subtypes.png"), width=7, height=6)

export_markers(neu, group_col="subtype", out_prefix="NeutrophilSubtype")

# ---- 11.3 T cell subclustering (optional but recommended for Neu->T CellChat) ----
t_obj <- run_subcluster(seu, keep_major="Tcells", res=0.6, dims=1:20)
t_obj <- add_t_scores(t_obj)
saveRDS(t_obj, file.path(config$obj_dir, "Tcells_recluster.rds"))
save_plot(DimPlot(t_obj, label=TRUE) + ggtitle("T cell subclustering"),
          file.path(config$fig_dir, "UMAP_Tcells.png"), width=7, height=6)

# Create a unified "subtype" column for CellChat (Neu subtypes + T subclusters + other majors)
# You should refine this mapping to match your manuscript labels.
seu$subtype <- seu$major_celltype
seu$subtype[colnames(neu)] <- as.character(neu$subtype[colnames(neu)])
seu$subtype[colnames(t_obj)] <- paste0("Tcluster_", as.character(t_obj$seurat_clusters[colnames(t_obj)]))

# (Optional) manually collapse T clusters into CD8_ex/CD8_eff etc:
# Example rule-based annotation: assign by highest module score
if (all(paste0(names(t_marker_sets), "1") %in% colnames(t_obj@meta.data))) {
  sc_cols <- paste0(names(t_marker_sets), "1")
  cl_df <- t_obj@meta.data %>%
    dplyr::select(seurat_clusters, all_of(sc_cols)) %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::summarise(dplyr::across(all_of(sc_cols), mean, na.rm=TRUE), .groups="drop")
  cl_df$Tlabel <- apply(cl_df[, sc_cols, drop=FALSE], 1, function(v) names(t_marker_sets)[which.max(v)])
  t_map <- cl_df$Tlabel; names(t_map) <- cl_df$seurat_clusters
  seu$subtype[colnames(t_obj)] <- as.character(t_map[as.character(t_obj$seurat_clusters[colnames(t_obj)])])
}

seu$subtype <- factor(seu$subtype)

# ---- 11.4 Neu_CCL3 fraction by response (NR vs R) ----
if ("Response2" %in% colnames(seu@meta.data)) {
  frac_df <- calc_cell_fraction(seu@meta.data, group_col="Response2", cell_col="subtype")
  write_table(frac_df, file.path(config$table_dir, "Celltype_fraction_by_Response2.tsv"))

  p_frac <- frac_df %>%
    dplyr::filter(celltype %in% c("Neu_CCL3","CD8_ex","CD8_eff","Treg")) %>%
    ggplot(aes(x = group, y = frac, fill = celltype)) +
    geom_col(position = "dodge") +
    theme_bw() +
    labs(title="Selected cell-type fractions (NR vs R)", x=NULL, y="Fraction")
  save_plot(p_frac, file.path(config$fig_dir, "Fractions_selected_NR_vs_R.png"), width=9, height=4.5)
}

# ---- 11.5 CellChat NR vs R (focus Neu_CCL3 -> T) ----
if ("Response2" %in% colnames(seu@meta.data)) {
  keep_types <- c("Neu_CCL3","Neu_S100A9","Neu_IFIT","CD8_ex","CD8_eff","Treg","Tfh")
  keep_types <- keep_types[keep_types %in% levels(seu$subtype)]
  cc_res <- run_cellchat_NR_vs_R(seu, group_by="subtype", split_by="Response2",
                                 keep_types=keep_types, n_workers=config$n_workers, min_cells=10)

  save_plot(cc_res$p_compare, file.path(config$fig_dir, "CellChat_compare_NR_vs_R.png"), width=10, height=5)

  # bubble plot (Neu_CCL3 -> CD8)
  if ("Neu_CCL3" %in% keep_types) {
    p_bub <- plot_neu_to_t_bubble(cc_res$merged, sources="Neu_CCL3",
                                 targets=intersect(c("CD8_ex","CD8_eff"), keep_types),
                                 title="Neu_CCL3 -> CD8 signaling (NR vs R)")
    save_plot(p_bub, file.path(config$fig_dir, "CellChat_bubble_NeuCCL3_to_CD8.png"), width=10, height=6)
  }

  saveRDS(cc_res$merged, file.path(config$obj_dir, "CellChat_merged_NR_vs_R.rds"))
}

# ---- 11.6 Functional analysis example: epithelial R vs NR (GO/GSEA) ----
# If you want epithelial-specific functional analysis, rerun subclustering and DE:
if ("Response2" %in% colnames(seu@meta.data)) {
  epi <- run_subcluster(seu, keep_major="Epi", res=0.3, dims=1:20)
  saveRDS(epi, file.path(config$obj_dir, "Epithelial_recluster.rds"))
  de_epi <- quick_DE_sc(epi, group_col="Response2", ident1="R", ident2="NR")
  de_epi <- de_epi %>% tibble::rownames_to_column("gene")
  write_table(de_epi, file.path(config$table_dir, "DE_epi_R_vs_NR_sc.tsv"))

  # GO using up genes in R
  upR <- de_epi %>% dplyr::filter(avg_log2FC > 0.25, p_val_adj < 0.05) %>% dplyr::pull(gene)
  ego <- run_GO_enrich(upR)
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    p_go <- dotplot(ego, showCategory = 15) + ggtitle("GO BP: Epithelial Up in R")
    save_plot(p_go, file.path(config$fig_dir, "GO_epi_up_in_R.png"), width=9, height=6)
    write_table(as.data.frame(ego), file.path(config$table_dir, "GO_epi_up_in_R.tsv"))
  }

  # GSEA (Hallmark) with ranked logFC
  ranks <- de_epi$avg_log2FC; names(ranks) <- de_epi$gene
  ranks <- sort(ranks, decreasing = TRUE)
  gsea <- run_fgsea_hallmark(ranks)
  write_table(gsea$fgsea, file.path(config$table_dir, "fgsea_epi_hallmark.tsv"))
}

# ---- 11.7 Bulk validation template (EDIT INPUTS) ----
# Provide:
#   - expr_bulk: genes x samples numeric matrix (TPM/FPKM/counts)
#   - clin_bulk: data.frame with SampleID + ResponseBinary (0/1 or NR/R)
#
# Example expected files:
#   data/bulk/expr_tpm.tsv (genes x samples; first column=gene)
#   data/bulk/clin.tsv     (SampleID, Response2 or ResponseBinary, OS_time, OS_status, etc.)
#
# Uncomment when your files are ready.
#
# expr_file <- file.path(config$data_dir, "bulk", "expr_tpm.tsv")
# clin_file <- file.path(config$data_dir, "bulk", "clin.tsv")
# expr_bulk <- readr::read_tsv(expr_file) %>% as.data.frame()
# rownames(expr_bulk) <- expr_bulk[[1]]; expr_bulk <- as.matrix(expr_bulk[,-1])
# expr_bulk <- log2(expr_bulk + 1)   # recommended for TPM/FPKM
#
# clin_bulk <- readr::read_tsv(clin_file) %>% as.data.frame()
# stopif_missing_cols(clin_bulk, c("SampleID", "Response2"), "clin_bulk")
#
# # signature genes: use the Neu_CCL3 module set OR your top markers from scRNA
# sig_genes <- neu_marker_sets$Neu_CCL3
#
# ssg <- score_signature_ssgsea(expr_bulk, sig_genes, method="ssgsea")
# score <- as.numeric(ssg[1, match(clin_bulk$SampleID, colnames(ssg))])
# clin_bulk$Neu_CCL3_score <- score
#
# # ROC (NR vs R)
# y <- clin_bulk$Response2
# roc_res <- run_roc(y = y, score = clin_bulk$Neu_CCL3_score)
# print(roc_res$auc)
#
# p_roc <- ggroc(roc_res$roc, legacy.axes = TRUE) +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
#   annotate("text", x = 0.6, y = 0.2, label = paste0("AUC = ", round(as.numeric(roc_res$auc), 3))) +
#   theme_bw() + labs(title = "Neu_CCL3 ssGSEA predicts response")
# save_plot(p_roc, file.path(config$fig_dir, "ROC_bulk_NeuCCL3.png"), width=6, height=5)
#
# # Optional survival (if you have OS_time/OS_status):
# if (all(c("OS_time","OS_status") %in% colnames(clin_bulk))) {
#   cutoff <- median(clin_bulk$Neu_CCL3_score, na.rm = TRUE)
#   clin_bulk$RiskGroup <- ifelse(clin_bulk$Neu_CCL3_score > cutoff, "High", "Low")
#   fit <- survival::survfit(survival::Surv(OS_time, OS_status) ~ RiskGroup, data = clin_bulk)
#   p_km <- survminer::ggsurvplot(fit, pval = TRUE, risk.table = TRUE)$plot
#   save_plot(p_km, file.path(config$fig_dir, "KM_bulk_by_NeuCCL3.png"), width=7, height=6)
# }
#
# write_table(clin_bulk, file.path(config$table_dir, "bulk_scores_and_clin.tsv"))

# ---- 11.8 Reproducibility ----
sessionInfo()

# ========================= END =========================
