suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(matrixStats)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(jsonlite)
  library(scTenifoldKnk)
})

# =========================================================
# 0. 全局参数
# =========================================================
BASE_ROOT <- Sys.getenv("MULTISTAD_DATA_ROOT", unset = "data")
DATASET <- "GSE183904"

SAME_SAMPLE_REFERENCE <- "GSM5573484_sample19"

MIN_CELLS_PER_SAMPLE <- 50
MIN_CELLS_REQUIRED <- 150

MAX_INPUT_CELLS <- 2000
MIN_GENE_PCT <- 0.05
MAX_GENES_FOR_KNK <- 1500

NC_NNET <- 10
TD_K <- 3
SEED <- 123

SKIP_IF_OUTPUT_EXISTS <- TRUE

# =========================================================
# 1. 路径
# =========================================================
SCRNA_ROOT <- file.path(BASE_ROOT, "scRNA")
PROCESSED_DIR <- file.path(BASE_ROOT, "processed_v2")
STEP13_DIR <- file.path(PROCESSED_DIR, "step13_scrna_celltype_anchor")
OUT_DIR <- file.path(Sys.getenv("MULTISTAD_OUTPUT_ROOT", unset = "results"), "virtual_knockout")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

META_FILE <- file.path(STEP13_DIR, DATASET, paste0(DATASET, "_coarse_celltype_metadata.csv"))

# =========================================================
# 2. 工具函数
# =========================================================
normalize_gene <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("\\.[0-9]+$", "", x)
  x <- toupper(x)
  x
}

safe_row_sums_detected <- function(mat) {
  if (inherits(mat, "sparseMatrix")) {
    Matrix::rowSums(mat > 0)
  } else {
    rowSums(mat > 0)
  }
}

merge_matrices_by_genes <- function(mat_list) {
  if (length(mat_list) == 1) return(mat_list[[1]])

  all_genes <- unique(unlist(lapply(mat_list, rownames)))
  merged_list <- vector("list", length(mat_list))

  for (i in seq_along(mat_list)) {
    mat <- mat_list[[i]]
    miss <- setdiff(all_genes, rownames(mat))
    if (length(miss) > 0) {
      zero_block <- matrix(
        0,
        nrow = length(miss),
        ncol = ncol(mat),
        dimnames = list(miss, colnames(mat))
      )
      mat <- rbind(mat, zero_block)
    }
    mat <- mat[all_genes, , drop = FALSE]
    merged_list[[i]] <- mat
  }

  do.call(cbind, merged_list)
}

load_metadata <- function(meta_file) {
  meta <- fread(meta_file, data.table = FALSE)
  if (ncol(meta) < 2) stop("Metadata file seems invalid: ", meta_file)

  colnames(meta)[1] <- "cell_key"
  meta$cell_key <- as.character(meta$cell_key)

  required_cols <- c("sample_name", "cell_id_raw", "predicted_celltype")
  miss <- setdiff(required_cols, colnames(meta))
  if (length(miss) > 0) {
    stop("Metadata missing required columns: ", paste(miss, collapse = ", "))
  }

  meta$sample_name <- as.character(meta$sample_name)
  meta$cell_id_raw <- as.character(meta$cell_id_raw)
  meta$predicted_celltype <- as.character(meta$predicted_celltype)
  meta
}

sample_file_exists <- function(sample_name) {
  f1 <- file.path(SCRNA_ROOT, "GSE183904", paste0(sample_name, ".csv.gz"))
  f2 <- file.path(SCRNA_ROOT, "GSE183904", paste0(sample_name, ".csv"))
  file.exists(f1) || file.exists(f2)
}

load_gse183904_sample_matrix <- function(sample_name) {
  f1 <- file.path(SCRNA_ROOT, "GSE183904", paste0(sample_name, ".csv.gz"))
  f2 <- file.path(SCRNA_ROOT, "GSE183904", paste0(sample_name, ".csv"))
  f <- if (file.exists(f1)) f1 else f2
  if (!file.exists(f)) stop("Cannot find GSE183904 matrix file for sample: ", sample_name)

  df <- tryCatch(
    fread(f, data.table = FALSE, check.names = FALSE),
    error = function(e) read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  )

  if (ncol(df) < 2) stop("Invalid GSE183904 matrix file: ", sample_name)

  gene_col <- colnames(df)[1]
  genes <- normalize_gene(df[[gene_col]])
  df[[gene_col]] <- NULL

  mat <- as.matrix(df)
  storage.mode(mat) <- "numeric"
  rownames(mat) <- genes

  keep <- !is.na(rownames(mat)) & rownames(mat) != ""
  mat <- mat[keep, , drop = FALSE]

  dup_n <- sum(duplicated(rownames(mat)))
  if (dup_n > 0) {
    cat("[INFO] GSE183904 duplicated genes merged in", sample_name, ":", dup_n, "\n")
    mat <- rowsum(mat, group = rownames(mat), reorder = FALSE)
  }

  mat
}

filter_and_select_genes <- function(counts_matrix, target_gene,
                                    min_gene_pct = 0.05,
                                    max_genes = 1500) {
  detected_n <- safe_row_sums_detected(counts_matrix)
  keep <- detected_n >= (ncol(counts_matrix) * min_gene_pct)
  mat <- counts_matrix[keep, , drop = FALSE]

  if (!(target_gene %in% rownames(mat)) && (target_gene %in% rownames(counts_matrix))) {
    mat <- rbind(mat, counts_matrix[target_gene, , drop = FALSE])
    mat <- mat[!duplicated(rownames(mat)), , drop = FALSE]
  }

  if (!(target_gene %in% rownames(mat))) {
    stop("Target gene not found after filtering: ", target_gene)
  }

  log_mat <- log1p(as.matrix(mat))
  vars <- matrixStats::rowVars(log_mat)
  names(vars) <- rownames(log_mat)
  vars <- sort(vars, decreasing = TRUE)

  top_genes <- names(vars)[seq_len(min(max_genes, length(vars)))]
  if (!(target_gene %in% top_genes)) {
    top_genes[length(top_genes)] <- target_gene
    top_genes <- unique(top_genes)
  }

  out <- as.matrix(mat[top_genes, , drop = FALSE])
  out
}

save_basic_plots <- function(dr_df, out_prefix) {
  if (nrow(dr_df) == 0) return(invisible(NULL))
  if (!("FC" %in% colnames(dr_df)) || !("p.adj" %in% colnames(dr_df))) return(invisible(NULL))

  dr_df$log2FC <- log2(dr_df$FC)
  dr_df$change <- "NOT"
  dr_df$change[dr_df$p.adj < 0.05 & dr_df$log2FC > 0.5] <- "UP"
  dr_df$change[dr_df$p.adj < 0.05 & dr_df$log2FC < -0.5] <- "DOWN"

  top_genes <- dr_df %>%
    filter(change != "NOT") %>%
    arrange(p.adj) %>%
    head(15) %>%
    pull(gene)

  pdf(paste0(out_prefix, "_volcano.pdf"), width = 8, height = 7)
  p1 <- ggplot(dr_df, aes(x = log2FC, y = -log10(p.adj))) +
    geom_point(aes(color = change), alpha = 0.8, size = 1.5) +
    scale_color_manual(values = c("UP" = "#DC0000FF", "DOWN" = "#4DBBD5FF", "NOT" = "grey85")) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "black", linewidth = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.5) +
    geom_text_repel(
      data = dr_df[dr_df$gene %in% top_genes, , drop = FALSE],
      aes(label = gene),
      size = 4,
      max.overlaps = 20
    ) +
    theme_classic() +
    labs(
      title = "scTenifoldKnk virtual KO volcano",
      x = "log2(Fold Change)",
      y = "-log10(FDR)"
    ) +
    theme(
      text = element_text(family = "serif"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      axis.title = element_text(face = "bold", size = 13),
      axis.text = element_text(size = 11, color = "black"),
      legend.position = "none"
    )
  print(p1)
  dev.off()
}

get_usable_samples <- function(meta, target_celltype, min_cells_per_sample = 50) {
  tab <- meta[meta$predicted_celltype == target_celltype, , drop = FALSE] %>%
    dplyr::count(sample_name, name = "n_target_cells") %>%
    dplyr::mutate(local_file_exists = vapply(sample_name, sample_file_exists, logical(1))) %>%
    dplyr::arrange(dplyr::desc(n_target_cells))

  usable <- tab %>%
    dplyr::filter(local_file_exists, n_target_cells >= min_cells_per_sample) %>%
    dplyr::pull(sample_name)

  list(table = tab, usable = usable)
}

run_one_task <- function(meta, mode, target_celltype, target_gene, target_samples = NULL, prefix) {
  cat("\n====================================================\n")
  cat("[TASK]", prefix, "\n")
  cat("mode =", mode, "| celltype =", target_celltype, "| gene =", target_gene, "\n")

  summary_file <- file.path(OUT_DIR, paste0(prefix, "_run_summary.json"))
  if (SKIP_IF_OUTPUT_EXISTS && file.exists(summary_file)) {
    cat("[SKIP] Existing output found:", summary_file, "\n")
    return(data.frame(
      prefix = prefix,
      mode = mode,
      target_celltype = target_celltype,
      target_gene = target_gene,
      status = "skipped_existing",
      message = summary_file,
      stringsAsFactors = FALSE
    ))
  }

  tryCatch({
    if (mode == "same_sample") {
      meta_sub <- meta[
        meta$predicted_celltype == target_celltype &
          meta$sample_name %in% target_samples,
        ,
        drop = FALSE
      ]
      samples_used <- unique(meta_sub$sample_name)
    } else if (mode == "pooled_auto") {
      res <- get_usable_samples(meta, target_celltype, MIN_CELLS_PER_SAMPLE)
      write.csv(
        res$table,
        file.path(OUT_DIR, paste0(prefix, "_sample_count_table.csv")),
        row.names = FALSE
      )

      samples_used <- res$usable
      meta_sub <- meta[
        meta$predicted_celltype == target_celltype &
          meta$sample_name %in% samples_used,
        ,
        drop = FALSE
      ]
    } else {
      stop("Unsupported mode: ", mode)
    }

    if (nrow(meta_sub) < MIN_CELLS_REQUIRED) {
      stop("Too few cells after metadata filtering: ", nrow(meta_sub))
    }

    cat("[INFO] Selected cells in metadata:", nrow(meta_sub), "\n")
    cat("[INFO] Samples used:", paste(samples_used, collapse = ", "), "\n")

    mat_list <- list()

    for (s in samples_used) {
      message("[INFO] Loading sample: ", s)

      raw_ids <- meta_sub$cell_id_raw[meta_sub$sample_name == s]
      raw_ids <- as.character(raw_ids)

      mat <- load_gse183904_sample_matrix(s)
      keep_cells <- intersect(colnames(mat), raw_ids)

      if (length(keep_cells) == 0) {
        cat("[WARN] No overlapping cells for sample:", s, "\n")
        next
      }

      mat <- mat[, keep_cells, drop = FALSE]
      mat_list[[s]] <- mat
    }

    if (length(mat_list) == 0) {
      stop("No cells matched between metadata and raw matrices.")
    }

    counts_matrix <- merge_matrices_by_genes(mat_list)
    cat("[INFO] Raw merged matrix:", nrow(counts_matrix), "genes x", ncol(counts_matrix), "cells\n")

    if (ncol(counts_matrix) < MIN_CELLS_REQUIRED) {
      stop("Too few matched cells in raw matrix after merge: ", ncol(counts_matrix))
    }

    set.seed(SEED)
    if (ncol(counts_matrix) > MAX_INPUT_CELLS) {
      sampled_cells <- sample(colnames(counts_matrix), MAX_INPUT_CELLS)
      counts_matrix <- counts_matrix[, sampled_cells, drop = FALSE]
      cat("[INFO] Downsampled to", ncol(counts_matrix), "cells\n")
    }

    target_gene <- normalize_gene(target_gene)
    filtered_counts <- filter_and_select_genes(
      counts_matrix = counts_matrix,
      target_gene = target_gene,
      min_gene_pct = MIN_GENE_PCT,
      max_genes = MAX_GENES_FOR_KNK
    )

    cat("[INFO] Final matrix for scTenifoldKnk:", nrow(filtered_counts), "genes x", ncol(filtered_counts), "cells\n")
    cat("[INFO] Target gene locked:", target_gene, "\n")

    write.csv(as.matrix(filtered_counts), file.path(OUT_DIR, paste0(prefix, "_input_matrix.csv")))

    set.seed(SEED)
    nc_nCells_use <- min(500, ncol(filtered_counts))
    if (nc_nCells_use < 100) {
      warning("Very few cells for scTenifoldKnk random subnetworks: ", nc_nCells_use)
    }

    ko_result <- scTenifoldKnk(
      countMatrix = as.matrix(filtered_counts),
      gKO = target_gene,
      nc_nNet = NC_NNET,
      nc_nCells = nc_nCells_use,
      td_K = TD_K
    )

    if (!("diffRegulation" %in% names(ko_result))) {
      stop("scTenifoldKnk result does not contain diffRegulation.")
    }

    diff_regulation <- as.data.frame(ko_result$diffRegulation)
    if (!("gene" %in% colnames(diff_regulation))) {
      diff_regulation$gene <- rownames(diff_regulation)
    }

    write.csv(diff_regulation, file.path(OUT_DIR, paste0(prefix, "_diffRegulation.csv")), row.names = FALSE)

    sig_changes <- subset(diff_regulation, p.adj < 0.05)
    write.csv(sig_changes, file.path(OUT_DIR, paste0(prefix, "_significant_genes.csv")), row.names = FALSE)

    save_basic_plots(diff_regulation, file.path(OUT_DIR, prefix))

    run_summary <- list(
      dataset = DATASET,
      mode = mode,
      target_celltype = target_celltype,
      target_gene = target_gene,
      target_samples = samples_used,
      n_samples_used = length(samples_used),
      samples_used = samples_used,
      n_cells_final = ncol(filtered_counts),
      n_genes_final = nrow(filtered_counts),
      nc_nNet = NC_NNET,
      nc_nCells = nc_nCells_use,
      td_K = TD_K,
      n_sig_genes = nrow(sig_changes)
    )

    write_json(run_summary, summary_file, pretty = TRUE, auto_unbox = TRUE)

    cat("[DONE]", prefix, "\n")

    data.frame(
      prefix = prefix,
      mode = mode,
      target_celltype = target_celltype,
      target_gene = target_gene,
      status = "success",
      message = paste0("n_cells=", ncol(filtered_counts), "; n_sig=", nrow(sig_changes)),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    cat("[ERROR]", prefix, ":", conditionMessage(e), "\n")
    data.frame(
      prefix = prefix,
      mode = mode,
      target_celltype = target_celltype,
      target_gene = target_gene,
      status = "failed",
      message = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })
}

# =========================================================
# 3. 读取 metadata
# =========================================================
meta <- load_metadata(META_FILE)

# =========================================================
# 4. 定义任务表
# =========================================================
task_plan <- list(
  list(
    mode = "same_sample",
    target_celltype = "Fibroblast",
    target_gene = "PKM",
    target_samples = c(SAME_SAMPLE_REFERENCE),
    prefix = "GSE183904__Fibroblast__PKM"
  ),
  list(
    mode = "same_sample",
    target_celltype = "Fibroblast",
    target_gene = "RPL36AL",
    target_samples = c(SAME_SAMPLE_REFERENCE),
    prefix = "GSE183904__Fibroblast__RPL36AL"
  ),
  list(
    mode = "pooled_auto",
    target_celltype = "Fibroblast",
    target_gene = "PKM",
    target_samples = NULL,
    prefix = "GSE183904__Fibroblast__PKM__pooledAuto"
  ),
  list(
    mode = "pooled_auto",
    target_celltype = "Fibroblast",
    target_gene = "RPL36AL",
    target_samples = NULL,
    prefix = "GSE183904__Fibroblast__RPL36AL__pooledAuto"
  ),
  list(
    mode = "pooled_auto",
    target_celltype = "Epithelial",
    target_gene = "PKM",
    target_samples = NULL,
    prefix = "GSE183904__Epithelial__PKM__pooledAuto"
  ),
  list(
    mode = "pooled_auto",
    target_celltype = "Epithelial",
    target_gene = "RPL36AL",
    target_samples = NULL,
    prefix = "GSE183904__Epithelial__RPL36AL__pooledAuto"
  ),
  list(
    mode = "pooled_auto",
    target_celltype = "Myeloid",
    target_gene = "PKM",
    target_samples = NULL,
    prefix = "GSE183904__Myeloid__PKM__pooledAuto"
  ),
  list(
    mode = "pooled_auto",
    target_celltype = "Myeloid",
    target_gene = "RPL36AL",
    target_samples = NULL,
    prefix = "GSE183904__Myeloid__RPL36AL__pooledAuto"
  )
)

task_table <- do.call(rbind, lapply(task_plan, function(x) {
  data.frame(
    prefix = x$prefix,
    mode = x$mode,
    target_celltype = x$target_celltype,
    target_gene = x$target_gene,
    target_samples = if (is.null(x$target_samples)) NA else paste(x$target_samples, collapse = ";"),
    stringsAsFactors = FALSE
  )
}))
write.csv(task_table, file.path(OUT_DIR, "task_plan.csv"), row.names = FALSE)

# =========================================================
# 5. 顺序执行全部任务
# =========================================================
status_list <- list()

for (i in seq_along(task_plan)) {
  task <- task_plan[[i]]
  res <- run_one_task(
    meta = meta,
    mode = task$mode,
    target_celltype = task$target_celltype,
    target_gene = task$target_gene,
    target_samples = task$target_samples,
    prefix = task$prefix
  )
  status_list[[i]] <- res
  status_df <- do.call(rbind, status_list)
  write.csv(status_df, file.path(OUT_DIR, "task_status.csv"), row.names = FALSE)
}

final_status <- do.call(rbind, status_list)
write.csv(final_status, file.path(OUT_DIR, "task_status.csv"), row.names = FALSE)

cat("\n====================================================\n")
cat("[ALL DONE] Master pipeline finished.\n")
print(final_status)
cat("Outputs saved in:\n", OUT_DIR, "\n")
