# Deconstructs a Seurat object into objects that can be read into anndata
# See: https://jr-leary7.github.io/quarto-site/tutorials/Seurat_AnnData_Conversion.html

# ----------------------------------------------------------------------
# Setup
# ----------------------------------------------------------------------

wd <- "C:/Users/aabassetc2/Documents/github/benchmark-integration"

suppressPackageStartupMessages(library('qs2'))
suppressPackageStartupMessages(library('dplyr'))
library('readr')
library('magrittr')
library('Matrix')
suppressPackageStartupMessages(library('Seurat'))
suppressPackageStartupMessages(library('biomaRt'))
library('optparse')
suppressPackageStartupMessages(library('logr'))

import::here(digest, 'sha1')

win_path <- function(...) {
  normalizePath(file.path(...), winslash = "/", mustWork = FALSE)
}

import::from(
  .from = win_path(wd, 'R', 'tools', 'text_tools.R'),
  'random_hash', .character_only = TRUE
)
import::from(
  .from = win_path(wd, 'R', 'tools', 'file_io.R'),
  'load_rdata', .character_only = TRUE
)

load_rdata <- function(filepath) {
  load(filepath)
  return(get(ls()[ls() != "filepath"]))
}

random_hash <- function(digits = 6) {
  return(substr(sha1(runif(1, 1, 2^31 - 1), digits = 14), 1, digits))
}

seurat_version <- substr(packageVersion("Seurat"), 1, 3)

# ----------------------------------------------------------------------
# Reduction map
# ----------------------------------------------------------------------

REDUCTION_MAP <- list(
  "cca"          = c("cca", "umap.cca"),
  "rpca"         = c("rpca", "umap.rpca"),
  "harmony"      = c("harmony", "umap.harmony"),
  "scvi"         = c("integrated.scvi", "umap"),
  "unintegrated" = c("pca", "umap")
)

detect_integration_method <- function(seurat_obj) {
  available <- names(seurat_obj@reductions)
  for (method_name in names(REDUCTION_MAP)) {
    expected <- REDUCTION_MAP[[method_name]]
    if (all(expected %in% available)) {
      if (method_name == "unintegrated") {
        integration_reductions <- c("cca", "integrated.cca", "rpca", "integrated.rpca",
                                    "harmony", "integrated.scvi", "scvi")
        if (any(integration_reductions %in% available)) next
      }
      return(method_name)
    }
  }
  if (any(c("integrated.cca", "cca") %in% available)) return("cca")
  if (any(c("integrated.rpca", "rpca") %in% available)) return("rpca")
  if ("harmony" %in% available) return("harmony")
  if (any(c("integrated.scvi", "scvi") %in% available)) return("scvi")
  return("unintegrated")
}

# ----------------------------------------------------------------------
# Args
# ----------------------------------------------------------------------

option_list <- list(
  make_option(c("-i", "--input-dir"), default = 'data/private',
              metavar = 'data/private', type = "character",
              help = "input h5ad files"),
  make_option(c("-o", "--output-dir"), default = 'data/output',
              metavar = 'data/output', type = "character",
              help = "set the output directory"),
  make_option(c("-t", "--troubleshooting"), default = FALSE, action = "store_true",
              metavar = "FALSE", type = "logical",
              help = "enable if troubleshooting to prevent overwriting your files")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
troubleshooting <- opt[['troubleshooting']]

# ----------------------------------------------------------------------
# Start Log
# ----------------------------------------------------------------------

start_time <- Sys.time()
log <- log_open(paste0(
  "convert_seurat_to_h5ad-", random_hash(), '-',
  strftime(start_time, format = "%Y%m%d_%H%M%S"), '.log'
))
log_print(paste('Script started at:', start_time))

# ----------------------------------------------------------------------
# Read & Process
# ----------------------------------------------------------------------

log_print(paste(Sys.time(), 'Seurat', seurat_version))

input_dir <- win_path(wd, opt[['input-dir']])
filenames <- list.files(input_dir,
                        pattern = "\\.(rds|rdata|qs)$",
                        recursive = FALSE, full.names = FALSE)
log_print(paste(Sys.time(), 'Files found...', paste(filenames, collapse = ', ')))

# ======================================================================
# BEGIN FOR LOOP — everything until the matching } is processed per file
# ======================================================================
for (filename in filenames) {
  
  filename_sans_ext <- tools::file_path_sans_ext(filename)
  ext <- tools::file_ext(filename)
  
  log_print(paste(Sys.time(), '====== Processing:', filename, '======'))
  
  input_filepath <- win_path(input_dir, filename)
  
  if (tolower(ext) == 'rds') {
    seurat_obj <- readRDS(input_filepath)
  } else if (tolower(ext) == 'rdata') {
    seurat_obj <- load_rdata(input_filepath)
  } else if (tolower(ext) == 'qs') {
    seurat_obj <- qs_read(input_filepath)
  } else {
    stop(paste('unknown filetype:', ext))
  }
  
  # Verify SCT assay exists
  if (!"SCT" %in% names(seurat_obj@assays)) {
    log_print(paste("WARNING: SCT assay not found in", filename))
    log_print(paste("Available assays:", paste(names(seurat_obj@assays), collapse = ", ")))
    next
  }
  DefaultAssay(seurat_obj) <- "SCT"
  
  output_dir <- win_path(wd, opt[['output-dir']], filename_sans_ext)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Detect integration method
  method <- detect_integration_method(seurat_obj)
  log_print(paste(Sys.time(), 'Detected integration method:', method))
  writeLines(method, con = win_path(output_dir, "integration_method.txt"))
  
  # ------------------------------------------------------------------
  # Metadata
  # ------------------------------------------------------------------
  
  log_print(paste(Sys.time(), 'Saving metadata...'))
  
  readr::write_csv(seurat_obj@meta.data,
                   file = win_path(output_dir, "cell_metadata.csv"), col_names = TRUE)
  
  readr::write_csv(data.frame(cell = colnames(seurat_obj)),
                   file = win_path(output_dir, "cells.csv"), col_names = TRUE)
  
  readr::write_csv(data.frame(gene = rownames(seurat_obj)),
                   file = win_path(output_dir, "genes.csv"), col_names = TRUE)
  
  # ------------------------------------------------------------------
  # Gene mapping
  # ------------------------------------------------------------------
  
  log_print(paste(Sys.time(), 'Generating gene mapping table...'))
  tryCatch({
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    
    gene_mapping_table <- getBM(
      filters = "hgnc_symbol",
      attributes = c("hgnc_symbol", "ensembl_gene_id", "gene_biotype"),
      values = rownames(seurat_obj),
      mart = mart, uniqueRows = TRUE)
    
    gene_mapping_table <- data.frame(hgnc_symbol = rownames(seurat_obj)) %>%
      left_join(gene_mapping_table, by = "hgnc_symbol") %>%
      with_groups(hgnc_symbol, mutate, R = row_number()) %>%
      filter(R == 1) %>%
      dplyr::select(-R)
    
    log_print(paste(Sys.time(), 'Saving gene mapping table...'))
    readr::write_csv(gene_mapping_table,
                     file = win_path(output_dir, "gene_mapping.csv"), col_names = TRUE)
  },
  error = function(condition) {
    log_print("WARNING: gene_ensemble not found!!!")
    log_print(paste("Error message:", conditionMessage(condition)))
  })
  
  # ------------------------------------------------------------------
  # Standard workflow (only if no reductions exist at all)
  # ------------------------------------------------------------------
  
  if (length(seurat_obj@reductions) == 0) {
    log_print(paste(Sys.time(), 'Reductions not found, performing SCT standard workflow...'))
    seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)
    ndim <- 30
    seurat_obj <- RunPCA(seurat_obj, npcs = ndim, verbose = FALSE, return.model = TRUE)
    seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:ndim, return.model = TRUE)
    seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:ndim, return.neighbor = TRUE)
    seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
    DefaultAssay(seurat_obj) <- "SCT"
  }
  
  # ------------------------------------------------------------------
  # Reductions
  # ------------------------------------------------------------------
  
  available_reductions <- names(seurat_obj@reductions)
  target_reductions <- REDUCTION_MAP[[method]]
  
  log_print(paste(Sys.time(), 'Available reductions:', paste(available_reductions, collapse = ", ")))
  log_print(paste(Sys.time(), 'Target reductions for', method, ':', paste(target_reductions, collapse = ", ")))
  
  reductions_dir <- win_path(output_dir, 'reductions')
  if (!dir.exists(reductions_dir)) {
    dir.create(reductions_dir, recursive = TRUE)
  }
  
  for (reduction in target_reductions) {
    if (reduction %in% available_reductions) {
      out_name <- gsub("\\.", "_", reduction)
      log_print(paste(Sys.time(), 'Saving reduction:', reduction, '->', paste0(out_name, ".csv")))
      readr::write_csv(
        as.data.frame(seurat_obj@reductions[[reduction]]@cell.embeddings),
        file = win_path(reductions_dir, paste0(out_name, ".csv")),
        col_names = TRUE)
    } else {
      log_print(paste("WARNING: Expected reduction", reduction, "not found!"))
    }
  }
  
  if (!"pca" %in% target_reductions && "pca" %in% available_reductions) {
    log_print(paste(Sys.time(), 'Also saving pre-integration PCA for reference'))
    readr::write_csv(
      as.data.frame(seurat_obj@reductions[["pca"]]@cell.embeddings),
      file = win_path(reductions_dir, "pca.csv"), col_names = TRUE)
  }
  
  # ------------------------------------------------------------------
  # K-nearest neighbors
  # ------------------------------------------------------------------
  
  neighbor_keys <- names(seurat_obj@neighbors)
  log_print(paste(Sys.time(), 'Available neighbors:', paste(neighbor_keys, collapse = ", ")))
  
  nn_key <- NULL
  if (length(neighbor_keys) > 0) {
    if ("SCT.nn" %in% neighbor_keys) {
      nn_key <- "SCT.nn"
    } else {
      nn_matches <- grep("\\.nn$", neighbor_keys, value = TRUE)
      if (length(nn_matches) > 0) {
        nn_key <- nn_matches[1]
      } else {
        nn_key <- neighbor_keys[1]
      }
    }
  }
  
  if (!is.null(nn_key)) {
    log_print(paste(Sys.time(), 'Using neighbor key:', nn_key))
    
    knn_param <- ncol(seurat_obj@neighbors[[nn_key]]@nn.idx)
    n_cells <- ncol(seurat_obj)
    row_idx <- col_idx <- X <- vector("numeric", length = n_cells * knn_param)
    
    for (i in seq(n_cells)) {
      row_idx[(knn_param * i - (knn_param - 1)):(knn_param * i)] <- i
      col_idx[(knn_param * i - (knn_param - 1)):(knn_param * i)] <- seurat_obj@neighbors[[nn_key]]@nn.idx[i, ]
      X[(knn_param * i - (knn_param - 1)):(knn_param * i)] <- seurat_obj@neighbors[[nn_key]]@nn.dist[i, ]
    }
    knn_dist_mat <- Matrix::sparseMatrix(i = row_idx, j = col_idx, x = X)
    
    log_print(paste(Sys.time(), 'Saving knn distances...'))
    readr::write_csv(as.data.frame(seurat_obj@neighbors[[nn_key]]@nn.idx),
                     file = win_path(output_dir, "KNN_indices.csv"), col_names = TRUE)
    Matrix::writeMM(knn_dist_mat, file = win_path(output_dir, "KNN_distances.mtx"))
    
    umap_key <- target_reductions[grep("^umap", target_reductions)]
    if (length(umap_key) > 0 && umap_key %in% available_reductions) {
      if (!is.null(seurat_obj@reductions[[umap_key]]@misc$model)) {
        embedding_key <- target_reductions[!grepl("^umap", target_reductions)]
        tryCatch({
          umap_reembed <- uwot::umap(
            X = seurat_obj@reductions[[embedding_key]]@cell.embeddings,
            n_neighbors = seurat_obj@reductions[[umap_key]]@misc$model$n_neighbors,
            n_components = 2, metric = "cosine",
            n_epochs = seurat_obj@reductions[[umap_key]]@misc$model$n_epochs,
            init = "spectral",
            learning_rate = seurat_obj@reductions[[umap_key]]@misc$model$alpha,
            min_dist = 0.1,
            local_connectivity = seurat_obj@reductions[[umap_key]]@misc$model$local_connectivity,
            nn_method = "annoy",
            negative_sample_rate = seurat_obj@reductions[[umap_key]]@misc$model$negative_sample_rate,
            ret_extra = c("fgraph"), n_threads = 2,
            a = seurat_obj@reductions[[umap_key]]@misc$model$a,
            b = seurat_obj@reductions[[umap_key]]@misc$model$b,
            search_k = seurat_obj@reductions[[umap_key]]@misc$model$search_k,
            approx_pow = seurat_obj@reductions[[umap_key]]@misc$model$approx_pow,
            verbose = TRUE)
          log_print(paste(Sys.time(), 'Saving umap connectivity matrix...'))
          Matrix::writeMM(umap_reembed$fgraph, file = win_path(output_dir, "UMAP_connectivity_matrix.mtx"))
        },
        error = function(condition) {
          log_print(paste("WARNING: UMAP re-embedding failed:", conditionMessage(condition)))
        })
      }
    }
  } else {
    log_print(paste(Sys.time(), 'No neighbors found, skipping KNN export'))
  }
  
  # ------------------------------------------------------------------
  # Count matrices — SeuratObject v5 API
  # ------------------------------------------------------------------
  
  counts_dir <- win_path(output_dir, 'sct-counts')
  if (!dir.exists(counts_dir)) {
    dir.create(counts_dir, recursive = TRUE)
  }
  
  get_sct_layer <- function(seurat_obj, layer_name) {
    tryCatch({
      mat <- GetAssayData(seurat_obj, assay = "SCT", layer = layer_name)
      return(mat)
    }, error = function(e1) {
      tryCatch({
        mat <- seurat_obj[["SCT"]][[layer_name]]
        return(mat)
      }, error = function(e2) {
        return(NULL)
      })
    })
  }
  
  counts_mat <- get_sct_layer(seurat_obj, "counts")
  if (!is.null(counts_mat) && nrow(counts_mat) > 0) {
    log_print(paste(Sys.time(), 'Saving SCT counts...',
                    nrow(counts_mat), 'genes x', ncol(counts_mat), 'cells'))
    Matrix::writeMM(counts_mat, file = win_path(counts_dir, "SCT_counts.mtx"))
  } else {
    log_print("WARNING: SCT counts not found or empty")
  }
  
  data_mat <- get_sct_layer(seurat_obj, "data")
  if (!is.null(data_mat) && nrow(data_mat) > 0) {
    log_print(paste(Sys.time(), 'Saving SCT data (normalized)...',
                    nrow(data_mat), 'genes x', ncol(data_mat), 'cells'))
    Matrix::writeMM(data_mat, file = win_path(counts_dir, "SCT_data.mtx"))
  } else {
    log_print("WARNING: SCT data not found or empty")
  }
  
  scale_mat <- get_sct_layer(seurat_obj, "scale.data")
  if (!is.null(scale_mat) && nrow(scale_mat) > 0) {
    log_print(paste(Sys.time(), 'Saving SCT scale.data...',
                    nrow(scale_mat), 'genes x', ncol(scale_mat), 'cells'))
    scale_df <- as.data.frame(as.matrix(scale_mat))
    readr::write_csv(scale_df,
                     file = win_path(counts_dir, "SCT_scale_data.csv"), col_names = TRUE)
    readr::write_csv(data.frame(gene = rownames(scale_mat)),
                     file = win_path(counts_dir, "SCT_scale_data_genes.csv"), col_names = TRUE)
  } else {
    log_print("WARNING: SCT scale.data not found or empty")
  }
  
  # Free memory before next file
  rm(seurat_obj)
  gc()
  
  log_print(paste(Sys.time(), 'Finished processing:', filename))
  log_print("------")
  
  # ======================================================================
  # END FOR LOOP
  # ======================================================================
}

# ----------------------------------------------------------------------
# End
# ----------------------------------------------------------------------

end_time <- Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()