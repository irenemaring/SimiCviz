# --- Generic function definitions ---

#' Compute activity scores for an object
#'
#' S4 generic for computing TF activity scores (AUC-like metrics).
#'
#' @param object An object that implements a \code{compute_auc} method.
#' @param ... Additional method-specific arguments.
#' @return Object with computed activity scores.
#' @export
setGeneric("compute_auc", function(object, ...) standardGeneric("compute_auc"))

#' Extract activity scores from an object
#'
#' S4 generic for retrieving computed TF activity scores.
#'
#' @param object An object that implements a \code{get_auc} method.
#' @param ... Additional method-specific arguments.
#' @return Extracted activity scores.
#' @export
setGeneric("get_auc", function(object, ...) standardGeneric("get_auc"))

# --- AUCProcessor class definition ---

#' AUCProcessor class
#'
#' Efficient calculator for TF activity scores (AUC) based on GRN weights
#' and cell-specific expression.
#'
#' @slot weights list with weight matrices per label (TFs x targets)
#' @slot expression data.frame or matrix (genes x cells)
#' @slot cell_labels data.frame with columns 'cell' and 'label'
#' @slot target_ids character vector of target gene IDs
#' @slot tf_ids character vector of TF gene IDs
#' @slot auc_results matrix to store computed activity scores (cells x TFs)
#' @slot parallel_params list with parallelization settings
#'
#' @export
setClass(
  "AUCProcessor",
  slots = c(
    weights           = "data.frame",
    expression        = "ANY",
    cell_labels       = "data.frame",
    target_ids        = "character",
    tf_ids            = "character",
    auc_results       = "ANY",  # Changed to "ANY" to allow initialziation
    parallel_params   = "list"
  ),
  prototype = list(
    weights           = data.frame(),
    expression        = NULL,
    cell_labels       = data.frame(cell = character(), label = integer()),
    target_ids        = character(),
    tf_ids            = character(),
    auc_results       = data.frame(), 
    parallel_params   = list(n_cores = 1, backend = "sequential")
  )
)

#' Initialize AUCProcessor
#'
#' @param weights list with weight matrices per label (if simic input) or data.frame with columns: tf, target, weight, label, tf (if data.frame input). adj_p_val (optional, for p-value filtering 'qc_type' = "p_value").
#' @param expression expression matrix (cells x cells) or path to file.
#'   **IMPORTANT**: Expression must be in genes x cells format (genes as rows, cells as columns),
#'   matching Seurat/SingleCellExperiment conventions.
#' @param cell_labels data.frame with cell-to-label mapping
#' @param qc_type character type of quality metric. Values: "adj_r2" or other column name (pvalue, adj_p_val, etc.) provided in weights in data.frame format. If "adj_r2", targets with R2 below threshold are filtered out. If other (i.e. "p_value"), targets with qc_metric above threshold are filtered out.
#' @param qc_threshold threshold for filtering targets based on qc_metric (default: 0.7 for R2)
#' @param adj_r2_list (optional) A list of R2 values per label per target. Only used if weights are provided in list format and qc_type is "adj_r2". 
#' @param n_cores number of cores for parallelization (default: 1)
#' @param backend parallelization backend ("sequential", "multicore", "multisession")
#'
#' @return AUCProcessor object
#' @export
AUCProcessor <- function(weights,
                        expression,
                        cell_labels,
                        qc_type = NULL,
                        qc_threshold = 0.7,
                        adj_r2_list = NULL,
                        n_cores = 1,
                        backend = "sequential"
 ) {
  
  if (is.null(weights)) {
    message("`weights` is required to construct a SimiCvizExperiment.")
    stop("Please provide a valid `weights` list or data.frame.")
  } else if (!class(weights) %in% c("data.frame" ,"list")) {
    stop("`weights` must be a dataframe or a list with one element per label.")
  } else {

    if (inherits(weights, "list")) {
      message("`weights` provided as a list, converting to list data.frame")
      weights_to_use <- weights

      if (!is.null(qc_type)){

        if (qc_type == "adj_r2") {
         # Filter weights by adjusted R2 if provided separatedly
         message("Filtering weights by adjusted R2 threshold...")
          weights_to_use <- lapply(names(weights), function(label) {
                  weights_label_mat <- weights[[label]]
                  adj_r2_label <- adj_r2_list[[label]]
                  weight_filtered <- weights_label_mat[, adj_r2_label >= qc_threshold, drop = FALSE]
                  return(weight_filtered)
          })

        names(weights_to_use) <- names(weights)
        }
      } else {
        # If no quality metric provided, use all weights
        message("No quality_metric provided in weights list; using all targets.")
      }
        weight_df <- .convert_weights_list_to_df(weights_to_use)
      
  } else if (inherits(weights, "data.frame")) {
      message("`weights` provided as a data.frame")
      weight_df <- weights
      # Filter weights by quality metric pvalue if provided in column.
      if (qc_type != "adj_r2" && qc_type %in% colnames(weight_df) ){
        message(sprintf("Filtering weights by %s < threshold %g",qc_type, qc_threshold))
        mask_pval <- which(weight_df[[qc_type]] <= qc_threshold) # Keep below threshold
        weight_df <- weight_df[mask_pval, , drop = FALSE]
      } else if (qc_type == "adj_r2" && "adj_r2" %in% colnames(weight_df)) {
        message("Filtering weights by adjusted R2 threshold...")
        mask_adj_r2 <- which(weight_df$adj_r2 >= qc_threshold)
        weight_df <- weight_df[mask_adj_r2, , drop = FALSE]
      } else {
        message(sprintf("No %s column found in weights data.frame; skipping filtering.",qc_type))
      }
  }
  }
  # Check cell_label input
  if (is.null(cell_labels) || !is.data.frame(cell_labels)) {
    stop("`cell_labels` must be a data.frame with columns 'cell' and 'label'.")
  }
  # Load expression if path provided
  expr_mat <- load_expression_matrix(expression)
  
  # Ensure expression is in cells x genes format and transpose if necessary
  expr_mat <- .ensure_genes_x_cells_format(expr_mat, cell_labels)
  
  # Subset both expression matrix and cell_labels to matching cells
  common_cells <- intersect(colnames(expr_mat), cell_labels$cell)
  if (length(common_cells) == 0) {
    stop(
      "No cells match between expression matrix (colnames) and cell_labels. ",
      "Expression matrix columns must contain cell identifiers matching cell_labels$cell."
    )
  }
  
  if (length(common_cells) < nrow(cell_labels)) {
    warning(sprintf(
      "%d cells in cell_labels not found in expression matrix. ",
      nrow(cell_labels) - length(common_cells)
    ))
  }
  
  if (length(common_cells) < ncol(expr_mat)) {
    warning(sprintf(
      "%d cells in expression matrix not found in cell_labels. Subsetting to matching cells.",
      ncol(expr_mat) - length(common_cells)
    ))
  }
  
  # Subset expression matrix to common cells
  expr_mat <- expr_mat[, common_cells, drop = FALSE]
  
  # Subset cell_labels to common cells
  cell_labels <- cell_labels[cell_labels$cell %in% common_cells, ]
  cell_labels <- cell_labels[match(colnames(expr_mat), cell_labels$cell), ]
  rownames(cell_labels) <- NULL
  
  # Validate inputs
  .validate_auc_processor_inputs(weight_df, expr_mat, cell_labels)
  
  # Extract TF and target IDs from weights
  tf_ids <- unique(weight_df$tf)
  target_ids <- unique(weight_df$target)
  
  # Validate that targets exist in expression matrix
  missing_targets <- setdiff(target_ids, rownames(expr_mat))
  if (length(missing_targets) > 0) {
    warning(sprintf(
      "Missing %d target genes in expression matrix. They will be excluded.",
      length(missing_targets)
    ))
    target_ids <- intersect(target_ids, rownames(expr_mat))
  }
  # Subset the expression matrix to only include target and tf genes
  genes_to_use <- unique(c(target_ids, tf_ids))
  missing_genes <- setdiff( rownames(expr_mat),genes_to_use)
  if (length(missing_genes) > 0) {
    message(sprintf(
      "Subsetting expression matrix to genes in GRNs. Number of tfs: %d and Number of targets: %d ",
      length(tf_ids),length(target_ids)))
      expr_mat <- expr_mat[genes_to_use, , drop = FALSE]
  }
  
  # Prepare parallelization params
  parallel_params <- list(
    n_cores = n_cores,
    backend = backend
  )
  
  # # Initialize empty quality_metric if not provided
  # if (is.null(quality_metric)) {
  #   quality_metric <- lapply(weights, function(w) rep(1.0, ncol(w)))
  #   warning("No quality_metric provided; using all targets.")
  # }
  
  new(
    "AUCProcessor",
    weights           = weight_df,
    expression        = expr_mat,
    cell_labels       = cell_labels,
    target_ids        = target_ids,
    tf_ids            = tf_ids,
    auc_results       = data.frame(),
    parallel_params   = parallel_params
  )
}

# --- Show method ---

#' @param object A \code{\linkS4class{AUCProcessor}} object.
#' @rdname AUCProcessor-class
#' @export
setMethod("show", "AUCProcessor", function(object) {
  n_labels <- length(object@weights)
  n_tfs <- length(object@tf_ids)
  n_targets <- length(object@target_ids)
  n_cells <- nrow(object@cell_labels)
  
  cat("An object of class AUCProcessor\n")
  cat(sprintf(" %d label(s), %d TF(s), %d target(s), %d cells\n",
              n_labels, n_tfs, n_targets, n_cells))
  
  if (nrow(object@auc_results) > 0) {
    cat(sprintf(" Computed AUC scores: %d cells x %d TFs\n",
                nrow(object@auc_results), ncol(object@auc_results)))
  } else {
    cat(" Computed AUC scores: none (run compute_auc() first)\n")
  }
  
  cat(sprintf(" Parallelization: %d cores (%s backend)\n",
              object@parallel_params$n_cores,
              object@parallel_params$backend))
  
  invisible(object)
})

# --- Main computation method ---

#' Compute activity scores for all cells
#'
#' Efficiently calculates TF activity scores for each cell using its
#' corresponding label-specific GRN. Uses BiocParallel for parallelization.
#'
#' @param object A \code{\linkS4class{AUCProcessor}} object.
#' @param sort_by sorting criterion ("expression" or "weight")
#' @param select_top_k keep only top K targets per TF (NULL = use all)
#' @param percent_of_target percentage of targets to use (0-1)
#' @param verbose print progress messages
#'
#' @return A \code{\linkS4class{AUCProcessor}} object with computed scores
#'   in the \code{auc_results} slot.
#' @export
setMethod(
  "compute_auc",
  signature(object = "AUCProcessor"),
  function(object,
           sort_by = "expression",
           select_top_k = NULL,
           percent_of_target = 1.0,
           verbose = TRUE) {
    
    if (!sort_by %in% c("expression", "weight")) {
      stop("sort_by must be one of: 'expression' or 'weight'")
    }
    ts <- Sys.time()
    
    # Map cells to labels for efficient lookup
    cell_to_label <- stats::setNames(object@cell_labels$label, object@cell_labels$cell)
    
    # Extract subset of expression matrix (targets only), transpose to cells x targets
    expr_targets <- Matrix::t(object@expression[object@target_ids, , drop = FALSE])
    
    # Precompute Euclidean norm of target expression
    target_norm <- sqrt(Matrix::colSums(expr_targets^2))
    names(target_norm) <- colnames(expr_targets)
    # Pre-split weights once per label
    weights_by_label <- split(object@weights, object@weights$label)

    # Work units are CELLS (not labels) -> works well even with 1 label
    all_cells <- rownames(expr_targets)
    if (is.null(all_cells)) stop("Cell IDs are missing in expression matrix.")
    n_cells <- length(all_cells)
    n_labels <- length(unique(object@cell_labels$label))
    
    # ---- Adaptive backend selection ----
    req_cores <- max(1L, as.integer(object@parallel_params$n_cores))
    user_backend <- object@parallel_params$backend
    is_windows <- identical(.Platform$OS.type, "windows")
    
    # Default behavior if backend was not explicitly valid
    if (is.null(user_backend) || !(user_backend %in% c("sequential", "multicore", "multisession"))) {
      user_backend <- "sequential"
    }

    # Auto strategy:
    # - if user asked sequential or cores==1 -> serial
    # - else prefer multicore on non-Windows
    # - on Windows, or if user explicitly asks multisession -> SnowParam
    if (req_cores <= 1L || user_backend == "sequential") {
      effective_backend <- "sequential"
    } else if (user_backend == "multisession") {
      effective_backend <- "multisession"
    } else if (user_backend == "multicore") {
      if (is_windows) {
        effective_backend <- "multisession"
      } else {
        effective_backend <- "multicore"
      }
    } else {
      # defensive fallback
      effective_backend <- if (is_windows) "multisession" else "multicore"
    }
    # Small-workload fallback to serial
    # heuristic: if very few cells, parallel overhead dominates
    if (n_cells < 100L || n_cells < (50L * req_cores)) {
      effective_backend <- "sequential"
    }
    # ---- Build BiocParallel param ----
    if (effective_backend == "multicore") {
      BPPARAM <- BiocParallel::MulticoreParam(
        workers = req_cores,
        progressbar = verbose
      )
    } else if (effective_backend == "multisession") {
      BPPARAM <- BiocParallel::SnowParam(
        workers = req_cores,
        type = "SOCK",
        progressbar = verbose
      )
    } else {
      BPPARAM <- BiocParallel::SerialParam(progressbar = verbose)
    }
    
    n_workers <- BiocParallel::bpnworkers(BPPARAM)
    # ---- Chunking strategy (parallelize over chunks of cells) ----
    # Bigger chunks reduce scheduler overhead; especially important for SnowParam
    if (effective_backend == "multisession") {
      chunk_size <- max(50L, ceiling(n_cells / (n_workers * 2L)))
    } else if (effective_backend == "multicore") {
      chunk_size <- max(20L, ceiling(n_cells / (n_workers * 4L)))
    } else {
      chunk_size <- n_cells
    }
    chunk_size <- min(chunk_size, n_cells)
    if (verbose) {
      message("Starting AUC computation...")
      message(sprintf("  Sorting by: %s", sort_by))
      message(sprintf("  Targets: %.0f%% of available", percent_of_target * 100))
      message(sprintf("  Labels: %d | Cells: %d", n_labels, n_cells))
      message(sprintf("  Backend requested: %s | backend used: %s | workers: %d",
                      object@parallel_params$backend, effective_backend, n_workers))
      if (effective_backend == "multisession") {
        message("  Note: multisession (SOCK) has worker startup/serialization overhead;")
        message("        best for larger datasets or Windows compatibility.")
      }
    }
    # Capture all required variables explicitly for SOCK workers
    .cell_to_label   <- cell_to_label
    .expr_targets    <- expr_targets
    .weights_by_label <- weights_by_label
    .target_norm     <- target_norm
    .tf_ids          <- object@tf_ids
    .sort_by         <- sort_by
    .select_top_k    <- select_top_k
    .percent_of_target <- percent_of_target
# ---- Cell scorer ----
    compute_one_cell <- function(i) {
      cell_id <- rownames(.expr_targets)[i]
      label   <- .cell_to_label[[cell_id]]

      full <- rep(NA_real_, length(.tf_ids))
      names(full) <- .tf_ids

      if (!is.null(label)) {
        weight_subset <- .weights_by_label[[as.character(label)]]
        if (!is.null(weight_subset) && nrow(weight_subset) > 0) {
          expr_row <- .expr_targets[i, , drop = TRUE]
          names(expr_row) <- colnames(.expr_targets)
          out <- .compute_cell_auc(
            expr_row         = expr_row,
            weight_df        = weight_subset,
            target_norm      = .target_norm,
            sort_by          = .sort_by,
            select_top_k     = .select_top_k,
            percent_of_target = .percent_of_target
          )
        full[names(out)] <- out
        }
      }
      
      data.frame(
        cell = cell_id,
        as.data.frame(t(full)),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
    
    # ---- Run + fallback ----
    all_idx <- seq_len(nrow(expr_targets))
    idx_chunks <- split(all_idx, ceiling(seq_along(all_idx) / chunk_size))

    auc_chunks <- tryCatch(
      {
        suppressWarnings(
          BiocParallel::bplapply(idx_chunks, function(idx) {
            library(Matrix)
            library(SimiCviz)
            rows <- lapply(idx, compute_one_cell)
            do.call(rbind, rows)
          }, BPPARAM = BPPARAM)
        )
      },
      error = function(e) {
        warning(sprintf(
          "Parallel processing failed: %s. Falling back to sequential processing...",
          e$message
        ))
        lapply(idx_chunks, function(idx) {
          rows <- lapply(idx, compute_one_cell)
          do.call(rbind, rows)
        })
      }
    )
    
    # ---- Combine + preserve order ----
    auc_df <- do.call(rbind, auc_chunks)
    rownames(auc_df) <- auc_df$cell
    auc_result <- auc_df[, -1, drop = FALSE]
    auc_result <- auc_result[object@cell_labels$cell, , drop = FALSE]
    
    object@auc_results <- auc_result
    
    te <- Sys.time()
    elapsed <- difftime(te, ts, units = "secs")
    
    if (verbose) {
      message(sprintf("AUC computation completed in %.2f seconds!", as.numeric(elapsed)))
      message(sprintf("  Result: %d cells x %d TFs", nrow(auc_result), ncol(auc_result)))
      if (effective_backend == "multisession") {
        message("  Tip: if runtime is dominated by setup overhead, try fewer workers (e.g., 2-4)")
        message("       or use backend='multicore' on Linux/macOS.")
      }
    }
    
    object
  }
)

# --- Extract results ---

#' Get computed activity scores
#'
#' @param object A \code{\linkS4class{AUCProcessor}} object.
#' @param format Output format: \code{"wide"} for cells x TF table, or
#'   \code{"long"} for long format with columns \code{cell}, \code{tf},
#'   \code{score}, and \code{label}.
#'
#' @return A data.frame with activity scores in the requested format.
#' @export
setMethod(
  "get_auc",
  signature(object = "AUCProcessor"),
  function(object, format = "matrix") {
    
    if (nrow(object@auc_results) == 0) {
      stop("No AUC results computed yet. Run compute_auc() first.")
    }
    
    if (format == "wide") {
      result <- object@auc_results
    } else if (format == "long") {
      # Long format: cell, tf, score, label
      result <- reshape2::melt(as.matrix(object@auc_results), 
                               varnames = c("cell", "tf"),
                               value.name = "score")
        label_map <- stats::setNames(object@cell_labels$label, object@cell_labels$cell)
        result$label <- label_map[result$cell]

    } else {
      stop("format must be 'matrix' or 'dataframe'")
    }
    result
}
)
# --- Internal helper functions ---

#' Validate expression matrix is in genes x cells format
#'
#' Strictly enforces that expression matrix is in genes x cells format
#' with column names matching cell identifiers in cell_labels.
#'
#' @param expr_mat expression matrix
#' @param cell_labels data.frame with cell identifiers in 'cell' column
#'
#' @return expression matrix (genes x cells) with validated dimensions
.ensure_genes_x_cells_format <- function(expr_mat, cell_labels) {
  # Load cell labels and extract cell IDs
  cell_labels <- load_cell_labels(cell_labels)
  cell_ids <- cell_labels$cell
  n_cells_expected <- length(cell_ids)
  
  # Check that expression matrix has column names
  if (is.null(colnames(expr_mat))) {
    stop(
      "Expression matrix must have column names corresponding to cell identifiers. ",
      "Column names are missing."
    )
  }
  
  # Check that column names match cell_labels
  expr_colnames <- colnames(expr_mat)
  expr_rownames <- rownames(expr_mat)
  n_cols <- length(expr_colnames)
  n_rows <- nrow(expr_mat)
  
  # Validate: colnames should match cells, not be cell count
  n_matching_cells <- sum(expr_colnames %in% cell_ids)
  matching_cells <- which(expr_colnames %in% cell_ids)
  if (n_matching_cells == 0) {
    # No matches found, likely wrong orientation or mismatched identifiers
    # Check if rownames match cell IDs (possible transposed format)
    if (!is.null(expr_rownames)) {
      n_rownames_match <- sum(expr_rownames %in% cell_ids)
      if (n_rownames_match > n_matching_cells) {
        # transposed format detected
        expr_mat <- Matrix::t(expr_mat)
        message("Transposing expression matrix to genes x cells format based on rownames matching cell IDs")
        n_matching_cells <- sum(expr_rownames %in% cell_ids)
        n_cols <- ncol(expr_mat)
        matching_cells <- which(expr_rownames %in% cell_ids)
        } else {
          stop(
            "Expression matrix column names do not match any cell identifiers in cell_labels and no rownames were found ",
            sprintf("Expression has %d columns, cell_labels has %d cells. ", n_cols, n_cells_expected)
          )
        }
      } else {
          stop(
      "Expression matrix column/row names do not match any cell identifiers in cell_labels. ",
          sprintf("Expression has %d columns, cell_labels has %d cells. ", n_cols, n_cells_expected),
          sprintf("First few column names: %s", paste(head(expr_colnames, 3), collapse = ", ")),
          sprintf("First few row names: %s", paste(head(expr_rownames, 3), collapse = ", ")),
          sprintf("First few cell IDs: %s", paste(head(cell_ids, 3), collapse = ", ")),
    )
  }
 }
  
  if (n_cols > n_matching_cells) {
    warning(sprintf(
      "Only %d/%d (%.1f%%) of cell_labels found in expression matrix",
      n_matching_cells, n_cols, 100 * n_matching_cells / n_cols
    ))
  }
  # Return the expression matrix in genes x cells format subset to matching cells
  expr_mat <- expr_mat[, matching_cells, drop = FALSE]
  
  return(expr_mat)
}




#' Validate AUCProcessor inputs
#'
#' @param weights list of weight matrices
#' @param expression expression matrix (genes x cells)
#' @param cell_labels cell-to-label mapping
.validate_auc_processor_inputs <- function(weights, expression, cell_labels) {
  # Check weights
  if (!is.data.frame(weights) || length(weights) == 0) {
    stop("weights must be a non-empty data.frame")
  }
  
  # Check all weights have at least tf, target, weight, label columns
  required_cols <- c("tf", "target", "weight", "label")
  if (!all(required_cols %in% colnames(weights))) {
    stop(sprintf("weights data.frame must contain columns: %s", paste(required_cols, collapse = ", ")))
  } else {
    # Check that weight values are numeric
    if (!is.numeric(weights$weight)) {
      stop("weights$weight column must be numeric")
    }
  }

  
  # Check expression is genes x cells
  if (!is.matrix(expression) && !is.data.frame(expression) && !class(expression)[1]== "dgCMatrix") {
    stop("expression must be a dgCMatrix, matrix or data.frame")
  }
  
  if (!is.null(colnames(expression))) {
    n_expr_cells <- ncol(expression)
  } else {
    stop("Expression matrix must have column names corresponding to cell identifiers.")
  }
  
  # Check cell_labels
  if (!is.data.frame(cell_labels)) {
    stop("cell_labels must be a data.frame")
  }
  if (!all(c("cell", "label") %in% colnames(cell_labels))) {
    stop("cell_labels must have columns 'cell' and 'label'")
  }
  # Check that weights label column matches cell_labels
  weight_labels <- unique(weights$label)
  cell_label_values <- unique(cell_labels$label)
  if (!all(cell_label_values %in% weight_labels)) {
    stop("All labels in cell_labels$label must be present in weights$label")
  }
  invisible(TRUE)
}

#' Compute AUC for a single cell
#'
#' @param expr_row expression values for target genes (vector)
#' @param weight_df weight data.frame subset for this cell's label (columns: tf, target, weight)
#' @param target_norm precomputed Euclidean norm of target expression (vector)
#' @param sort_by sorting criterion ("expression" or "weight")
#' @param select_top_k keep only top K targets per TF (NULL = all)
#' @param percent_of_target percentage of targets to use (0-1)
#'
#' @return named vector of AUC scores (one per TF)
.compute_cell_auc <- function(expr_row,
                              weight_df,
                              target_norm,
                              sort_by,
                              select_top_k = NULL,
                              percent_of_target = 1.0) {
  
  # Get unique TFs for this label
  tfs <- unique(weight_df$tf)
  n_tfs <- length(tfs)
  auc_scores <- numeric(n_tfs)
  names(auc_scores) <- tfs
  # For each TF
  for (tf in tfs) {
    
    # Extract targets and weights for this TF
    tf_data <- weight_df[weight_df$tf == tf, ]
    
    # Match targets to expression row
    targets_in_tf <- intersect(tf_data$target, names(target_norm))
    
    if (length(targets_in_tf) == 0) {
      auc_scores[tf] <- NA_real_
      next
    }
    
    # Get expression and weights for these targets
    expr_vals <- expr_row[targets_in_tf]
    weight_vals <- tf_data$weight[match(targets_in_tf, tf_data$target)]
    target_norm_regulon  <- target_norm[match(targets_in_tf, names(target_norm))]
    # Normalize weights
    weights_norm <- weight_vals / target_norm_regulon
    # Determine sorting order
    if (sort_by == "expression") {
      sort_idx <- order(expr_vals, decreasing = TRUE)
    } else {
      sort_idx <- order(weight_vals, decreasing = TRUE)
    }
    
    # Apply sorting
    expr_sorted <- expr_vals[sort_idx]
    weight_sorted <- abs(weights_norm[sort_idx])
    
    # Select percentage of targets
    n_selected <- max(1, round(length(weight_sorted) * percent_of_target))
    
    # Apply top-K if specified
    if (!is.null(select_top_k)) {
      n_selected <- min(n_selected, select_top_k)
    }
    
    # Calculate cumulative sum of top targets
    weight_top <- weight_sorted[seq_len(n_selected)]
    cum_weights <- cumsum(weight_top)
    auc_score <- sum(cum_weights) / (sum(weight_top) * n_selected)
    
    auc_scores[tf] <- auc_score
  }
  
  auc_scores
}
