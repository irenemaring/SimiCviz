#' SimiCviz: Visualization Tools for SimiC Regulatory Network Outputs
#'
#' SimiCviz provides utilities to:
#' \itemize{
#'   \item Read SimiC or SimiCPipeline output tables from CSV and pickle files.
#'   \item Organize weights, AUCs, and metadata into a coherent object.
#'   \item Visualize regulatory networks and performance metrics.
#'   \item Export plots (PDF) and processed data tables (CSV) into a
#'         reproducible directory hierarchy.
#' }
#'
#' The package is designed to be lightweight and focused on visualization,
#' leaving modeling and inference to SimiCPipeline itself.
#'
#' @docType package
#' @name SimiCviz
NULL

# --- S4 class definition -------------------------------------------------

#' SimiCvizExperiment class
#'
#' Lightweight S4 container for SimiC visualization data.
#'
#' @slot weights list with weight matrices and adjusted R2 per label.
#' @slot auc list containing the collected AUC data. The canonical element is
#'   \code{$collected}: a data.frame with cells in rows, TFs in columns, and
#'   an optional \code{label} column. This format is shared across all GRN
#'   methods that produce per-cell TF activity scores.
#' @slot cell_labels data.frame with columns \code{cell} (character) and
#'   \code{label} (integer).
#'   Provides the mapping between cell identifiers and phenotype labels.
#'   Required for label-specific AUC visualizations.
#' @slot tf_ids character vector of TF identifiers.
#' @slot target_ids character vector of target gene identifiers.
#' @slot label_names named character vector mapping integer labels to display names.
#' @slot colors named character vector mapping integer labels to colors.
#' @slot meta list with arbitrary metadata.
#' @export
setClass(
  "SimiCvizExperiment",
  slots = c(
    weights     = "list",
    auc         = "list",
    cell_labels = "data.frame",
    tf_ids      = "character",
    target_ids  = "character",
    label_names = "character",
    colors      = "character",
    meta        = "list"
  ),
  prototype = list(
    weights     = list(),
    auc         = list(),
    cell_labels = data.frame(cell = character(), label = integer(),
                             stringsAsFactors = FALSE),
    tf_ids      = character(),
    target_ids  = character(),
    label_names = character(),
    colors      = character(),
    meta        = list()
  )
)

# --- show method ----------------------------------------------------------

#' @rdname SimiCvizExperiment-class
#' @export
setMethod("show", "SimiCvizExperiment", function(object) {
  n_labels   <- length(object@weights)
  n_tfs      <- length(object@tf_ids)
  n_targets  <- length(object@target_ids)
  n_cells    <- nrow(object@cell_labels)
  has_auc_df <- !is.null(object@auc$collected)
  has_perlab <- !is.null(object@auc$per_label)

  cat("An object of class SimiCvizExperiment\n")
  cat(sprintf(" %d label(s), %d TF(s), %d target(s)\n", n_labels, n_tfs, n_targets))

  # Weights
  if (n_labels > 0L) {
    dims <- vapply(object@weights, function(m) paste(nrow(m), "x", ncol(m)), character(1))
    cat(sprintf(" Weights: %d matri%s [%s]\n",
                n_labels, ifelse(n_labels == 1L, "x", "ces"),
                paste(names(object@weights), dims, sep = ": ", collapse = ", ")))
  } else {
    cat(" Weights: none\n")
  }

  # AUC
  if (has_auc_df) {
    auc_df <- object@auc$collected
    auc_tfs <- setdiff(colnames(auc_df), "label")
    cat(sprintf(" AUC: collected (%d cells x %d TFs)", nrow(auc_df), length(auc_tfs)))
    if ("label" %in% colnames(auc_df)) cat(" + label column")
    cat("\n")
  } else if (has_perlab) {
    n_pl <- length(object@auc$per_label)
    cat(sprintf(" AUC: %d per-label raw matri%s (not yet collected)\n",
                n_pl, ifelse(n_pl == 1L, "x", "ces")))
  } else {
    cat(" AUC collected: none\n")
  }
  # if (!has_auc_df && has_perlab){
  #   auc_list <- object@auc$labelled
  #   auc_df <- auc_list[[1]]
  #   auc_tfs <- setdiff(colnames(auc_list[[1]]), "label")
  #   cat(sprintf(" AUC: labelled list of %d (%d cells x %d TFs)", length(auc_list), nrow(auc_df), length(auc_tfs)))
  #   cat("\n")
  # }

  # Cell labels
  if (n_cells > 0L) {
    ul <- sort(unique(object@cell_labels$label))
    cat(sprintf(" Cell labels: %d cells across %d label(s) [%s]\n",
                n_cells, length(ul), paste(ul, collapse = ", ")))
  } else {
    cat(" Cell labels: none\n")
  }

  # Label names
  if (length(object@label_names) > 0L) {
    ln <- paste(names(object@label_names), object@label_names, sep = " = ", collapse = ", ")
    cat(sprintf(" Label names: %s\n", ln))
  }

  # Colors
  if (length(object@colors) > 0L) {
    cc <- paste(names(object@colors), object@colors, sep = " = ", collapse = ", ")
    cat(sprintf(" Colors: %s\n", cc))
  }

  # TFs & targets (first few)
  if (n_tfs > 0L) {
    shown <- if (n_tfs > 6L) paste0(paste(head(object@tf_ids, 6), collapse = ", "), ", ...") else paste(object@tf_ids, collapse = ", ")
    cat(sprintf(" TFs: %s\n", shown))
  }
  if (n_targets > 0L) {
    shown <- if (n_targets > 6L) paste0(paste(head(object@target_ids, 6), collapse = ", "), ", ...") else paste(object@target_ids, collapse = ", ")
    cat(sprintf(" Targets: %s\n", shown))
  }

  # Meta
  if (length(object@meta) > 0L) {
    cat(sprintf(" Meta keys: %s\n", paste(names(object@meta), collapse = ", ")))
  }

  invisible(object)
})

# --- cell-label loader ---------------------------------------------------

#' Load cell-to-label mapping
#'
#' Reads or constructs a two-column \code{data.frame} (\code{cell}, \code{label})
#' from various input formats.
#'
#' @param x One of:
#'   \describe{
#'     \item{character (length 1)}{Path to a CSV / TSV file.
#'       If the file has columns named \code{cell} and \code{label} they are used
#'       directly. Otherwise the first column is treated as \code{cell} and the
#'       second as \code{label}. A single-column file (or headerless file with one
#'       field) is treated as a label-only vector — see next bullet.}
#'     \item{character / integer / numeric vector}{A vector of labels, one per
#'       cell, in the \strong{same order} as rows of the AUC matrix.
#'       Cell IDs are generated as \code{cell_1, cell_2, \ldots}.}
#'     \item{named vector}{Names are used as cell IDs and values as labels.}
#'     \item{data.frame}{Must contain columns \code{cell} and \code{label},
#'       or exactly two columns (first = cell, second = label).}
#'   }
#'
#' @return A \code{data.frame} with columns \code{cell} (character) and
#'   \code{label} (integer), sorted by \code{cell}.
#' @export
load_cell_labels <- function(x) {
  if (is.data.frame(x)) {
    df <- .parse_cell_labels_df(x)
  } else if (is.character(x) && length(x) == 1L && file.exists(x)) {
    df <- .parse_cell_labels_file(x)
  } else if (is.vector(x) && !is.list(x)) {
    df <- .parse_cell_labels_vector(x)
  } else {
    stop(
      "`cell_labels` must be a file path, a vector of labels, a named vector, ",
      "or a data.frame with columns 'cell' and 'label'."
    )
  }
  df$label <- as.integer(df$label)
  df$cell  <- as.character(df$cell)
  df
}

# --- internal parsers -----------------------------------------------------

.parse_cell_labels_df <- function(df) {
  if (ncol(df) > 2L){
    message("Cell labels file contains extra columns")
    message(paste(colnames(df), collapse = " / "))
  }
  if (all(c("cell", "label") %in% tolower(colnames(df)))) {
    idx_cols <- which( tolower(colnames(df))%in% c("cell", "label"))
    df <- df[, idx_cols, drop = FALSE]
    colnames(df) <- c("cell", "label")
    return(df)
  }
  if (ncol(df) == 1L) {
    # single column → labels only; use rownames or generate IDs
    ids <- if (!is.null(rownames(df)) &&
               !all(rownames(df) == as.character(seq_len(nrow(df))))) {
      rownames(df)
    } else {
      paste0("cell_", seq_len(nrow(df)))
    }
    return(data.frame(cell = ids, label = df[[1]], stringsAsFactors = FALSE))
  }
}

.parse_cell_labels_file <- function(path) {
  # try comma first, then tab
  df <- tryCatch(
    utils::read.csv(path, stringsAsFactors = FALSE),
    error = function(e) utils::read.delim(path, stringsAsFactors = FALSE)
  )
  .parse_cell_labels_df(df)
}

.parse_cell_labels_vector <- function(x) {
  if (!is.null(names(x))) {
    # named vector: names → cell IDs, values → labels
    data.frame(cell = names(x), label = as.integer(x), stringsAsFactors = FALSE)
  } else {
    # unnamed: positional labels
    data.frame(
      cell  = paste0("cell_", seq_along(x)),
      label = as.integer(x),
      stringsAsFactors = FALSE
    )
  }
}

# --- constructor ----------------------------------------------------------

#' Construct a SimiCvizExperiment object
#'
#' @param weights list with weight matrices and adjusted R2.
#' @param auc list containing AUC data. The preferred format is
#'   \code{list(collected = df)} where \code{df} is a cells × TF data.frame.
#'   A plain data.frame is automatically wrapped into this structure.
#'   This canonical format is compatible with any GRN method that produces
#'   per-cell TF activity scores.
#' @param cell_labels Cell-to-label mapping. Accepts any format supported by
#'   \code{\link{load_cell_labels}}: a file path, a vector, a named vector,
#'   or a data.frame.
#'   \strong{If a plain vector is provided} (no names), it must be in the same
#'   row order as the AUC matrix.
#' @param tf_ids character vector of TF identifiers.
#' @param target_ids character vector of target gene identifiers.
#' @param label_names optional named character vector mapping labels to display
#'   names. Length must match the number of unique labels in \code{cell_labels}
#'   (if provided) or \code{weights}.
#' @param colors optional named character vector mapping labels to colors.
#'   Same length requirement as \code{label_names}.
#' @param meta optional list of additional metadata.
#'
#' @return An object of class \code{"SimiCvizExperiment"}.
#' @export
SimiCvizExperiment <- function(weights = NULL,
                               auc = NULL,
                               cell_labels = NULL,
                               tf_ids = NULL,
                               target_ids = NULL,
                               label_names = character(),
                               colors = character(),
                               meta = list()) {
  # --- weights ---
  if (is.null(weights)) {
    weights <- list()
  } else if (!is.list(weights)) {
    stop("`weights` must be a list (e.g. one element per label).")
  }

  # --- auc: normalize to list(collected = df) or list(per_label = auc_list) ---
  if (is.null(auc)) {
    auc <- list()
  } else if (is.data.frame(auc)) {
    # bare data.frame → wrap into canonical structure
    auc <- list(collected = auc)
  } else if (all(class(auc) == "list" & !names(auc) %in% c("collected", "per_label") & identical(names(auc), names(weights)))){
    auc <-  list(per_label = auc)
  } else if (all(class(auc) == "list" & names(auc) == "per_label")){
    message("AUC in list format")
  } else if (all(class(auc) == "list" & names(auc) =="collected")){
    message("AUC in collected format")
  } else {
    stop("`auc` must be a list or a data.frame (cells x TF collected format).")
  }

  # --- cell_labels ---
  if (is.null(cell_labels)) {
    cl_df <- data.frame(cell = character(), label = integer(),
                        stringsAsFactors = FALSE)
  } else {
    cl_df <- load_cell_labels(cell_labels)
  }

  tf_ids     <- if (is.null(tf_ids))     character() else as.character(tf_ids)
  target_ids <- if (is.null(target_ids)) character() else as.character(target_ids)

  # --- infer n_labels ---
  n_labels_cl <- length(unique(cl_df$label))
  n_labels_w  <- length(weights)
  n_labels <- if (n_labels_cl > 0L) n_labels_cl else n_labels_w

  if (n_labels_cl > 0L && n_labels_w > 0L && n_labels_cl != n_labels_w) {
    warning(sprintf(
      "Number of unique labels in `cell_labels` (%d) differs from `weights` length (%d). Using `cell_labels` as reference.",
      n_labels_cl, n_labels_w
    ))
  }

  # --- label_names / colors ---
  norm_lab <- .normalize_label_names(label_names, n_labels = n_labels)
  norm_col <- .normalize_label_colors(colors,     n_labels = n_labels)

  new(
    "SimiCvizExperiment",
    weights     = weights,
    auc         = auc,
    cell_labels = cl_df,
    tf_ids      = tf_ids,
    target_ids  = target_ids,
    label_names = norm_lab,
    colors      = norm_col,
    meta        = meta
  )
}

#' @export
is.SimiCvizExperiment <- function(x) {
  inherits(x, "SimiCvizExperiment")
}

# --- label utilities (mirrors Python behaviour) ---------------------------

.normalize_label_names <- function(x, n_labels = NULL) {
  if (is.null(x) || length(x) == 0) {
    if (!is.null(n_labels) && n_labels > 0L) {
      message("No `label_names` provided; using default 'Label i' naming.")
    }
    return(character())
  }

  labs <- as.character(x)

  if (!is.null(n_labels) && n_labels > 0L && length(labs) != n_labels) {
    stop(sprintf(
      "`label_names` length (%d) must match the number of labels in experiment (%d).",
      length(labs), n_labels
    ))
  }

  if (is.null(names(labs))) {
    names(labs) <- as.character(seq_len(length(labs)) - 1L)
  }

  lab_ids <- as.integer(names(labs))
  if (any(is.na(lab_ids))) {
    stop("`label_names` names must be coercible to integer labels (e.g. '0','1',...).")
  }

  stats::setNames(labs, as.character(lab_ids))
}

.normalize_label_colors <- function(x, n_labels = NULL) {
  if (is.null(x) || length(x) == 0) {
    if (!is.null(n_labels) && n_labels > 0L) {
      message("No `colors` provided; using internal default palette.")
    }
    return(character())
  }

  cols <- as.character(x)

  if (!is.null(n_labels) && n_labels > 0L && length(cols) != n_labels) {
    stop(sprintf(
      "`colors` length (%d) must match the number of labels in experiment (%d).",
      length(cols), n_labels
    ))
  }

  if (is.null(names(cols))) {
    names(cols) <- as.character(seq_len(length(cols)) - 1L)
  }

  lab_ids <- as.integer(names(cols))
  if (any(is.na(lab_ids))) {
    stop("`colors` names must be coercible to integer labels (e.g. '0','1',...).")
  }

  stats::setNames(cols, as.character(lab_ids))
}

#' Set or update label names and colors
#'
#' Validates that the provided vectors match the number of labels
#' in \code{cell_labels} (preferred) or \code{weights}.
#'
#' @param x \code{SimiCvizExperiment}.
#' @param label_names named character vector mapping labels to display names.
#' @param colors optional named character vector mapping labels to colors.
#'
#' @return Modified \code{SimiCvizExperiment} object.
#' @export
setLabelNames <- function(x,
                          label_names,
                          colors = NULL) {
  if (!is.SimiCvizExperiment(x)) {
    stop("setLabelNames: 'x' must be a SimiCvizExperiment.")
  }

  # authoritative label count: cell_labels > weights
  n_cl <- length(unique(x@cell_labels$label))
  n_w  <- length(x@weights)
  n_labels <- if (n_cl > 0L) n_cl else n_w

  x@label_names <- .normalize_label_names(label_names, n_labels = n_labels)

  if (!is.null(colors)) {
    x@colors <- .normalize_label_colors(colors, n_labels = n_labels)
  } else if (length(x@colors) == 0L) {
    x@colors <- character()
  }

  x
}

#' Set or update cell-to-label mapping
#'
#' @param x \code{SimiCvizExperiment}.
#' @param cell_labels Any format accepted by \code{\link{load_cell_labels}}.
#'
#' @return Modified \code{SimiCvizExperiment} object.
#' @export
setCellLabels <- function(x, cell_labels) {
  if (!is.SimiCvizExperiment(x)) {
    stop("setCellLabels: 'x' must be a SimiCvizExperiment.")
  }
  x@cell_labels <- load_cell_labels(cell_labels)

  # re-validate existing label_names / colors against new label count
  n_labels <- length(unique(x@cell_labels$label))
  if (length(x@label_names) > 0L && length(x@label_names) != n_labels) {
    warning(sprintf(
      "Existing `label_names` length (%d) no longer matches the number of labels (%d). Clearing.",
      length(x@label_names), n_labels
    ))
    x@label_names <- character()
  }
  if (length(x@colors) > 0L && length(x@colors) != n_labels) {
    warning(sprintf(
      "Existing `colors` length (%d) no longer matches the number of labels (%d). Clearing.",
      length(x@colors), n_labels
    ))
    x@colors <- character()
  }

  x
}