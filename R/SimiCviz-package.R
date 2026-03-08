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
load_cell_labels <- function(x, ...) {
  if (is.data.frame(x)) {
    df <- .parse_cell_labels_df(x)
  } else if (is.character(x) && length(x) == 1L && file.exists(x)) {
    df <- .parse_cell_labels_file(x, ...)
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
  } else if (ncol(df) == 1L){
    stop("Check cell labels input format!")
  }
}

.parse_cell_labels_file <- function(path, ...) {
  if(!file.exists(path)){
    stop(sprintf("File not found: %s", path))
  }
  # try comma first, then tab
  df <- tryCatch(
    utils::read.table(path, stringsAsFactors = FALSE, ...),
    error = function(e) utils::read.table(path, stringsAsFactors = FALSE, ...)
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
#' @param auc optional list containing AUC data. The preferred format is
#'   \code{list(collected = df)} where \code{df} is a cells × TF data.frame.
#'   A plain data.frame is automatically wrapped into this structure.
#'   This canonical format is compatible with any GRN method that produces
#'   per-cell TF activity scores.
#' @param cell_labels optional data.frame or vector with cell-to-label mapping. 
#'   \strong{If a plain vector is provided} (no cell ids), it must be in the same
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
                               label_names = NULL,
                               colors = NULL,
                               meta = list()) {
  # --- weights ---

  if (is.null(weights)) {
    message("`weights` is required to construct a SimiCvizExperiment.")
    stop("Please provide a valid `weights` list or data.frame.")
  } else if (!class(weights) %in% c("data.frame" ,"list")) {
    stop("`weights` must be a dataframe or a list with one element per label.")
  }else{
    weights_input <- weights
    if (is.data.frame(weights_input)) {
      weights <- .convert_weights_df_to_list(weights_input)
      message("`weights` provided as a data.frame, converting to list format.")
    } else {
      message("`weights` provided as a list.")
      weights <- weights_input
  }
  }

  cl_df <- data.frame(cell = character(), label = integer(),
                             stringsAsFactors = FALSE)
  n_labels_w  <- length(weights)
  n_labels_cl <- 0L
  # --- auc ---
  if (is.null(auc)) {
    auc <- list()
    message("No `auc` provided; AUC-related visualizations will be unavailable.")
  } else if (!class(auc) %in% c("data.frame" ,"list")) {
    stop("`auc` must be a dataframe or a list with one element per label.")
  } else {
    # IF AUC is valid class, THEN CELL LABELS MUST BE PROVIDED
    flag <- FALSE

    # --- cell_labels ---
    if (is.null(cell_labels)) {
      message("`cell_labels` are required if auc is provided to construct a SimiCvizExperiment.")
      stop("Please provide a valid `cell_labels` data.frame or label vector.")
    } else if (any(is.vector(cell_labels), is.data.frame(cell_labels))){

      cl_df <- load_cell_labels(cell_labels)
      n_labels_cl <- length(unique(cl_df$label))

      if (is.vector(cell_labels) && identical(cl_df$cell, paste0("cell_", seq_along(cell_labels)))){
        # Mark flag TRUE if cell_labels was a vector, then cell is cell_1, cell_2, to change afterwards with rownames in auc
        flag <- is.vector(cell_labels) 
        }
        
    } else{
      stop("`cell_labels` must be a data.frame or a vector.")
    }

    auc_input <- auc
    if (is.data.frame(auc_input)) {
    # Check if it is in long format (cell, tf, score) or wide format (cells x TFs)
      if (all(c("cell", "tf", "score") %in% tolower(colnames(auc_input)))) {
        message("AUC in long format (cell, tf, score), converting to wide format (collected).")
        auc_tmp <-.convert_auc_long_to_wide(auc_input)
        if (flag){
          warning("`cell_labels` were provided as a vector without cell IDs!")
          stop("Cannot map cell ID's if auc is in wide format")
        }
      } else{
        message("AUC in wide format (cells x TFs).") 
        auc_tmp <- auc_input
      }
      if (flag){
        message("`cell_labels` were provided as a vector without cell IDs; assuming order matches AUC rows.")
        cl_df$cell <- rownames(auc_tmp) # assign cell IDs from AUC rownames
      }
    } else if (all(class(auc_input) == "list" & identical(names(auc), names(weights)))) {
      # Assume it's a per-label list of AUC data.frames (not yet in wide format / collected)
      message("AUC in per-label list format, converting to wide format using cell_labels (collected).")
      if (flag){
        message("`cell_labels` were provided as a vector without cell IDs; assuming order matches AUC rows.")
        cl_df$cell <- rownames(auc_input[[1]]) # assign cell IDs from AUC rownames
      }
      # Check if auc has the cells in cell_labels
      if (!any(rownames(auc_input[[1]]) %in% (cl_df$cell))){
          warning("No `auc` cell IDs match those in cell_labels. Distribution curves will not be displayed.")
          auc_tmp <- NULL
      } else {
        if(all(rownames(auc_input[[1]]) %in% (cl_df$cell))){
          message("AUC cell IDs match those in cell_labels.")
        } else if (any(rownames(auc_input[[1]]) %in% (cl_df$cell))){
          warning("Some auc cell IDs do not match those in cell_labels. Auc results will be displayed for matching cells only.")
        }
        # Convert to wide format using cell_labels to subset cells.
        auc_tmp <- auc_list_to_df(auc_input, cell_labels = cl_df)
      }
    } else{
      auc_tmp <- NULL
    }

   if (is.null(auc_tmp)) {
      message("AUC data could not be processed; AUC-related visualizations will be unavailable.")
      auc <- NULL
    }else{
      auc <- list(collected = auc_tmp)
    }

  if (!is.null(auc) && !any(rownames(auc$collected) %in% cl_df$cell)){
    stop("AUC cell rownames don't match those in cell_labels.")
  }
  }
  # --- infer n_labels ---
  n_labels <- if (n_labels_w > 0L) n_labels_w else n_labels_cl
  if (n_labels_cl == 0L){
    message("No `cell_labels` provided; using `weights` length as reference for number of labels.")
  }else if (n_labels_cl > 0L && n_labels_w > 0L && n_labels_cl != n_labels_w) {
    warning(sprintf(
      "Number of unique labels in `cell_labels` (%d) differs from `weights` length (%d). Using `weights` as reference.",
      n_labels_cl, n_labels_w
    ))
  }
  # If TFs or targets are not provided, infer from WEIGHT colnames and rownames (weight is mandatory, auc si not)
  tf_ids     <-  unique(c(sapply(weights,rownames)))
  target_ids <-  unique(c(sapply(weights,colnames)))  

  # --- label_names / colors ---
  
  norm_lab <- .normalize_label_names(label_names, n_labels = n_labels, keys = names(weights))
  norm_col <- .normalize_label_colors(colors,     n_labels = n_labels, keys = names(weights))

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

# --- label utilities ---------------------------

.normalize_label_names <- function(label_names, n_labels, keys) {
  if (is.null(label_names) || length(label_names) == 0) {
    if (n_labels > 0L) {
      message("No `label_names` provided; using default 'Label i' naming.")
      label_names <- paste("Label", seq_len(n_labels) - 1L,sep = "_")
    }
  }

  labs <- as.character(label_names)

  if (!is.null(n_labels) && n_labels > 0L && length(labs) != n_labels) {
    stop(sprintf(
      "`label_names` length (%d) must match the number of labels in experiment (%d).",
      length(labs), n_labels
    ))
  }

  if (is.null(names(label_names))) {
    names(labs) <- keys
  } else{
    names(labs) <- names(label_names)
  }

  lab_ids <- as.integer(names(labs))
  if (any(is.na(lab_ids))) {
    stop("`label_names` names must be coercible to integer labels (e.g. '0','1',...).")
  }
  if (!all(lab_ids %in% keys)){
    stop("`label_names` names must match the label identifiers used in `weights`.")
  }

  stats::setNames(labs, as.character(lab_ids))
}

.normalize_label_colors <- function(x, n_labels, keys) {

  default_colors <- c("#5e82bd", "#ed9900", "#008b00","#cd9b9b","#800080","#ff676f")
  
  if (is.null(x) || length(x) == 0) {
    if (!is.null(n_labels) && n_labels > 0L) {
      message("No `colors` provided; using internal default palette.")
      x <- default_colors[seq_len(n_labels)]
    }
  }

  cols <- as.character(unique(x))
  if (n_labels > 0L && length(cols) > n_labels) {
    warning(sprintf(
      "`colors` length (%d) does not match the number of labels in experiment (%d). Selecting first colors: (%s)",
      length(cols), n_labels, paste(cols[1:n_labels], collapse=", ")
    ))
    cols <- cols[1:n_labels]
  } else if (n_labels > 0L && length(cols) < n_labels){
    warning(sprintf(
      "`colors` length (%d) is less than the number of labels in experiment (%d). Adding default extra colors.",
      length(cols), n_labels
    ))
    default_colors_allowed  <- setdiff(default_colors, cols)
    cols  <- c(cols, default_colors_allowed[seq_len(n_labels - length(cols))])
  }

  if (is.null(names(cols))) {
    names(cols) <- keys
  } else{
    names(cols) <- names(x)
  }

  col_ids <- as.integer(names(cols))
  if (any(is.na(col_ids))) {
    stop("`colors` names must be coercible to integer labels (e.g. '0','1',...).")
  }
  if (!all(col_ids %in% keys)){
    stop("`color` names must match the label identifiers used in `weights`.")
  }

  stats::setNames(cols, as.character(col_ids))
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
  keys <- names(x@weights)
  x@label_names <- .normalize_label_names(label_names, n_labels = n_labels, keys = keys)

  if (!is.null(colors)) {
    x@colors <- .normalize_label_colors(colors, n_labels = n_labels, keys = keys)
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