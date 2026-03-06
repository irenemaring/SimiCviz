##  Simple Import Functions ---------------------------------------------- 

#' Read pickle files
#' 
#' @param file path to the pickle file
#' @return object containing the contents of the pickle file
#' @export
read_pickle <- function(file) {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required to read pickle files. Please install it.")
  }
  if (all(sapply(c("pickle$", "pkl$"), grepl, file), USE.NAMES = FALSE)) {
     stop("Expected a pickle file with extension .pickle or .pkl: ")
  }
  if (!file.exists(file)) {
    stop(sprintf("File % does not exist.", file))
  }
  .ensure_python()
  pickle <- reticulate::import("pickle")
  con <- reticulate::import("builtins")$open(file, "rb")
  out <- pickle$load(con)
  con$close()
  results_pickle <- reticulate::py_to_r(out)
  
  return(results_pickle)
}

#' Read SimiCPipeline output weights (pickle file)
#'
#' @param file path to the pickle file containing SimiC weights.
#' @return A list of length n_labels (phenotypes) with GRN weight matrices (TFs x targets).
#' @export
read_weights_pickle <- function(file) {
  weight_dict <- read_pickle(file)  
  if (!all(names(weight_dict) %in% c("weight_dic", "adjusted_r_squared", 
                                   "standard_error", "TF_ids", "query_targets"))) {
    stop("Unexpected structure in weights pickle file: %s. Make sure the weights file is correct", file)
  }
  weight_mat_list <- weight_dict$weight_dic
  weight_mat_out <- list()
  # Remove bias term and add rownames and colnames to each weight matrix
  for (key in names(weight_mat_list)) {
    tmp_weights <-  as.data.frame(weight_mat_list[[key]])
    bias_idx <- nrow(tmp_weights)
    tmp_weights <- tmp_weights[-bias_idx, , drop = FALSE]
    rownames(tmp_weights) <- weight_dict$TF_ids
    colnames(tmp_weights) <- weight_dict$query_targets
    
    weight_mat_out[[key]] <- tmp_weights
  }
  
  return(weight_mat_out)
}

#' Read SimiCPipeline output AUC matrices (pickle file)
#'
#' Returns a list of matrices (one per label), each with cells in rows and
#' TFs in columns.
#'
#' @param file path to a pickle file containing per-label AUC matrices.
#' @return A named list of data.frames / matrices (one per label).
#' @export
read_auc_pickle <- function(file) {
  auc_list <- read_pickle(file)
  if (!is.list(auc_list)) {
    stop("Expected a list of AUC matrices in pickle file: ", file)
  }
  auc_list
}

#' Load a collected AUC CSV into a data.frame
#'
#' Reads a cells × TF CSV file (as written by \code{\link{save_collected_auc}}
#' or SimiCPipeline). This is the canonical import for any GRN method that
#' produces per-cell TF activity scores.
#'
#' @param file Path to a CSV file with cells in rows and TFs in columns.
#' @param ... Additional arguments passed to \code{\link[utils]{read.csv}}.
#' @return A data.frame with cell IDs as row names.
#' @export
load_collected_auc <- function(file, ...) {
  if (!file.exists(file)) {
    stop("File does not exist: ", file)
  }
  df <- utils::read.csv(file, header = TRUE, row.names = 1,
                         stringsAsFactors = FALSE, ...)
  df
}

# ---- Legacy / CSV readers -----------------------------------------------

#' Read Long-style TF-target weights from CSV (compatibility)
#'
#' @param file path to a CSV file containing TF-target weights.
#'   Expected minimal columns: \code{tf}, \code{target}, \code{weight}.
#'   Additional columns such as \code{label}, etc. are preserved.
#' @param ... additional arguments passed to [utils::read.csv()].
#'
#' @return data.frame with weights.
#' @export
read_weights_csv <- function(file, ...) {
  df <- utils::read.csv(file, header = TRUE,
                         stringsAsFactors = FALSE, ...)
  required <- c("tf", "target", "weight")
  missing <- setdiff(required, tolower(colnames(df)))
  if (length(missing) > 0) {
    stop("Missing required columns in weights CSV: ",
         paste(missing, collapse = ", "))
  }
  if (!"label" %in% tolower(colnames(df))) {
    warning("'label' column not found in weights CSV. Phenotype-specific weights will not be supported.")
  }
  df
}

#' Read SimiC-style AUC results from CSV (compatibility)
#'
#' @param file path to a CSV file containing AUC metrics.
#'   Expected columns include at least \code{cell}, \code{tf}, and \code{score}.
#' @param ... additional arguments passed to [utils::read.csv()].
#'
#' @return data.frame with AUC metrics.
#' @export
read_auc_csv <- function(file, ...) {
  df <- utils::read.csv(file, stringsAsFactors = FALSE, ...)
  required <- c("cell", "tf", "score")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop("Missing required columns in AUC CSV: ",
         paste(missing, collapse = ", "))
  }
  df
}

# -----WEIGHTS I/O pipeline ---------------------------------------------------

#' Convert a long-format weights data.frame to a list of TF × target matrices
#'
#' Used internally by \code{\link{load_SimiCviz_from_csv}} to transform
#' a data.frame with columns \code{tf}, \code{target}, \code{weight}
#' (and optionally \code{label}) into the named-list-of-matrices format
#' expected by \code{\link{SimiCvizExperiment}}.
#'
#' @param df A data.frame with at least columns \code{tf}, \code{target},
#'   \code{weight}. If a \code{label} column is present, one matrix is
#'   created per unique label; otherwise a single-element list is returned.
#' @return A named list of data.frames (TFs in rows, targets in columns).
#' @keywords internal
#' 
.convert_weights_df_to_list <- function(df) {

  # Split by label if present
  split_df <- if ("label" %in% tolower(colnames(df))) {
    split(df[, setdiff(names(df), "label")], df$label)
  } else {
    list(`0` = df)
  }

  result <- lapply(split_df, function(sub) {

    # Check duplicates within subset
    if (any(duplicated(sub[c("tf", "target")]))) {
      stop("Duplicate tf-target pairs found in weights data.frame.
           Please ensure each tf-target pair has a unique weight value.")
    }

    sub %>%
      tidyr::pivot_wider(
        names_from = target,
        values_from = weight,
        values_fill = list(weight = 0)
      ) %>%
      tibble::column_to_rownames("tf") %>%
      as.data.frame()
  })

  result
}

# ---- AUC I/O pipeline ---------------------------------------------------

#' Collect per-label AUC matrices using cell labels
#'
#' Given a list of AUC matrices (one per label) and a cell-labels data.frame,
#' subsets each matrix to include only the cells belonging to that label and collects them into a single data.frame. 
#' This is a helper function to transform raw per-label AUC outputs into the "collected" format used by \code{SimiCvizExperiment@auc}.
#'
#' @param auc_list Named list of AUC matrices (cells × TFs). Names should
#'   correspond to label identifiers (e.g. "0", "1").
#' @param cell_labels A data.frame with columns \code{cell} and \code{label},
#'   as returned by \code{\link{load_cell_labels}}.
#' @return A data.frame with cells in rows and TFs in columns, containing the activity scores for the specific labels.
#' @export

auc_list_to_df <- function(auc_list, cell_labels_df) {
  if (!is.list(auc_list) || length(auc_list) == 0L) {
    stop("`auc_list` must be a non-empty list of AUC matrices.")
  }
  if (!is.data.frame(cell_labels_df) ||
      !all(c("cell", "label") %in% colnames(cell_labels_df))) {
    stop("`cell_labels_df` must be a data.frame with columns 'cell' and 'label'.")
  }

  labels <- sort(unique(cell_labels_df$label))
  result <- list()

  for (lab in labels) {
    lab_key <- as.character(lab)
    cells_in_label <- cell_labels_df$cell[cell_labels_df$label == lab]

    # Find matching matrix in auc_list
    if (lab_key %in% names(auc_list)) {
      mat <- auc_list[[lab_key]]
    } else if (as.integer(lab) + 1L <= length(auc_list)) {
      # fallback: use positional index (0-based label → 1-based index) asumes labels are in order in auc_list
      mat <- auc_list[[as.integer(lab) + 1L]]
    } else {
      warning(sprintf("No AUC matrix found for label '%s'. Skipping.", lab_key))
      next
    }

    mat <- as.data.frame(mat)

    
    auc_cell_ids <- rownames(mat)
    if (all(cells_in_label %in% auc_cell_ids)) {
      
      # direct subsetting the auc matrix by cell labels for label i
      
      mat_sub <- mat[cells_in_label, , drop = FALSE]
    
    }  else if (length(intersect(cells_in_label, auc_cell_ids)) > 0L){

      warning(sprintf(
        "Label '%s': some cells in cell_labels_df are not found in AUC matrix rownames. Subsetting by intersection.",
        lab_key
      ))

      matched <- intersect(cells_in_label, auc_cell_ids)
      mat_sub <- mat[matched, , drop = FALSE]

    }  else if (length(intersect(cells_in_label, auc_cell_ids)) == 0L ) {
      # If cell_labels_df has only label info (i.e. loaded with .parse_cell_labels_vector -> cell_1, cell_2 names), match by position
      warning(sprintf(
        "Label '%s': cells in cell_labels_df are not found in AUC matrix rownames. Attempting positional matching.",
        lab_key
      ))
      
      cell_idx <- as.integer(gsub("cell_", "", cells_in_label))
      
      # Check that index is within bounds of the AUC matrix
      
      if (max(cell_idx) > nrow(mat)){
        stop(sprintf(
          "Label '%s': cell index derived from cell_labels_df exceeds AUC matrix rows. Cannot subset.",
          lab_key
        ))
      }
      
      mat_sub <- mat[cell_idx, , drop = FALSE]
      
      } else {
        stop("Label '%s': Double check input auc or cell_labels. Cannot subset.", lab_key)
      }


    result[[lab_key]] <- mat_sub
  }

  result

  # Output collected format
  auc_collected = dplyr::bind_rows(result)

  return(auc_collected)
}

# ---- Main pipeline loader -----------------------------------------------

#' Create a SimiCvizExperiment from SimiCPipeline output path
#'
#' Loads weights and AUC data following this priority chain for AUC:
#' \enumerate{
#'   \item Collected CSV (cells × TF, all labels merged)
#'         load_collected_auc
#'   \item Per-label pickle → collect → save CSV
#'         read_auc_pickle → auc_list_to_df → save_collected_auc
#' }
#'
#' The canonical AUC format stored in \code{SimiCvizExperiment@auc} is always
#' a single data.frame with cells in rows and TFs in columns (the "collected"
#' format), compatible with any GRN method producing per-cell TF activity scores.
#'
#' @param project_dir root directory of a SimiCPipeline project.
#' @param run_name name of the run experiment.
#' @param lambda1 L1 regularization parameter.
#' @param lambda2 L2 regularization parameter.
#' @param meta optional named list with metadata.
#'
#' @return \code{\link{SimiCvizExperiment}} object.
#' @export
load_SimiCPipeline <- function(project_dir,
                               run_name = NULL,
                               lambda1 = NULL,
                               lambda2 = NULL,
                               meta = list()) {
  if (!dir.exists(project_dir)) {
    stop("Project directory does not exist: ", project_dir)
  }
  if (is.null(run_name)) {
    stop("run_name must be provided to identify the experiment within the project directory.")
  }

  matrix_dir <- file.path(project_dir, "outputSimic/matrices", run_name)
  
  if (!dir.exists(matrix_dir)) {
    stop("Matrix directory does not exist for the specified project_dir and run_name: ", matrix_dir)
  }

  base_name <- paste0(run_name, "_L1_", lambda1, "_L2_", lambda2, "_")
  input_files  <- list.files(file.path(project_dir, "inputFiles"), full.names = TRUE) # For cell labels
  output_files <- list.files(matrix_dir, pattern = base_name, full.names = TRUE) # For weights and AUC

  # ---- Load weights ----
  # Work with filtered weights
  weights_file <- output_files[grepl("simic_matrices_filtered_BIC\\.pickle$", output_files)]

  if (length(weights_file) == 0L) {
    message("No filtered weights file found. Attempting to load unfiltered weights.")
    weights_file <- output_files[grepl("simic_matrices\\.pickle$", output_files)]
  }
  if (length(weights_file) > 1L) {
    print("Multiple weights files found:")
    print(weights_file)
    stop("Please review the directory outputs. Only one file per run_name + L1 + L2 should be found.")
  }
  if (length(weights_file) == 0L) {
    stop("No weights pickle file found for: ", base_name)
  }

  weights <- read_weights_pickle(weights_file)
  out     <- read_pickle(weights_file)
  r2_file <- out$adjusted_r_squared
  TF_ids        <- out$TF_ids
  query_targets <- out$query_targets

  meta$adjusted_r_squared <- r2_file
  meta$run_name <- run_name
  meta$lambda1  <- lambda1
  meta$lambda2  <- lambda2

  # ---- Load AUC (priority: collected CSV → pickle → NULL) ----
  auc_collected <- NULL

  # Step 1: Try collected CSV

  auc_csv_file <- output_files[grepl("*_collected\\.csv$", output_files)]

  if (length(auc_csv_file) == 1L) {
    message("Found collected AUC CSV file: ", auc_csv_file)
    auc_collected <- load_collected_auc(auc_csv_file)
  } else if (length(auc_csv_file) > 1L) {
    print("Multiple AUC CSV files found:")
    print(auc_csv_file)
    stop("Multiple AUC CSV files found. Please make sure only one collected.csv file by run_name and lambda parameters is found.")
  } else {
    # Step 2: Try per-label pickle
    message("No collected AUC CSV file found. Attempting to load per-label pickle.")
    auc_pickle_file <- output_files[grepl(".*wAUC.*\\.pickle$", output_files)]

    if (length(auc_pickle_file) == 0L) {
      message("No AUC pickle file found. AUC will be set to NULL.")
    } else if (length(auc_pickle_file) > 1L) {
      print("Multiple AUC pickle files found:")
      print(auc_pickle_file)
      stop("Multiple AUC pickle files found. Please make sure only one pickle file by run_name and lambda parameters is found.")
    } else {
      message("Found AUC pickle file: ", auc_pickle_file)
      auc_list_raw <- read_auc_pickle(auc_pickle_file)

      # Step 3: Subset and collect only if cell labels are available
      if (!is.null(cell_labels_df)) {
        # auc_list_subset <- subset_auc_by_labels(auc_list_raw, cell_labels_df)

        # Step 4: Collect into single cells × TF data.frame
        # auc_collected <- collect_auc(auc_list_subset, cell_labels = cell_labels_df)
        auc_collected  <- auc_list_to_df(auc_list_raw, cell_labels_df)
        # Step 5: Save collected CSV for future fast loading
        collected_csv_path <- file.path(
          matrix_dir,
          paste0(base_name, "wAUC_matrics_filtered_BIC_collected.csv")
        )
        tryCatch({
          save_collected_auc(auc_collected, collected_csv_path, overwrite = FALSE)
          message("Saved collected AUC CSV: ", collected_csv_path)
        }, error = function(e) {
          message("Could not save collected AUC CSV: ", conditionMessage(e))
        })
      } else {
        warning(
          "Cell labels are not available. Cannot produce a collected AUC ",
          "data.frame with unique cell identifiers. AUC will be stored as a ",
          "raw per-label list. Provide cell labels via `setCellLabels()` and ",
          # "then use subset_auc_by_labels`()`and `collect_auc()` to generate the collected format."
            "then use `auc_list_to_df()` to generate the collected format."
        )
        auc_collected <- NULL
      }
    }
  }
  
  # ---- Validate TFs match between weights and AUC ----
  if (!is.null(auc_collected)) {
    # Remove 'label' column for TF comparison if present
    tf_cols_auc <- setdiff(colnames(auc_collected), "label")
    if (!identical(sort(tf_cols_auc), sort(TF_ids))) {
      warning("TFs in AUC (", paste(tf_cols_auc, collapse = ", "),
              ") do not exactly match TFs in weights (",
              paste(TF_ids, collapse = ", "), ").")
    }
  }
  
  # ---- Load cell labels ----
  cell_labels_file <- input_files[grepl("annotation.csv$", input_files)]
  message("Looking for cell labels in input files...")
  
  cell_labels <- NULL
  cell_labels_df <- NULL
  if (length(cell_labels_file) == 0L) {
    message("No cell labels file found. Cell labels will be set to NULL.")
  } else if (length(cell_labels_file) > 1L) {
    print("Multiple cell labels files found:")
    print(cell_labels_file)
    stop("Multiple cell labels files found. Please ensure only one .txt file is present in the inputFiles directory.")
  } else {
    message("Found cell labels file: ", cell_labels_file)
    cell_labels_df <- load_cell_labels(cell_labels_file)
  }
  # Likely cell_labels_file was a vector and automatic cell names were generated cell_1, cell_2
  # if (cell_labels_df$cell[1]=="cell_1")
  # ---- Construct experiment ----
  # auc slot: collected df when available, otherwise raw per-label list from pickle
  if (!is.null(auc_collected)) {
    auc_slot <- list(collected = auc_collected)
  } else if (exists("auc_list_raw", inherits = FALSE) && !is.null(auc_list_raw)) {
    auc_slot <- list(per_label = lapply(auc_list_raw, as.data.frame))
  } else {
    auc_slot <- list()
  }

  simic <- SimiCvizExperiment(
    weights     = weights,
    auc         = auc_slot,
    cell_labels = cell_labels_df,
    meta        = meta
  )

  simic
}

# ---- CSV-based loader ----------------------------------------------------

#' Create a SimiCvizExperiment from CSV files
#'
#' @param weights_file CSV with weights (long format: tf, target, weight,
#'   and optionally label).
#' @param auc_file optional CSV with AUC metrics (collected cells × TF format).
#' @param meta optional named list with metadata.
#'
#' @return \code{\link{SimiCvizExperiment}} object.
#' @export
load_from_csv <- function(weights_file,
                          auc_file = NULL,
                          cell_labels = NULL,
                          meta = list()) {
  weight_df <- read_weights_csv(weights_file)
  weights   <- .convert_weights_df_to_list(weight_df)

  auc <- if (!is.null(auc_file)) {
    if (is.null(cell_labels)) {
      message("`cell_labels` are required if auc is provided to construct a SimiCvizExperiment.")
      stop("Please provide a valid `cell_labels` data.frame or file path.")
    }
    list(collected = read.csv(auc_file))
  } else {
    auc  <- NULL
  }


  SimiCvizExperiment(
    weights    = weights,
    auc        = auc,
    cell_labels = NULL,
    meta       = meta
  )
}

# ---- Export --------------------------------------------------------------

#' Save collected AUC data.frame to CSV
#'
#' Writes the cells × TF collected AUC data.frame to a CSV file,
#' preserving cell IDs as row names.
#'
#' @param df A data.frame as returned by \code{\link{auc_list_to_df}}.
#' @param file Output CSV file path.
#' @param overwrite Logical; overwrite if file exists (default FALSE).
#' @return Invisibly, the file path.
#' @export
save_collected_auc <- function(df, file, overwrite = FALSE) {
  if (!overwrite && file.exists(file)) {
    stop("File already exists and overwrite = FALSE: ", file)
  }
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(df, file, row.names = TRUE)
  message("Collected AUC saved to: ", file)
  invisible(file)
}

#' Export SimiCvizExperiment tables to CSV
#'
#' Exports weights and AUC tables into an organized directory structure:
#' \describe{
#'   \item{data/weights}{weights tables}
#'   \item{data/auc}{AUC tables (collected cells × TF format)}
#' }
#'
#' @param x \code{\link{SimiCvizExperiment}} object.
#' @param out_dir root output directory.
#' @param prefix filename prefix (default: "SimiCviz").
#' @param overwrite logical; overwrite existing files.
#'
#' @return Invisibly, a list of file paths.
#' @export
export_SimiCviz_csv <- function(x,
                                out_dir,
                                prefix = "SimiCviz",
                                overwrite = FALSE) {
  if (!is.SimiCvizExperiment(x)) {
    stop("x must be a SimiCvizExperiment")
  }
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  data_dir <- file.path(out_dir, "data")
  w_dir <- file.path(data_dir, "weights")
  a_dir <- file.path(data_dir, "auc")
  dir.create(w_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(a_dir, recursive = TRUE, showWarnings = FALSE)

  res <- list()

  if (length(x@weights) > 0L) {
    f_w <- file.path(w_dir, paste0(prefix, "_weights.csv"))
    if (!overwrite && file.exists(f_w)) {
      stop("File already exists and overwrite = FALSE: ", f_w)
    }
    utils::write.csv(x@weights, f_w, row.names = FALSE)
    res$weights <- f_w
  }

  if (length(x@auc) > 0L && !is.null(x@auc$collected)) {
    f_a <- file.path(a_dir, paste0(prefix, "_auc_collected.csv"))
    if (!overwrite && file.exists(f_a)) {
      stop("File already exists and overwrite = FALSE: ", f_a)
    }
    save_collected_auc(x@auc$collected, f_a, overwrite = overwrite)
    res$auc <- f_a
  }

  invisible(res)
}
