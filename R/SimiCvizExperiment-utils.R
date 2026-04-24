#' Calculate activity scores from SimiCvizExperiment weights
#'
#' Efficiently computes TF activity scores for all cells using label-specific
#' GRNs. Avoids redundant computation by calculating only for each cell's
#' corresponding label.
#'
#' @param simic SimiCvizExperiment object with weights and expression data
#' @param expression expression matrix (cells x genes). If NULL, uses data from
#'   metadata or file path if available.
#' @param adj_r2_threshold minimum adjusted R2 for target filtering (default: 0.7)
#' @param sort_by sorting criterion ("expression", "weight", or "adj_r2")
#' @param select_top_k keep only top K targets per TF (NULL = all)
#' @param percent_of_target percentage of targets to use (0-1)
#' @param n_cores number of cores for parallelization (1 = sequential)
#' @param backend parallelization backend ("sequential", "multicore", "multisession")
#' @param verbose print progress messages
#'
#' @return SimiCvizExperiment with computed AUC scores added to @auc slot
#' @export
calculate_activity_scores <- function(simic,
                                      expression = NULL,
                                      adj_r2_threshold = NULL,
                                      sort_by = "expression",
                                      select_top_k = NULL,
                                      percent_of_target = 1.0,
                                      n_cores = 1,
                                      backend = "sequential",
                                      verbose = TRUE) {
  
  if (!is.SimiCvizExperiment(simic)) {
    stop("simic must be a SimiCvizExperiment object")
  }
  
  # Validate weights structure
  if (length(simic@weights) == 0) {
    stop("simic object has no weights. Please provide weights to SimiCvizExperiment.")
  }
  qc_threshold <- NULL
  qc_type <- NULL
  adjusted_r_squared <- NULL
  adj_r2_list <- NULL
  
  # Load expression if not provided
  if (is.null(expression)) {
    stop("Expression matrix is required but not provided. Please provide an expression matrix (genes x cells) or path to expression matrix.")
  }
  # Set parameters to filter weights by adjusted R2 if available in metadata
  if (!is.null(adj_r2_threshold)){
    if (!is.null(simic@meta$adjusted_r_squared)) {
      qc_threshold <- adj_r2_threshold
      qc_type <- "adj_r2"
      adj_r2_list <- simic@meta$adjusted_r_squared
     }else{
      message("No adjusted R2 list found in simic metadata. Proceeding without adjusted R2 filtering.")
     }
  }
  # Initialize processor
  if (verbose) message("Initializing AUCProcessor...")
  processor <- AUCProcessor(
    weights = simic@weights,
    expression = expression,
    cell_labels = simic@cell_labels,
    qc_type = qc_type,
    qc_threshold = qc_threshold,
    adj_r2_list = adj_r2_list,
    n_cores = n_cores,
    backend = backend
  )
  
  # Compute activity scores
  if (verbose) message("Computing activity scores...")
  processor <- compute_auc(
    processor,
    sort_by = sort_by,
    select_top_k = select_top_k,
    percent_of_target = percent_of_target,
    verbose = verbose
  )
  
  # Extract results and add to simic object
  auc_matrix <- get_auc(processor, format = "wide")
  simic@auc$collected <- as.data.frame(auc_matrix)
  
  if (verbose) message("Activity scores computed and added to simic object")
  
  simic
}

# --- Helper operators ---

`%||%` <- function(x, y) if (is.null(x)) y else x
