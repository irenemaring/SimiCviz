#' Plot adjusted R² distributions
#'
#' Plots histograms of adjusted R² values per label, similar to the Python
#' \code{SimiCVisualization$plot_r2_distribution}.
#'
#' @param adjusted_r_squared A \strong{named list} of numeric vectors, one per
#'   label.
#'   Names must be label identifiers (e.g. \code{"0"}, \code{"1"}, …).
#'   This metric is specific to SimiC's regression step; other methods may not
#'   produce it, so it is kept as an explicit input rather than extracted from
#'   the experiment object.
#' @param x Optional \code{SimiCvizExperiment} object used solely to resolve
#'   display names and colors for labels.
#'   If \code{NULL}, default \code{"Label i"} naming and a built-in palette are
#'   used.
#' @param labels Optional vector of labels to plot (subset of
#'   \code{names(adjusted_r_squared)}). Defaults to all.
#' @param threshold Numeric R² threshold line and summary-statistic cutoff
#'   (default \code{0.7}).
#' @param grid Optional numeric vector of length 2: \code{c(nrow, ncol)}.
#'   If \code{NULL}, one label per row is used.
#' @param save logical; save to PDF (default \code{FALSE}).
#' @param filename PDF filename (default \code{"R2_distributions.pdf"}).
#' @param out_dir output directory (default \code{getwd()}).
#' @param width,height page dimensions in inches (defaults 10 x 5*nrow).
#'
#' @return Called for side effects (plots). Returns \code{invisible(NULL)}.
#' @examples
#'   simic <- readRDS(system.file("extdata", "simic_full.rds", package = "SimiCviz"))
#'   plot_r2_distribution(simic@meta$adjusted_r_squared, simic, grid = c(2, 2))
#'   plot_r2_distribution(simic@meta$adjusted_r_squared, simic, threshold = 0.9, labels = c(0, 1))
#' @export
plot_r2_distribution <- function(adjusted_r_squared,
                                 x = NULL,
                                 labels = NULL,
                                 threshold = 0.7,
                                 grid = NULL,
                                 save = FALSE,
                                 filename = NULL,
                                 out_dir = getwd(),
                                 width = 10,
                                 height = NULL) {

  default_colors <- c("#5e82bd", "#ed9900", "#008b00", "#cd9b9b", "#800080", "#ff676f")

  # --- validate adjusted_r_squared ---
  if (missing(adjusted_r_squared) || is.null(adjusted_r_squared) || !is.list(adjusted_r_squared)) {
    stop("`adjusted_r_squared` must be a named list of numeric vectors (one per label).")
  }

  obj <- adjusted_r_squared

  all_labels <- names(obj)
  if (is.null(all_labels) || any(all_labels == "")) {
    stop("`adjusted_r_squared` must be a *named* list with label identifiers as names.")
  }

  # --- optional experiment object for display names / colors ---
  sim_obj <- NULL
  if (!is.null(x)) {
    if (!is.SimiCvizExperiment(x)) {
      stop("`x` must be a SimiCvizExperiment or NULL.")
    }
    sim_obj <- x
  }

  # --- label subsetting ---
  if (!is.null(labels)) {
    lab_chr <- as.character(labels)
    missing_labs <- setdiff(lab_chr, all_labels)
    if (length(missing_labs)) {
      warning("Some requested labels not found and will be ignored: ",
              paste(missing_labs, collapse = ", "))
    }
    use_labels <- intersect(lab_chr, all_labels)
    if (!length(use_labels)) {
      stop("None of the requested labels are present in adjusted_r_squared.")
    }
  } else {
    use_labels <- all_labels
  }

  n_labels <- length(use_labels)

  # --- color and name maps from object (same approach as plot_auc.R) ---
  lab_names_map <- if (!is.null(x) && length(x@label_names) > 0L) x@label_names else character()
  col_map       <- if (!is.null(x) && length(x@colors)      > 0L) x@colors      else character()

  # --- grid layout ---
  if (is.null(grid)){
    nrow <- n_labels
    ncol <- 1L
  } else if (length(grid) != 2 || !is.numeric(grid) || any(grid < 0)) {
    stop("`grid` must be a numeric positive vector of length 2 (nrow, ncol).")
  } else{
    # take the min number between the provided grid and the number of labels
    nrow <- ifelse(as.integer(grid[1]) == 0, n_labels, min(as.integer(grid[1]), n_labels)) 
    ncol <- ifelse(as.integer(grid[2]) == 0, ceiling(n_labels / nrow), min(as.integer(grid[2]), ceiling(n_labels / nrow))) 
  }

  if (is.null(height)) height <- 5 * nrow

  # ---- plotting helper (draws into current device) ----
  .draw_all <- function() {
    op <- par(no.readonly = TRUE)
    on.exit(par(op), add = TRUE)
    par(mfrow = c(nrow, ncol))

    for (i in seq_len(n_labels)) {
      lab_name   <- use_labels[[i]]
      r2_vals    <- obj[[lab_name]]

      if (is.null(r2_vals) || !length(r2_vals)) {
        plot.new(); title(main = paste("No R\u00B2 data for label", lab_name)); next
      }

      sel        <- r2_vals > threshold
      n_selected <- sum(sel, na.rm = TRUE)
      mean_r2    <- if (n_selected > 0) mean(r2_vals[sel], na.rm = TRUE) else 0

      lab_id_num <- tryCatch(as.integer(lab_name),
                             warning = function(.) NA_integer_,
                             error   = function(.) NA_integer_)
      if (is.na(lab_id_num)) lab_id_num <- i
      
      # Display name: use label_names slot if available, else "Label X"
      lab_display <- if (lab_name %in% names(lab_names_map))
        unname(lab_names_map[lab_name])
      else
        paste0("Label ", lab_name)

      # Color: use colors slot if available, else cycle default palette
      lab_color <- if (lab_name %in% names(col_map))
        unname(col_map[lab_name])
      else
        default_colors[(i - 1L) %% length(default_colors) + 1L]

      n_tfs <- if (!is.null(sim_obj)) max(100, length(sim_obj@tf_ids)) else 100L

      main_title <- sprintf(
        "%s\nTargets selected: %d, Mean R\u00b2: %.3f",
        lab_display, n_selected, mean_r2
      )

      r2_histogram(
        adjusted_r_squared = r2_vals,
        threshold          = threshold,
        n_tfs              = n_tfs,
        main               = main_title,
        col                = lab_color
      )
    }
  }

  if (save) {
    fname <- if (is.null(filename)) "R2_distributions.pdf" else filename
    fpath <- file.path(out_dir, fname)
    dir.create(dirname(fpath), recursive = TRUE, showWarnings = FALSE)
    grDevices::pdf(fpath, width = width, height = height, onefile = TRUE)
    .draw_all()
    grDevices::dev.off()
    message("Saved R\u00B2 distributions to: ", fpath)
  } else {
    .draw_all()
  }

  invisible(NULL)
}

# internal helper for a single histogram, similar to Python implementation
r2_histogram <- function(adjusted_r_squared,
                         threshold = 0.7,
                         n_tfs = 100,
                         main = "",
                         col = "grey") {
  hist(
    adjusted_r_squared,
    col      = col,
    breaks   = n_tfs,
    xlab     = "Adjusted R\u00B2",
    main     = main,
    border   = "black"
  )
  abline(v = threshold, col = "red", lwd = 2, lty = 2)
  grid(nx = NA, ny = NULL, col = "grey80", lty = "dotted")
}
