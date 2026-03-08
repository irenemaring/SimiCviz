#' @keywords internal
.resolve_labels <- function(x, labels = NULL) {
  if (!is.null(labels)) {
    if (is.character(labels) && !is.null(x@label_names) && all(labels %in% x@label_names)) {
      labels <- names(x@label_names)[which(x@label_names %in% labels)]
    } else if (is.character(labels) && !all(labels %in% x@label_names)) {
      show <- labels[!labels %in% x@label_names]
      stop("Some labels not found in SimiCvizExperiment@label_names: \n\t missing -> ",
           paste(show, collapse = ", "))
    }
    return(as.integer(labels))
  }
  auc_df <- x@auc$collected
  if (!is.null(auc_df) && "label" %in% colnames(auc_df)) return(sort(unique(auc_df$label)))
  if (nrow(x@cell_labels) > 0L) return(sort(unique(x@cell_labels$label)))
  as.integer(names(x@weights))
}

#' @keywords internal
.resolve_tf_names <- function(x, tf_names = NULL) {
  all_tfs <- .auc_tf_names(x)
  if (is.null(tf_names)) return(sort(all_tfs))
  tf_names <- as.character(tf_names)
  valid <- intersect(tf_names, all_tfs)
  if (length(valid) == 0L) stop("None of the specified TFs found in AUC data.")
  valid
}

#' @keywords internal
.paginate_plots <- function(plot_list, per_page, nrow, ncol) {
  n     <- length(plot_list)
  pages <- list()
  for (start in seq(1, n, by = per_page)) {
    end   <- min(start + per_page - 1L, n)
    batch <- plot_list[start:end]
    while (length(batch) < nrow * ncol) batch[[length(batch) + 1L]] <- grid::nullGrob()
    pages[[length(pages) + 1L]] <- gridExtra::arrangeGrob(grobs = batch, nrow = nrow, ncol = ncol)
  }
  pages
}

#' RStudio-safe page drawing
#' @keywords internal
.draw_pages <- function(pages) {
  for (pg in pages) {
    gridExtra::grid.arrange(grobs   = pg$grobs,
                            nrow    = pg$nrow,
                            ncol    = pg$ncol,
                            widths  = pg$widths,
                            heights = pg$heights)
  }
}


#' RStudio-safe page saving to PDF
#' @keywords internal
.save_pages <- function(pages, fpath, width, height) {
  grDevices::pdf(fpath, width = width, height = height, onefile = TRUE)
  for (pg in pages) {
    gridExtra::grid.arrange(grobs   = pg$grobs,
                            nrow    = pg$nrow,
                            ncol    = pg$ncol,
                            widths  = pg$widths,
                            heights = pg$heights)
  }
  grDevices::dev.off()
}