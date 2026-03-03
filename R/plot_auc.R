# ---- Internal helpers ----------------------------------------------------

#' Subset collected AUC for a specific label
#' @keywords internal
.subset_auc_for_label <- function(x, label) {
  if (!is.SimiCvizExperiment(x)) stop("x must be a SimiCvizExperiment")
  lab_int <- as.integer(label)
  lab_chr <- as.character(lab_int)
  if (nrow(x@cell_labels) == 0L)
    stop(sprintf("No cell labels found in the experiment; cannot subset AUC for label '%s'.", lab_chr))
  auc_df <- tryCatch(x@auc$collected, error = function(e) NULL)
  if (!is.null(auc_df) && "label" %in% colnames(x@cell_labels)) {
    cell_idx <- x@cell_labels$cell[x@cell_labels$label == lab_int]
    return(auc_df[cell_idx, , drop = FALSE])
  }
  pl <- tryCatch(x@auc$per_label, error = function(e) NULL)
  if (!is.null(pl) && "label" %in% colnames(x@cell_labels)) {
    if (lab_chr %in% names(pl)) {
      auc_tmp  <- as.data.frame(pl[[lab_chr]])
      cell_idx <- x@cell_labels$cell[x@cell_labels$label == lab_int]
      return(auc_tmp[cell_idx, , drop = FALSE])
    }
  }
  stop(sprintf("No AUC data found for label '%s'.", lab_chr))
}

#' @keywords internal
.auc_tf_names <- function(x) {
  auc_df <- x@auc$collected
  if (!is.null(auc_df)) return(setdiff(colnames(auc_df), "label"))
  pl <- x@auc$per_label
  if (!is.null(pl) && length(pl) > 0L) return(unique(unlist(lapply(pl, colnames))))
  character()
}

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

# ---- ECDF metrics -------------------------------------------------------

#' @keywords internal
.compute_ecdf_metrics <- function(values, x_lower = 0, x_upper = 1, percentile = 0.5) {
  n           <- length(values)
  sorted_vals <- sort(values)
  in_range    <- sorted_vals[sorted_vals > x_lower & sorted_vals < x_upper]
  breakpoints <- c(x_lower, in_range, x_upper)

  total_area <- 0
  area_at_p  <- NULL
  x_at_p     <- NULL

  for (i in seq_len(length(breakpoints) - 1L)) {
    x_left  <- breakpoints[i]
    x_right <- breakpoints[i + 1L]
    width   <- x_right - x_left
    ecdf_left  <- sum(sorted_vals <= x_left)  / n
    ecdf_right <- sum(sorted_vals <= x_right) / n
    if (is.null(area_at_p) && ecdf_right >= percentile) {
      if (ecdf_left >= percentile) {
        x_at_p    <- x_left
        area_at_p <- total_area
      } else {
        total_area <- total_area + ecdf_left * width
        x_at_p     <- x_right
        area_at_p  <- total_area
        next
      }
    }
    total_area <- total_area + ecdf_left * width
  }
  if (is.null(area_at_p)) { x_at_p <- x_upper; area_at_p <- total_area }
  range_width <- x_upper - x_lower
  if (range_width > 0) total_area <- total_area / range_width
  list(ecdf_auc = total_area, auc_at_percentile = area_at_p, x_at_percentile = x_at_p)
}

# ---- Public: calculate_ecdf_auc -----------------------------------------

#' @export
calculate_ecdf_auc <- function(x, tf_names = NULL, labels = NULL,
                               integration_range = c(0, 1), percentile = 0.5) {
  if (!is.SimiCvizExperiment(x)) stop("x must be a SimiCvizExperiment")
  labels   <- .resolve_labels(x, labels)
  tf_names <- .resolve_tf_names(x, tf_names)
  x_lower  <- integration_range[1]
  x_upper  <- integration_range[2]

  rows <- list()
  for (tf in tf_names) {
    row <- list(TF = tf)
    for (lab in labels) {
      lab_name <- x@label_names[as.character(lab)]
      sub  <- tryCatch(.subset_auc_for_label(x, lab), error = function(e) NULL)
      if (is.null(sub) || !tf %in% colnames(sub)) {
        row[[paste0(lab_name, "_ecdf_auc")]] <- NA_real_
        row[[paste0(lab_name, "_auc50")]]    <- NA_real_
        row[[paste0(lab_name, "_x_at_p50")]] <- NA_real_
        next
      }
      vals <- stats::na.omit(sub[[tf]])
      if (length(vals) < 1L) {
        row[[paste0(lab_name, "_ecdf_auc")]] <- NA_real_
        row[[paste0(lab_name, "_auc50")]]    <- NA_real_
        row[[paste0(lab_name, "_x_at_p50")]] <- NA_real_
        next
      }
      m <- .compute_ecdf_metrics(vals, x_lower, x_upper, percentile)
      row[[paste0(lab_name, "_ecdf_auc")]] <- m$ecdf_auc
      row[[paste0(lab_name, "_auc50")]]    <- m$auc_at_percentile
      row[[paste0(lab_name, "_x_at_p50")]] <- m$x_at_percentile
    }
    rows[[length(rows) + 1L]] <- row
  }

  df <- do.call(rbind, lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
  rownames(df) <- df$TF
  df$TF <- NULL

  if (length(labels) > 1L) {
    lab_names <- x@label_names
    for (suffix in c("_ecdf_auc", "_auc50", "_x_at_p50")) {
      cols <- intersect(paste0(lab_names, suffix), colnames(df))
      if (length(cols) > 1L) {
        vals_mat <- df[, cols, drop = FALSE]
        df[[paste0("delta", suffix)]] <- apply(vals_mat, 1, max, na.rm = TRUE) -
                                         apply(vals_mat, 1, min, na.rm = TRUE)
      }
    }
  }
  df
}

# ---- Public: plot_auc_distributions --------------------------------------

#' @export
plot_auc_distributions <- function(x, tf_names = NULL, labels = NULL,
                                   fill = TRUE, alpha = 0.5, bw_adjust = 1,
                                   rug = FALSE, grid = c(4L, 2L), save = FALSE,
                                   filename = NULL, out_dir = getwd(),
                                   width = 14, height = NULL) {
  if (!requireNamespace("ggplot2",   quietly = TRUE)) stop("ggplot2 required")
  if (!requireNamespace("gridExtra", quietly = TRUE)) stop("gridExtra required")

  labels        <- .resolve_labels(x, labels)
  lab_names     <- x@label_names
  tf_names      <- .resolve_tf_names(x, tf_names)
  n_tfs         <- length(tf_names)
  col_map       <- x@colors
  col_map_named <- stats::setNames(unname(col_map), lab_names[names(col_map)])
  diss_score    <- calculate_dissimilarity(x, tf_names = tf_names, labels = labels, verbose = FALSE)

  if (is.null(grid)) {
    grid_cols      <- 2L
    grid_rows      <- ceiling(n_tfs / grid_cols)
    plots_per_page <- n_tfs
  } else {
    grid_rows      <- as.integer(grid[1])
    grid_cols      <- as.integer(grid[2])
    plots_per_page <- grid_rows * grid_cols
  }
  if (is.null(height)) height <- 5 * grid_rows

  plot_list <- list()
  for (tf in tf_names) {
    long <- .build_density_long(x, tf, labels)
    if (nrow(long) == 0L) {
      p <- ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = paste("No data for", tf), size = 5) +
        ggplot2::theme_void() + ggplot2::ggtitle(tf)
      plot_list[[length(plot_list) + 1L]] <- p
      next
    }
    label_counts <- table(long$label_name)
    new_labels   <- paste0(names(label_counts), " (n=", label_counts, ")")

    p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$value,
                                             fill = .data$label_name,
                                             colour = .data$label_name))
    if (fill) {
      p <- p + ggplot2::geom_density(alpha = alpha, adjust = bw_adjust)
    } else {
      p <- p + ggplot2::geom_density(alpha = 0, adjust = bw_adjust, linewidth = 1)
    }
    if (rug) p <- p + ggplot2::geom_rug(alpha = 0.4, colour = "grey40")

    p <- p +
      ggplot2::scale_fill_manual(values = col_map_named, labels = new_labels) +
      ggplot2::scale_colour_manual(values = col_map_named, labels = new_labels) +
      ggplot2::labs(x = "Activity Score", y = "Density", fill = NULL, colour = NULL) +
      ggplot2::xlim(0, 1) +
      ggplot2::ggtitle(tf, subtitle = paste("Dissimilarity score:",
                                             signif(diss_score[tf, "MinMax_score"], 3))) +
      ggplot2::theme_classic() +
      ggplot2::theme(legend.position = c(1, 1), legend.justification = c(1, 1),
                     legend.direction = "vertical",
                     plot.title = ggplot2::element_text(face = "bold", size = 12))
    plot_list[[length(plot_list) + 1L]] <- p
  }

  pages <- .paginate_plots(plot_list, plots_per_page, grid_rows, grid_cols)

  if (save) {
    fname <- filename %||% "AUC_distributions.pdf"
    fpath <- file.path(out_dir, fname)
    dir.create(dirname(fpath), recursive = TRUE, showWarnings = FALSE)
    .save_pages(pages, fpath, width, height)
    message("Saved AUC distributions to: ", fpath)
  } else {
    .draw_pages(pages)
  }
  invisible(pages)
}

# ---- Internal: ECDF table grob ------------------------------------------

#' @keywords internal
.build_ecdf_table_grob <- function(ecdf_df, tf, labels, label_names, col_map) {
  if (!requireNamespace("gridExtra", quietly = TRUE)) return(NULL)
  if (is.null(ecdf_df) || !tf %in% rownames(ecdf_df)) return(NULL)

  lab_names  <- unname(label_names[as.character(labels)])
  ecdf_vals  <- vapply(lab_names, function(ln) {
    col <- paste0(ln, "_ecdf_auc"); if (col %in% colnames(ecdf_df)) ecdf_df[tf, col] else NA_real_
  }, numeric(1))
  auc50_vals <- vapply(lab_names, function(ln) {
    col <- paste0(ln, "_auc50");    if (col %in% colnames(ecdf_df)) ecdf_df[tf, col] else NA_real_
  }, numeric(1))
  median_vals <- vapply(lab_names, function(ln) {
    col <- paste0(ln, "_x_at_p50"); if (col %in% colnames(ecdf_df)) ecdf_df[tf, col] else NA_real_
  }, numeric(1))

  fmt    <- function(v) ifelse(is.na(v), "\u2014", sprintf("%.4f", v))
  tab_df <- data.frame(Metric = c("ECDF-AUC", "AUC50", "Median AS"), stringsAsFactors = FALSE)
  for (i in seq_along(lab_names))
    tab_df[[lab_names[i]]] <- c(fmt(ecdf_vals[i]), fmt(auc50_vals[i]), fmt(median_vals[i]))

  header_col_fills <- vapply(as.character(labels), function(cid)
    if (cid %in% names(col_map)) col_map[[cid]] else "grey90", character(1), USE.NAMES = FALSE)

  tt <- gridExtra::ttheme_minimal(
    core    = list(fg_params = list(fontsize = 8, hjust = 0.5, x = 0.5),
                   bg_params = list(fill = "white",         col = "grey80")),
    colhead = list(fg_params = list(fontsize = 8, fontface = "bold", hjust = 0.5, x = 0.5),
                   bg_params = list(fill = header_col_fills, col = "grey80")),
    rowhead = list(fg_params = list(fontsize = 8, fontface = "bold", hjust = 0.5, x = 0.5),
                   bg_params = list(fill = "grey85",         col = "grey80"))
  )
  gridExtra::tableGrob(tab_df[, -1, drop = FALSE], rows = tab_df$Metric, theme = tt)
}

# ---- Public: plot_auc_cumulative -----------------------------------------

#' @export
plot_auc_cumulative <- function(x, tf_names = NULL, labels = NULL, alpha = 0.8,
                                rug = FALSE, grid = c(4L, 2L), percentile = 0.5,
                                include_table = TRUE, save = FALSE, filename = NULL,
                                out_dir = getwd(), width = 14, height = NULL) {
  if (!requireNamespace("ggplot2",   quietly = TRUE)) stop("ggplot2 required")
  if (!requireNamespace("gridExtra", quietly = TRUE)) stop("gridExtra required")

  labels        <- .resolve_labels(x, labels)
  tf_names      <- .resolve_tf_names(x, tf_names)
  n_tfs         <- length(tf_names)
  lab_names     <- x@label_names
  col_map       <- x@colors
  col_map_named <- stats::setNames(unname(col_map), lab_names[names(col_map)])

  if (is.null(grid)) {
    grid_cols      <- 2L
    grid_rows      <- ceiling(n_tfs / grid_cols)
    plots_per_page <- n_tfs
  } else {
    grid_rows      <- as.integer(grid[1])
    grid_cols      <- as.integer(grid[2])
    plots_per_page <- grid_rows * grid_cols
  }
  if (is.null(height)) height <- ifelse(include_table, 6.5, 5) * grid_rows

  ecdf_df <- tryCatch(
    calculate_ecdf_auc(x, tf_names = tf_names, labels = labels, percentile = percentile),
    error = function(e) NULL
  )

  plot_list <- list()
  for (tf in tf_names) {
    long <- .build_density_long(x, tf, labels)
    if (nrow(long) == 0L) {
      p <- ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5, label = paste("No data for", tf), size = 5) +
        ggplot2::theme_void() + ggplot2::ggtitle(tf)
      plot_list[[length(plot_list) + 1L]] <- p
      next
    }
    label_counts <- table(long$label_name)
    new_labels   <- paste0(names(label_counts), " (n=", label_counts, ")")

    p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$value, colour = .data$label_name)) +
      ggplot2::stat_ecdf(geom = "step", linewidth = 1, alpha = alpha) +
      ggplot2::scale_colour_manual(values = col_map_named, labels = new_labels) +
      ggplot2::labs(x = "Activity Score", y = "Cumulative Probability", colour = NULL) +
      ggplot2::xlim(0, 1) +
      ggplot2::theme_classic() +
      ggplot2::theme(legend.position = c(1, 0), legend.justification = c(1, 0),
                     legend.direction = "vertical",
                     plot.title = ggplot2::element_text(face = "bold", size = 12))

    if (rug) p <- p + ggplot2::geom_rug(alpha = 0.4, colour = "grey40", sides = "b")

    subtitle_parts <- NULL
    if (!is.null(ecdf_df) && tf %in% rownames(ecdf_df)) {
      delta_cols <- grep("^delta_", colnames(ecdf_df), value = TRUE)
      sub_parts  <- character()
      if ("delta_ecdf_auc" %in% delta_cols)
        sub_parts <- c(sub_parts, sprintf("Delta*AUC: %.4f",  ecdf_df[tf, "delta_ecdf_auc"]))
      if ("delta_auc50"    %in% delta_cols)
        sub_parts <- c(sub_parts, sprintf("Delta*AUC50: %.4f", ecdf_df[tf, "delta_auc50"]))
      if (length(sub_parts) > 0L) subtitle_parts <- paste(sub_parts, collapse = "  |  ")
    }
    p <- p + ggplot2::ggtitle(tf, subtitle_parts)

    if (include_table) {
      tbl_grob <- .build_ecdf_table_grob(ecdf_df, tf, labels, lab_names, col_map)
      if (!is.null(tbl_grob)) {
        combined <- gridExtra::arrangeGrob(p, tbl_grob,
                                           nrow = 2, heights = grid::unit(c(3, 1), "null"))
        plot_list[[length(plot_list) + 1L]] <- combined
      } else {
        plot_list[[length(plot_list) + 1L]] <- p
      }
    } else {
      plot_list[[length(plot_list) + 1L]] <- p
    }
  }

  pages <- .paginate_plots(plot_list, plots_per_page, grid_rows, grid_cols)

  if (save) {
    fname <- filename %||% "AUC_cumulative.pdf"
    fpath <- file.path(out_dir, fname)
    dir.create(dirname(fpath), recursive = TRUE, showWarnings = FALSE)
    .save_pages(pages, fpath, width, height)
    message("Saved AUC cumulative to: ", fpath)
  } else {
    .draw_pages(pages)
  }
  invisible(pages)
}

# ---- Public: plot_auc_summary_statistics ---------------------------------

#' @export
plot_auc_summary_statistics <- function(x, labels = NULL, high_threshold = 0.5,
                                        save = FALSE, filename = NULL,
                                        out_dir = getwd(), width = 14, height = 10) {
  if (!requireNamespace("ggplot2",   quietly = TRUE)) stop("ggplot2 required")
  if (!requireNamespace("gridExtra", quietly = TRUE)) stop("gridExtra required")

  labels <- .resolve_labels(x, labels)

  all_long <- do.call(rbind, lapply(labels, function(lab) {
    sub  <- tryCatch(.subset_auc_for_label(x, lab), error = function(e) NULL)
    if (is.null(sub) || ncol(sub) == 0L) return(NULL)
    vals <- unlist(sub, use.names = FALSE)
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0L) return(NULL)
    data.frame(value = vals,
               label_name = unname(x@label_names[as.character(lab)]),
               stringsAsFactors = FALSE)
  }))
  if (is.null(all_long) || nrow(all_long) == 0L)
    stop("No AUC data available for the requested labels.")

  lab_order             <- unname(x@label_names)
  all_long$label_name   <- factor(all_long$label_name, levels = lab_order)
  col_map               <- stats::setNames(unname(x@colors), unname(x@label_names))

  p_box <- ggplot2::ggplot(all_long, ggplot2::aes(x = .data$label_name, y = .data$value,
                                                    fill = .data$label_name)) +
    ggplot2::geom_boxplot(alpha = 0.6, outlier.size = 0.5) +
    ggplot2::scale_fill_manual(values = col_map) +
    ggplot2::labs(x = NULL, y = "Activity Score", title = "Activity Score Distribution (Boxplot)") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   plot.title  = ggplot2::element_text(face = "bold"))

  p_violin <- ggplot2::ggplot(all_long, ggplot2::aes(x = .data$label_name, y = .data$value,
                                                       fill = .data$label_name)) +
    ggplot2::geom_violin(alpha = 0.6, draw_quantiles = c(0.25, 0.5, 0.75)) +
    ggplot2::scale_fill_manual(values = col_map) +
    ggplot2::labs(x = NULL, y = "Activity Score", title = "Activity Score Distribution (Violin)") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   plot.title  = ggplot2::element_text(face = "bold"))

  summary_df <- do.call(rbind, lapply(labels, function(lab) {
    sub  <- .subset_auc_for_label(x, lab)
    vals <- unlist(sub, use.names = FALSE); vals <- vals[!is.na(vals)]
    data.frame(label_name = unname(x@label_names[as.character(lab)]),
               mean_val = mean(vals), sd_val = stats::sd(vals), stringsAsFactors = FALSE)
  }))
  summary_df$label_name <- factor(summary_df$label_name, levels = lab_order)

  p_mean <- ggplot2::ggplot(summary_df, ggplot2::aes(x = .data$label_name, y = .data$mean_val,
                                                       fill = .data$label_name)) +
    ggplot2::geom_col(alpha = 0.7, colour = "black", linewidth = 0.5) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$mean_val - .data$sd_val,
                                         ymax = .data$mean_val + .data$sd_val), width = 0.3) +
    ggplot2::scale_fill_manual(values = col_map) +
    ggplot2::labs(x = NULL, y = "Mean Activity Score", title = "Mean Activity Score by Phenotype") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   plot.title  = ggplot2::element_text(face = "bold"))

  ha_df <- do.call(rbind, lapply(labels, function(lab) {
    sub          <- .subset_auc_for_label(x, lab)
    mean_per_tf  <- colMeans(sub, na.rm = TRUE)
    data.frame(label_name = unname(x@label_names[as.character(lab)]),
               n_high = sum(mean_per_tf > high_threshold, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
  ha_df$label_name <- factor(ha_df$label_name, levels = lab_order)

  p_high <- ggplot2::ggplot(ha_df, ggplot2::aes(x = .data$label_name, y = .data$n_high,
                                                  fill = .data$label_name)) +
    ggplot2::geom_col(alpha = 0.7, colour = "black", linewidth = 0.5) +
    ggplot2::scale_fill_manual(values = col_map) +
    ggplot2::labs(x = NULL, y = "Number of TFs",
                  title = sprintf("TFs with High Activity (Mean AS > %.1f)", high_threshold)) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   plot.title  = ggplot2::element_text(face = "bold"))

  combined <- gridExtra::arrangeGrob(p_box, p_violin, p_mean, p_high, ncol = 2, nrow = 2,
                                      top = grid::textGrob("Activity Scores Summary Statistics",
                                        gp = grid::gpar(fontsize = 16, fontface = "bold")))
  if (save) {
    fname <- filename %||% "AUC_summary_statistics.pdf"
    fpath <- file.path(out_dir, fname)
    dir.create(dirname(fpath), recursive = TRUE, showWarnings = FALSE)
    .save_pages(list(combined), fpath, width, height)
    message("Saved AUC summary statistics to: ", fpath)
  } else {
    .draw_pages(list(combined))
  }
  invisible(combined)
}

# ---- Public: plot_auc_heatmap -------------------------------------------

#' @export
plot_auc_heatmap <- function(x, tf_names = NULL, labels = NULL, top_n = NULL,
                             save = FALSE, filename = NULL, out_dir = getwd(),
                             width = 10, height = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")

  labels   <- .resolve_labels(x, labels)
  tf_names <- .resolve_tf_names(x, tf_names)

  rows <- list()
  for (tf in tf_names) {
    row <- list(tf = tf)
    for (lab in labels) {
      sub <- tryCatch(.subset_auc_for_label(x, lab), error = function(e) NULL)
      val <- if (!is.null(sub) && tf %in% colnames(sub)) mean(sub[[tf]], na.rm = TRUE) else NA_real_
      row[[x@label_names[as.character(lab)]]] <- val
    }
    rows[[length(rows) + 1L]] <- as.data.frame(row, stringsAsFactors = FALSE)
  }
  mat_df       <- do.call(rbind, rows)
  rownames(mat_df) <- mat_df$tf
  num_mat      <- mat_df[, setdiff(colnames(mat_df), "tf"), drop = FALSE]

  if (!is.null(top_n) && top_n < nrow(num_mat)) {
    rng  <- apply(num_mat, 1, function(r) diff(range(r, na.rm = TRUE)))
    keep <- names(sort(rng, decreasing = TRUE))[seq_len(top_n)]
    num_mat <- num_mat[keep, , drop = FALSE]
  }
  if (is.null(height)) height <- max(6, nrow(num_mat) * 0.35 + 2)

  long <- data.frame(
    tf        = rep(rownames(num_mat), times = ncol(num_mat)),
    condition = rep(colnames(num_mat), each  = nrow(num_mat)),
    auc       = unlist(num_mat, use.names = FALSE),
    stringsAsFactors = FALSE
  )
  long$tf        <- factor(long$tf,        levels = rev(rownames(num_mat)))
  long$condition <- factor(long$condition, levels = colnames(num_mat))

  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$condition, y = .data$tf, fill = .data$auc)) +
    ggplot2::geom_tile(colour = "grey80") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", .data$auc)), size = 3) +
    ggplot2::scale_fill_viridis_c(option = "C", na.value = "grey90") +
    ggplot2::labs(x = "Phenotype", y = "TF", fill = "Mean AUC", title = "Mean Activity Score per TF") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   plot.title  = ggplot2::element_text(face = "bold", size = 14))

  if (save) {
    fname <- filename %||% "AUC_heatmap.pdf"
    fpath <- file.path(out_dir, fname)
    dir.create(dirname(fpath), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(fpath, p, width = width, height = height)
    message("Saved AUC heatmap to: ", fpath)
  } else {
    print(p)
  }
  invisible(p)
}

# ---- Public: export_auc_pdfs --------------------------------------------

#' @export
export_auc_pdfs <- function(x, out_dir, prefix = "SimiCviz", labels = NULL,
                            overwrite = FALSE, ...) {
  plot_dir <- file.path(out_dir, "plots", "auc")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  res <- list()
  tryCatch({ plot_auc_distributions(x, labels = labels, save = TRUE,
    filename = paste0(prefix, "_auc_distributions.pdf"), out_dir = plot_dir, ...)
    res$distributions <- file.path(plot_dir, paste0(prefix, "_auc_distributions.pdf"))
  }, error = function(e) message("Skipping distributions: ", conditionMessage(e)))
  tryCatch({ plot_auc_cumulative(x, labels = labels, save = TRUE,
    filename = paste0(prefix, "_auc_cumulative.pdf"), out_dir = plot_dir, ...)
    res$cumulative <- file.path(plot_dir, paste0(prefix, "_auc_cumulative.pdf"))
  }, error = function(e) message("Skipping cumulative: ", conditionMessage(e)))
  tryCatch({ plot_auc_summary_statistics(x, labels = labels, save = TRUE,
    filename = paste0(prefix, "_auc_summary.pdf"), out_dir = plot_dir, ...)
    res$summary <- file.path(plot_dir, paste0(prefix, "_auc_summary.pdf"))
  }, error = function(e) message("Skipping summary: ", conditionMessage(e)))
  tryCatch({ plot_auc_heatmap(x, labels = labels, save = TRUE,
    filename = paste0(prefix, "_auc_heatmap.pdf"), out_dir = plot_dir, ...)
    res$heatmap <- file.path(plot_dir, paste0(prefix, "_auc_heatmap.pdf"))
  }, error = function(e) message("Skipping heatmap: ", conditionMessage(e)))
  invisible(res)
}

# ---- Internal helpers ---------------------------------------------------

#' @keywords internal
.build_density_long <- function(x, tf, labels) {
  rows <- list()
  for (lab in labels) {
    sub  <- tryCatch(.subset_auc_for_label(x, lab), error = function(e) NULL)
    if (is.null(sub) || !tf %in% colnames(sub)) next
    vals <- stats::na.omit(sub[[tf]])
    if (length(vals) < 1L) next
    rows[[length(rows) + 1L]] <- data.frame(
      value      = vals,
      label      = lab,
      label_name = unname(x@label_names[as.character(lab)]),
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0L)
    return(data.frame(value = numeric(), label = integer(), label_name = character(),
                      stringsAsFactors = FALSE))
  long <- do.call(rbind, rows)
  long$label_name <- factor(long$label_name,
                             levels = intersect(unname(x@label_names), unique(long$label_name)))
  long
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



#' @keywords internal
`%||%` <- function(a, b) if (!is.null(a)) a else b
