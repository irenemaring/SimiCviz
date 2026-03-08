# ---- Internal helpers ----------------------------------------------------

#' Convert weight list to long-format data frame
#'
#' @param x SimiCvizExperiment
#' @return data.frame with columns: tf, target, weight, label, label_name
#' @keywords internal
.weights_to_long <- function(x) {
  if (length(x@weights) == 0L) stop("No weight matrices found in the experiment.")
  W_long_list <- list()
  for (lab in names(x@weights)) {
    W <- as.matrix(x@weights[[lab]])
    W_long <- as.data.frame(as.table(W))
    colnames(W_long) <- c("tf", "target", "weight")
    W_long$label <- as.integer(lab)
    W_long$label_name <- unname(x@label_names[lab])
    W_long$tf <- as.character(W_long$tf)
    W_long$target <- as.character(W_long$target)
    W_long_list[[lab]] <- W_long
  }
  dplyr::bind_rows(W_long_list)
}

#' Render a single weight barplot on a ggplot
#'
#' @param df long-format data.frame (already subsetted) with columns:
#'   gene, weight, label_name
#' @param gene_name character
#' @param gene_type character; "tf" or "target"
#' @param col_map named character vector of colours (keyed by label_name)
#' @param top_n integer or NULL; keep only top N partners by mean abs weight
#' @param allowed_genes character vector or NULL; restrict partner genes to this set
#' @return ggplot object
#' @keywords internal
.make_weight_barplot <- function(df, gene_name, gene_type,
                                 col_map, top_n = NULL,
                                 allowed_genes = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")

  x_label <- if (gene_type == "tf") "Target Genes" else "Transcription Factors"
  title    <- if (gene_type == "tf") paste("TF:", gene_name) else paste("Target:", gene_name)

  if (nrow(df) == 0L) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5,
                          label = paste("No non-zero weights\nfor", gene_name),
                          size = 4, hjust = 0.5, vjust = 0.5) +
        ggplot2::theme_void() +
        ggplot2::ggtitle(title)
    )
  }

  # The partner gene column (x-axis)
  partner_col <- if (gene_type == "tf") "target" else "tf"

  # Restrict to allowed partner genes if provided
  if (!is.null(allowed_genes)) {
    df_lab_list <- list()
    for (lab in unique(df$label)) {
        lab_key <- as.character(lab)
        df_lab <- df[df$label == lab, , drop = FALSE]
        df_lab <- df_lab[df_lab[[partner_col]] %in% allowed_genes[[lab_key]], , drop = FALSE]
        df_lab_list[[lab_key]] <- df_lab
       }
    df <- dplyr::bind_rows(df_lab_list)

  if (nrow(df) == 0L) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.5,
                          label = paste("No weights after filtering\nfor", gene_name),
                          size = 4, hjust = 0.5, vjust = 0.5) +
        ggplot2::theme_void() +
        ggplot2::ggtitle(title)
    )
  }
  }

  # Rename partner column to "gene" for uniform downstream code
  df$gene <- df[[partner_col]]

  # Optionally keep top N partners by mean absolute weight
  if (!is.null(top_n)) {
    mean_abs <- tapply(abs(df$weight), df$gene, mean)
    top_genes <- names(sort(mean_abs, decreasing = TRUE))[
      seq_len(min(top_n, length(mean_abs)))
    ]
    df <- df[df$gene %in% top_genes, , drop = FALSE]
  }

  # Order partner genes by max absolute weight (descending)
  max_abs   <- tapply(abs(df$weight), df$gene, max)
  gene_order <- names(sort(max_abs, decreasing = TRUE))

  # Complete the grid: every gene x label_name combination must exist 
  all_label_names <- names(col_map)
  full_grid <- expand.grid(gene = gene_order, label_name = all_label_names,
                           stringsAsFactors = FALSE)
  df_slim <- df[, c("gene", "label_name", "weight"), drop = FALSE]
  df_plot <- merge(full_grid, df_slim, by = c("gene", "label_name"), all.x = TRUE)
  df_plot$weight[is.na(df_plot$weight)] <- 0

  df_plot$gene       <- factor(df_plot$gene,       levels = gene_order)
  df_plot$label_name <- factor(df_plot$label_name, levels = all_label_names)

    # In target mode we do not want the weights that are 0 for all phenotypes, so we filter them out
  if (gene_type == "target") {
    gene_zero <- tapply(df_plot$weight, df_plot$gene, function(w) all(w == 0))
    genes_to_keep <- names(gene_zero)[!gene_zero]
    df_plot <- df_plot[df_plot$gene %in% genes_to_keep, , drop = FALSE]
    }
  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = .data$gene,
                                              y = .data$weight,
                                              fill = .data$label_name)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8),
                      colour = "black", linewidth = 0.3, width = 0.8) +
    ggplot2::scale_fill_manual(values = col_map, drop = FALSE) +
    ggplot2::geom_hline(yintercept = 0, colour = "black", linewidth = 0.5) +
    ggplot2::labs(x = x_label, y = "Weight", fill = NULL, title = title) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1,
                                           size = if (gene_type == "tf") 7 else 8),
      legend.position.inside = c(0.5, 0.5),
      legend.justification   = c(0.5, 0.5),
      legend.direction       = "vertical",
      plot.title = ggplot2::element_text(face = "bold", size = 12)
    )

  p
}


# ---- Shared paginated engine --------------------------------------------

#' Shared paginated weight barplot engine
#'
#' @param x SimiCvizExperiment
#' @param gene_names character vector or NULL (all genes of that type)
#' @param gene_type "tf" or "target"
#' @param labels integer vector of labels (default: all)
#' @param top_n integer or NULL; max partners shown per gene (TF mode only)
#' @param allowed_targets character vector or NULL; restrict targets to this set
#'   (TF mode only; use to pass pre-filtered target lists)
#' @param grid integer vector c(nrow, ncol) per page; NULL = single page
#' @param save logical
#' @param filename PDF filename
#' @param out_dir output directory
#' @param width,height page dimensions in inches
#' @return Invisibly, list of page grobs.
#' @keywords internal
.plot_weights_paginated <- function(x,
                                    gene_names      = NULL,
                                    gene_type       = c("tf", "target"),
                                    labels          = NULL,
                                    top_n           = NULL,
                                    allowed_targets = NULL,
                                    grid            = c(4L, 1L),
                                    save            = FALSE,
                                    filename        = NULL,
                                    out_dir         = getwd(),
                                    width           = 16,
                                    height          = NULL) {
  if (!requireNamespace("ggplot2",   quietly = TRUE)) stop("ggplot2 required")
  if (!requireNamespace("gridExtra", quietly = TRUE)) stop("gridExtra required")

  gene_type <- match.arg(gene_type)
  labels    <- .resolve_labels(x, labels)

  # Determine gene list
  all_genes <- if (gene_type == "tf") x@tf_ids else x@target_ids
  if (is.null(gene_names)) {
    genes <- sort(all_genes)
  } else {
    genes <- intersect(as.character(gene_names), all_genes)
    if (length(genes) == 0L)
      stop(sprintf("None of the specified %ss found in the experiment.", gene_type))
  }

  n_genes <- length(genes)
  message(sprintf("Plotting %d %s weight barplot(s)...", n_genes,
                  if (gene_type == "tf") "TF" else "target"))

  # Color map: label_name -> color (ordered by label integer)
  labels_chr <- as.character(labels)
  lab_names  <- unname(x@label_names[labels_chr])
  col_map    <- stats::setNames(unname(x@colors[labels_chr]), lab_names)

  # Convert all weights to long format once, then subset
  long_df <- .weights_to_long(x)

  # Subset to requested labels
  long_df <- long_df[long_df$label %in% labels , , drop = FALSE]

  # Grid setup
  if (is.null(grid)) {
    grid_cols <- 1L
    grid_rows <- n_genes
    plots_per_page <- n_genes
  } else {
    grid_rows <- as.integer(grid[1])
    grid_cols <- as.integer(grid[2])
    plots_per_page <- grid_rows * grid_cols
  }

  fig_height_per_row <- if (gene_type == "tf") 5 else 4
  if (is.null(height)) height <- fig_height_per_row * grid_rows

  # Build all plots
  plot_list <- list()
  for (gene in genes) {
    # Subset long_df for this gene
    if (gene_type == "tf") {
      df_gene <- long_df[long_df$tf == gene, , drop = FALSE]
    } else {
      df_gene <- long_df[long_df$target == gene, , drop = FALSE]
    }

    p <- tryCatch(
      .make_weight_barplot(
        df          = df_gene,
        gene_name   = gene,
        gene_type   = gene_type,
        col_map     = col_map,
        top_n       = if (gene_type == "tf") top_n else NULL,
        allowed_genes = if (gene_type == "tf") allowed_targets else NULL
      ),
      error = function(e) {
        message(sprintf("  Skipping %s: %s", gene, conditionMessage(e)))
        ggplot2::ggplot() +
          ggplot2::annotate("text", x = 0.5, y = 0.5,
                            label = paste("Error plotting\n", gene),
                            size = 4, hjust = 0.5, vjust = 0.5) +
          ggplot2::theme_void()
      }
    )
    plot_list[[length(plot_list) + 1L]] <- p
  }

  # Paginate
  pages <- .paginate_plots(plot_list, plots_per_page, grid_rows, grid_cols)

  # Output
  if (save) {
    default_fname <- sprintf("%s_weights.pdf",
                             if (gene_type == "tf") "TF" else "target")
    fname <- filename %||% default_fname
    fpath <- file.path(out_dir, fname)
    dir.create(dirname(fpath), recursive = TRUE, showWarnings = FALSE)
    .save_pages(pages, fpath, width, height)
    message("Saved weight barplots to: ", fpath)
  } else {
    .draw_pages(pages)
  }

  invisible(pages)
}

# ---- Public: plot_tf_weights --------------------------------------------

#' Plot weight barplots for transcription factors
#'
#' For each TF, shows a grouped bar chart of regulatory weights across all
#' target genes, coloured by phenotype label. Mirrors
#' \code{SimiCVisualization.plot_tf_weights} from the Python package.
#'
#' @param x \code{SimiCvizExperiment}
#' @param tf_names character vector of TFs to plot (default: all)
#' @param labels integer or character vector of labels (default: all)
#' @param top_n integer; show only the top N targets by mean absolute weight
#'   (default 50). Applied after \code{allowed_targets} filtering.
#' @param allowed_targets character vector or NULL; restrict targets to this
#'   pre-filtered set before applying \code{top_n}. Useful for passing targets
#'   that pass an external significance filter.
#' @param grid integer vector \code{c(nrow, ncol)} per page (default
#'   \code{c(4, 1)}). Set \code{NULL} for a single page.
#' @param save logical; save to PDF (default FALSE)
#' @param filename PDF filename (default \code{"TF_weights.pdf"})
#' @param out_dir output directory (default \code{getwd()})
#' @param width,height page dimensions in inches
#'
#' @return Invisibly, a list of page grobs.
#' @export
plot_tf_weights <- function(x,
                            tf_names        = NULL,
                            labels          = NULL,
                            top_n           = 50L,
                            allowed_targets = NULL,
                            grid            = c(4L, 1L),
                            save            = FALSE,
                            filename        = NULL,
                            out_dir         = getwd(),
                            width           = 16,
                            height          = NULL) {
  .plot_weights_paginated(x,
                          gene_names      = tf_names,
                          gene_type       = "tf",
                          labels          = labels,
                          top_n           = top_n,
                          allowed_targets = allowed_targets,
                          grid            = grid,
                          save            = save,
                          filename        = filename %||% "TF_weights.pdf",
                          out_dir         = out_dir,
                          width           = width,
                          height          = height)
}

# ---- Public: plot_target_weights ----------------------------------------

#' Plot weight barplots for target genes
#'
#' For each target gene, shows a grouped bar chart of incoming regulatory
#' weights from all TFs, coloured by phenotype label. Mirrors
#' \code{SimiCVisualization.plot_target_weights} from the Python package.
#'
#' @param x \code{SimiCvizExperiment}
#' @param target_names character vector of targets to plot (default: all)
#' @param labels integer or character vector of labels (default: all)
#' @param grid integer vector \code{c(nrow, ncol)} per page (default
#'   \code{c(4, 1)}). Set \code{NULL} for a single page.
#' @param save logical; save to PDF (default FALSE)
#' @param filename PDF filename (default \code{"target_weights.pdf"})
#' @param out_dir output directory (default \code{getwd()})
#' @param width,height page dimensions in inches
#'
#' @return Invisibly, a list of page grobs.
#' @export
plot_target_weights <- function(x,
                                target_names = NULL,
                                labels       = NULL,
                                grid         = c(4L, 1L),
                                top_n      = NULL,
                                save         = FALSE,
                                filename     = NULL,
                                out_dir      = getwd(),
                                width        = 12,
                                height       = NULL) {
  .plot_weights_paginated(x,
                          gene_names = target_names,
                          gene_type  = "target",
                          labels     = labels,
                          top_n      = top_n,
                          grid       = grid,
                          save       = save,
                          filename   = filename %||% "target_weights.pdf",
                          out_dir    = out_dir,
                          width      = width,
                          height     = height)
}
