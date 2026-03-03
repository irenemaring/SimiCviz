# ---- Internal helpers ----------------------------------------------------

#' Collect weight data for a single gene (TF or target)
#'
#' @param x SimiCvizExperiment
#' @param gene_name character; TF or target gene name
#' @param gene_type character; "tf" or "target"
#' @param labels integer vector of labels
#' @return data.frame with columns: gene, weight, abs_weight, label, label_name
#' @keywords internal
.collect_weight_data <- function(x, gene_name, gene_type, labels) {
  if (length(x@weights) == 0L) stop("No weight matrices found in the experiment.")

  tf_ids     <- x@tf_ids
  target_ids <- x@target_ids

  if (gene_type == "tf") {
    gene_index  <- which(tf_ids == gene_name)
    partner_ids <- target_ids
  } else {
    gene_index  <- which(target_ids == gene_name)
    partner_ids <- tf_ids
  }

  if (length(gene_index) == 0L) {
    stop(sprintf("Gene '%s' not found in %s_ids.", gene_name,
                 if (gene_type == "tf") "tf" else "target"))
  }
  if (length(gene_index) > 1L) {
    warning(sprintf("Gene '%s' matches multiple indices; using first.", gene_name))
    gene_index <- gene_index[1L]
  }

  rows <- list()
  for (lab in labels) {
    lab_chr <- as.character(lab)
    W <- x@weights[[lab_chr]]
    if (is.null(W)) next

    # W is tfs x targets (rows = TFs, cols = targets)
    # mirrors Python: weight_dic[label][tf_index, :] for TF
    #                 weight_dic[label][:, target_index] for target
    if (gene_type == "tf") {
      weights_vec <- as.numeric(W[gene_index, ])   # row = TF → all targets
    } else {
      weights_vec <- as.numeric(W[, gene_index])   # col = target → all TFs
    }

    if (length(weights_vec) != length(partner_ids)) {
      warning(sprintf(
        "Length mismatch for '%s' label %s: weights (%d) vs partners (%d). Skipping.",
        gene_name, lab_chr, length(weights_vec), length(partner_ids)
      ))
      next
    }

    nonzero <- which(weights_vec != 0)
    if (length(nonzero) == 0L) next

    rows[[length(rows) + 1L]] <- data.frame(
      gene       = partner_ids[nonzero],
      weight     = weights_vec[nonzero],
      abs_weight = abs(weights_vec[nonzero]),
      label      = lab,
      label_name = unname(x@label_names[lab_chr]),
      stringsAsFactors = FALSE
    )
  }

  if (length(rows) == 0L) {
    return(data.frame(gene = character(), weight = numeric(),
                      abs_weight = numeric(), label = integer(),
                      label_name = character(), stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}

#' Render a single weight barplot on a ggplot
#'
#' @param df data.frame from .collect_weight_data
#' @param gene_name character
#' @param gene_type character; "tf" or "target"
#' @param labels integer vector of labels used
#' @param col_map named character vector of colours (keyed by label_name)
#' @param top_n integer or NULL; keep only top N partners by mean abs weight
#' @return ggplot object
#' @keywords internal
.make_weight_barplot <- function(df, gene_name, gene_type, labels,
                                 col_map, top_n = NULL) {
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

  # Optionally keep top N partners
  if (!is.null(top_n)) {
    top_genes <- names(sort(tapply(df$abs_weight, df$gene, mean),
                            decreasing = TRUE))[seq_len(min(top_n, length(unique(df$gene))))]
    df <- df[df$gene %in% top_genes, , drop = FALSE]
  }

  # Order genes by max abs weight
  gene_order <- names(sort(tapply(df$abs_weight, df$gene, max), decreasing = TRUE))

  # ---- Complete the grid: every gene x label_name combination must exist ----
  all_label_names <- names(col_map)
  full_grid <- expand.grid(gene = gene_order, label_name = all_label_names,
                           stringsAsFactors = FALSE)
  df <- merge(full_grid, df[, c("gene", "label_name", "weight", "abs_weight")],
              by = c("gene", "label_name"), all.x = TRUE)
  df$weight[is.na(df$weight)]         <- 0
  df$abs_weight[is.na(df$abs_weight)] <- 0

  df$gene       <- factor(df$gene,       levels = gene_order)
  df$label_name <- factor(df$label_name, levels = all_label_names)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$gene,
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
#' @param grid integer vector c(nrow, ncol) per page; NULL = single page
#' @param save logical
#' @param filename PDF filename
#' @param out_dir output directory
#' @param width,height page dimensions in inches
#' @return Invisibly, list of page grobs.
#' @keywords internal
.plot_weights_paginated <- function(x,
                                    gene_names = NULL,
                                    gene_type  = c("tf", "target"),
                                    labels     = NULL,
                                    top_n      = NULL,
                                    grid       = c(4L, 1L),
                                    save       = FALSE,
                                    filename   = NULL,
                                    out_dir    = getwd(),
                                    width      = 16,
                                    height     = NULL) {
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

  # Color map: label_name -> color
  lab_names  <- unname(x@label_names[as.character(labels)])
  col_map    <- stats::setNames(unname(x@colors[as.character(labels)]), lab_names)

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
    df <- tryCatch(
      .collect_weight_data(x, gene, gene_type, labels),
      error = function(e) {
        message(sprintf("  Skipping %s: %s", gene, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(df)) df <- data.frame(gene = character(), weight = numeric(),
                                       abs_weight = numeric(), label = integer(),
                                       label_name = character(),
                                       stringsAsFactors = FALSE)
    p <- .make_weight_barplot(df, gene, gene_type, labels, col_map,
                              top_n = if (gene_type == "tf") top_n else NULL)
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
    grDevices::pdf(fpath, width = width, height = height, onefile = TRUE)
    for (pg in pages) {
      grid::grid.newpage()
      grid::grid.draw(pg)
    }
    grDevices::dev.off()
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
#'   (default 50)
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
                            tf_names = NULL,
                            labels   = NULL,
                            top_n    = 50L,
                            grid     = c(4L, 1L),
                            save     = FALSE,
                            filename = NULL,
                            out_dir  = getwd(),
                            width    = 16,
                            height   = NULL) {
  .plot_weights_paginated(x,
                          gene_names = tf_names,
                          gene_type  = "tf",
                          labels     = labels,
                          top_n      = top_n,
                          grid       = grid,
                          save       = save,
                          filename   = filename %||% "TF_weights.pdf",
                          out_dir    = out_dir,
                          width      = width,
                          height     = height)
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
                                save         = FALSE,
                                filename     = NULL,
                                out_dir      = getwd(),
                                width        = 12,
                                height       = NULL) {
  .plot_weights_paginated(x,
                          gene_names = target_names,
                          gene_type  = "target",
                          labels     = labels,
                          top_n      = NULL,
                          grid       = grid,
                          save       = save,
                          filename   = filename %||% "target_weights.pdf",
                          out_dir    = out_dir,
                          width      = width,
                          height     = height)
}
