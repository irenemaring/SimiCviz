#' Extract stacked TF-target weight matrix across labels
#'
#' For a given TF, extracts the regulatory weight for every target gene in each
#' label. Targets whose adjusted R² falls below \code{r2_threshold} in a given
#' label receive \code{NA} (analogous to Python's \code{get_TF_network(...,
#' stacked=TRUE)}).
#'
#' @param x A \code{SimiCvizExperiment} object.
#' @param tf_name Character; name of the transcription factor.
#' @param labels Integer vector of labels to include (default: all).
#' @param r2_threshold Numeric or \code{NULL}. If non-NULL, weights for
#'   targets with adjusted R² below this value are set to \code{NA}.
#' @return A data.frame with target genes as rows and one column per label
#'   (named by \code{label_names}).
#' @examples
#'  # Example usage
#'   simic <- readRDS(system.file("extdata", "simic_full.rds", package = "SimiCviz"))
#'   network <- get_tf_network(simic, simic@tf_ids[1], labels = c(1, 2), r2_threshold = 0.7)
#'   print(network)
#'
#' @export

get_tf_network <- function(x, tf_name, labels = NULL, r2_threshold = NULL) {

  if (length(x@weights) == 0L)
    stop("No weight matrices found in the experiment.")

  labels <- .resolve_labels(x, labels)

  # Validate TF

  if (!tf_name %in% x@tf_ids)
    stop(sprintf("TF '%s' not found in the experiment.", tf_name))

  tf_idx <- which(x@tf_ids == tf_name)
  target_ids <- x@target_ids

  # Build result data.frame
  result <- data.frame(row.names = target_ids)

  for (lab in labels) {
    lab_chr  <- as.character(lab)
    lab_name <- unname(x@label_names[lab_chr])
    W <- as.matrix(x@weights[[lab_chr]])

    # Extract the TF row (TFs are rows, targets are columns)
    weights_vec <- W[tf_idx, ]

    # Apply R2 filtering if requested
    if (!is.null(r2_threshold) && length(x@meta$adjusted_r_squared) > 0L) {
      r2_vals <- x@meta$adjusted_r_squared[[lab_chr]]
      if (!is.null(r2_vals)) {
        below_thresh <- r2_vals < r2_threshold
        weights_vec[below_thresh] <- NA_real_
      }
    } else if (!is.null(r2_threshold) && length(x@meta$adjusted_r_squared) == 0L) {
        warning(sprintf("No adjusted R\u00B2 values found for label '%s'; skipping R\u00B2 filtering.", lab_name))
      }

    result[[lab_name]] <- weights_vec
  }

  result
}


# ---- Public: plot_tf_network_heatmap ------------------------------------

#' Plot heatmap of a TF regulatory network across phenotypes
#'
#' Displays a heatmap showing the regulatory weights of a single transcription
#' factor across its top target genes, with one column per phenotype label.
#' Targets that fail the adjusted-R² filter are shown as grey cells labelled
#' \emph{"< R² threshold"}, mirroring the Python
#' \code{SimiCVisualization.plot_tf_network_heatmap()} method.
#'
#' @param x A \code{SimiCvizExperiment} object.
#' @param tf_name Character; name of the transcription factor to plot.
#' @param labels Integer or character vector of labels to include (default:
#'   all).
#' @param top_n Integer; number of top target genes to display, ranked by
#'   maximum absolute weight across labels (default: 10).
#' @param r2_threshold Numeric or \code{NULL}. Targets with adjusted R² below
#'   this value receive \code{NA} and are rendered as grey tiles. Passed to
#'   \code{.get_tf_network()}.
#' @param cmap Colour palette specification, identical to
#'   \code{\link{plot_dissimilarity_heatmap}}. Recommended diverging palettes
#'   such as \code{"RdBu_r"} (red-white-blue, reversed) or a two-colour
#'   vector like \code{c("blue", "white", "red")}. Default: \code{c("blue",
#'   "white", "red")}.
#' @param show_values Logical; annotate each cell with its weight or
#'   \emph{"< R² threshold"} text (default \code{TRUE}).
#' @param save Logical; save the plot to a PDF file (default \code{TRUE}).
#' @param filename Custom filename (default: auto-generated from
#'   \code{tf_name}).
#' @param out_dir Output directory for the PDF (default: working directory).
#' @param width,height PDF dimensions in inches (auto-calculated if
#'   \code{NULL}).
#'
#' @return Invisibly, a list with components:
#'   \describe{
#'     \item{plot}{The \code{ggplot} object.}
#'     \item{data}{The (filtered) data.frame of weights used for the plot.}
#'   }
#'
#' @examples
#'   simic <- readRDS(system.file("extdata", "simic_full.rds", package = "SimiCviz"))
#'   plot_tf_network_heatmap(simic, "Pms1")
#'
#'   # Custom palette and top targets
#'   plot_tf_network_heatmap(simic, "Ets1", top_n = 15,
#'                           cmap = "RdBu_r", r2_threshold = 0.7)
#'
#'   # Only specific labels, no saving
#'   plot_tf_network_heatmap(simic, "Gli3", labels = c(0, 3),
#'                           save = FALSE)
#'
#' @import ggplot2
#' @export
plot_tf_network_heatmap <- function(x,
                                    tf_name,
                                    labels        = NULL,
                                    top_n         = 10L,
                                    r2_threshold  = NULL,
                                    cmap          = c("blue", "white", "red"),
                                    show_values   = TRUE,
                                    save          = FALSE,
                                    filename      = NULL,
                                    out_dir       = getwd(),
                                    width         = NULL,
                                    height        = NULL) {

  if (!is.SimiCvizExperiment(x))
    stop("plot_tf_network_heatmap: 'x' must be a SimiCvizExperiment.")
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required.")

  # ---- Extract network data ------------------------------------------------
  network <- get_tf_network(x,
                             tf_name      = tf_name,
                             labels       = labels,
                             r2_threshold = r2_threshold)

  if (nrow(network) == 0L)
    stop(sprintf("No network data found for TF '%s'.", tf_name))

  # ---- Select top N targets by max absolute weight --------------------------
  max_abs <- apply(abs(network), 1, max, na.rm = TRUE)
  # Replace -Inf (all-NA rows) with 0

  max_abs[!is.finite(max_abs)] <- 0

  top_n <- min(top_n, nrow(network))
  top_targets <- names(sort(max_abs, decreasing = TRUE))[seq_len(top_n)]
  network <- network[top_targets, , drop = FALSE]

  # ---- Convert to long format -----------------------------------------------
  network$target <- factor(rownames(network), levels = rev(rownames(network)))
  long_df <- reshape2::melt(network, id.vars = "target",
                            variable.name = "condition",
                            value.name    = "weight")
  long_df$condition <- factor(long_df$condition, levels = colnames(network)[colnames(network) != "target"])
  long_df$is_na     <- is.na(long_df$weight)

  # ---- Symmetric colour limits -----------------------------------------------
  non_na_vals <- long_df$weight[!long_df$is_na]
  if (length(non_na_vals) == 0L) {
    non_na_vals <- 0
  }

  # ---- Build fill scale (diverging, with NA colour) -------------------------
  fill_scale <- .build_ggplot_fill_scale(cmap)

  # ---- Adaptive text colours -------------------------------------------------
  pal_fun <- .build_palette_function(cmap, domain = non_na_vals)
  long_df$fill_color <- ifelse(long_df$is_na, "#A9A9A9",
                               pal_fun(long_df$weight))
  long_df$text_color <- ifelse(long_df$is_na, "black",
                               .get_text_color(long_df$fill_color))

  # ---- Cell label text -------------------------------------------------------
  long_df$label_text <- ifelse(
    long_df$is_na,
    "< R\u00B2 threshold",
    sprintf("%.2f", long_df$weight)
  )

  # ---- Auto dimensions -------------------------------------------------------
  n_rows <- length(top_targets)
  n_cols <- ncol(network) - 1L   # exclude the "target" column
  if (is.null(width))  width  <- max(8, n_cols * 1.5 + 3)
  if (is.null(height)) height <- max(6, n_rows * 0.4 + 2)

  # ---- Build plot -------------------------------------------------------------
  title_txt <- sprintf(" %s Regulatory Network", tf_name)
  subtitle_txt <- sprintf( "Top %d targets across phenotypes", top_n)

  p <- ggplot2::ggplot(long_df,
                       ggplot2::aes(x = .data$condition,
                                    y = .data$target,
                                    fill = .data$weight)) +
    # Grey tiles for NA values (drawn first, underneath)
    ggplot2::geom_tile(data = long_df[long_df$is_na, , drop = FALSE],
                       ggplot2::aes(x = .data$condition, y = .data$target),
                       fill = "darkgrey", colour = "grey80", linewidth = 0.5,
                       inherit.aes = FALSE) +
    # Coloured tiles for non-NA values
    ggplot2::geom_tile(data = long_df[!long_df$is_na, , drop = FALSE],
                       colour = "grey80", linewidth = 0.5) +
    fill_scale +
    ggplot2::labs(x = NULL, y = "Target Genes",
                  fill = "Weight",
                  title = title_txt,
                  subtitle = subtitle_txt) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x   = ggplot2::element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y   = ggplot2::element_text(size = 9),
      plot.title     = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      panel.grid     = ggplot2::element_blank(),
      legend.position = "right"
    )

  if (show_values) {
    p <- p +
      ggplot2::geom_text(
        ggplot2::aes(label = .data$label_text,
                     colour = .data$text_color),
        size = 2.8,
        show.legend = FALSE
      ) +
      ggplot2::scale_colour_identity()
  }

  # ---- Save / display ---------------------------------------------------------
  if (save) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    fname <- if (is.null(filename)) {
      sprintf("network_%s_heatmap.pdf", tf_name)
    } else {
      filename
    }
    fpath <- file.path(out_dir, fname)
    ggplot2::ggsave(fpath, plot = p, width = width, height = height)
    message(sprintf("Saved TF network heatmap to %s", fpath))
  } else {
    print(p)
  }

  # Return data without the helper "target" column
  network$target <- NULL

  invisible(p)
}
