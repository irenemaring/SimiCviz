#' Plot dissimilarity heatmap
#'
#' Displays a heatmap of TF regulatory dissimilarity scores, optionally
#' broken down by cell group. Uses \pkg{ggplot2} tile geometry.
#'
#' @param x A \code{SimiCvizExperiment} object.
#' @param top_n Integer; number of top TFs to display (default: all).
#' @param sort_by Column name to sort TFs by.
#'   Default: \code{"MinMax_score"} (no groups) or \code{"mean_score"} (groups).
#' @param dissim_df Optional pre-computed dissimilarity data.frame. If
#'   \code{NULL}, \code{\link{calculate_dissimilarity}} is called internally.
#' @param cmap Colour palette specification. Can be:
#'   \itemize{
#'     \item A viridis palette name: \code{"viridis"}, \code{"magma"},
#'       \code{"plasma"}, \code{"inferno"}, or \code{"cividis"}
#'       (requires \pkg{viridisLite}).
#'     \item A single colour string (gradient from white to that colour).
#'     \item A character vector of 2+ colours for a custom gradient.
#'     \item \code{NULL} (default): built-in viridis-like gradient.
#'   }
#' @param show_values Logical; annotate cells with numeric values (default \code{TRUE}).
#' @param save Logical; save the plot to a PDF file (default \code{TRUE}).
#' @param out_dir Output directory for the PDF (default: working directory).
#' @param filename Custom filename (default: auto-generated).
#' @param width,height PDF dimensions in inches.
#' @param ... Additional arguments passed to \code{calculate_dissimilarity}.
#'
#' @return Invisibly, a list with \code{plot} (the ggplot object) and
#'   \code{data} (the dissimilarity data.frame).
#'
#' @examples
#' \dontrun{
#'   plot_dissimilarity_heatmap(simic, top_n = 20)
#'   plot_dissimilarity_heatmap(simic, top_n = 15, cmap = "magma", save = FALSE)
#'   plot_dissimilarity_heatmap(simic, cmap = c("white", "red", "darkred"))
#' }
#' @import ggplot2
#' @rdname plot_dissimilarity_heatmap
#' @export
plot_dissimilarity_heatmap <- function(x,
                                       top_n = NULL,
                                       sort_by = NULL,
                                       dissim_df = NULL,
                                       cmap = NULL,
                                       show_values = TRUE,
                                       save = TRUE,
                                       out_dir = getwd(),
                                       filename = NULL,
                                       width = NULL,
                                       height = NULL,
                                       ...
                                       ) {

  if (!is.SimiCvizExperiment(x)) {
    stop("plot_dissimilarity_heatmap: 'x' must be a SimiCvizExperiment.")
  }

  # Compute or use pre-computed scores
  if (is.null(dissim_df)) {
    dissim_df <- calculate_dissimilarity(x, verbose = FALSE, ...)
  }
  if (is.null(dissim_df) || nrow(dissim_df) == 0L) {
    stop("No dissimilarity scores to plot.")
  }

  # Sort
  if (is.null(sort_by)) {
    sort_by <- if ("mean_score" %in% colnames(dissim_df)) "mean_score" else colnames(dissim_df)[1]
  }
  if (!(sort_by %in% colnames(dissim_df))) {
    warning(sprintf("sort_by='%s' not found; using first column.", sort_by))
    sort_by <- colnames(dissim_df)[1]
  }
  dissim_df <- dissim_df[order(dissim_df[[sort_by]], decreasing = FALSE), , drop = FALSE]

  # Top N (after sorting descending, take top_n then re-sort for plot)
  if (!is.null(top_n) && top_n < nrow(dissim_df)) {
    # Take the top_n highest (they are at the bottom after ascending sort)
    dissim_df <- utils::tail(dissim_df, top_n)
  }

  n_rows <- nrow(dissim_df)
  n_cols <- ncol(dissim_df)

  # Prettify column names for display
  col_display <- colnames(dissim_df)
  col_display[col_display == "mean_score"]   <- "Mean"
  col_display[col_display == "MinMax_score"] <- "MinMax Score"
  col_map <- stats::setNames(col_display, colnames(dissim_df))

  # Convert to long format for ggplot
  dissim_df$TF <- factor(rownames(dissim_df), levels = rownames(dissim_df))
  long_df <- reshape2::melt(dissim_df, id.vars = "TF",
                            variable.name = "Metric", value.name = "Score")
  # Apply display names
  long_df$Metric <- factor(col_map[as.character(long_df$Metric)],
                           levels = col_display)

  # Build colour scale
  fill_scale <- .build_ggplot_fill_scale(cmap)

  # Auto dimensions
  if (is.null(width))  width  <- max(5, 1.5 + n_cols * 1.5)
  if (is.null(height)) height <- max(5, n_rows * 0.35 + 2)

  # Title
  title_txt <-  "Regulatory Dissimilarity Scores"


  # Build palette function (must match fill scale!)
  pal_fun <- .build_palette_function(cmap, domain = long_df$Score)

  # Compute actual fill colors per cell
  long_df$fill_color <- pal_fun(long_df$Score)

  # Compute adaptive text colors using colorspace
  long_df$text_color <- .get_text_color(long_df$fill_color)
  p <- ggplot2::ggplot(long_df, ggplot2::aes(x = Metric, y = TF, fill = Score)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    fill_scale +
    ggplot2::labs(title = title_txt,
                  x = NULL, y = "Transcription Factors",
                  fill = "Dissimilarity\nScore") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y  = ggplot2::element_text(size = 8),
      plot.title    = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
      panel.grid    = ggplot2::element_blank(),
      legend.position = "right"
    )


  if (show_values) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.4f", Score),
                  colour = text_color),
      size = 2.8,
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_identity()
  }
  
  if (save) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    fname <- if (is.null(filename)) "dissimilarity_heatmap.pdf" else filename
    fpath <- file.path(out_dir, fname)
    ggplot2::ggsave(fpath, plot = p, width = width, height = height)
    message(sprintf("Saved dissimilarity heatmap to %s", fpath))
  } else {
    print(p)
  }

  # Remove the TF column we added for melting
  dissim_df$TF <- NULL

  invisible(list(plot = p, data = dissim_df))
}