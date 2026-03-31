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

# PALETTES HEATMAPS  ------------
#' Build a ggplot2 fill scale from a cmap specification
#'
#' @param cmap Palette spec (NULL, viridis name, single colour, or vector).
#' @return A ggplot2 scale object.
#' @keywords internal
#' @noRd
.build_ggplot_fill_scale <- function(cmap) {
  viridis_names <- c("viridis", "magma", "plasma", "inferno", "cividis",
                     "rocket", "mako", "turbo")

  # Default: viridis-like

  if (is.null(cmap)) {
    return(ggplot2::scale_fill_gradientn(
      colours = c("#440154", "#31688e", "#35b779", "#fde725")
    ))
  }

  # Named viridis palette
  if (is.character(cmap) && length(cmap) == 1L &&
      tolower(cmap) %in% viridis_names) {
    palette_name <- tolower(cmap)
    if (requireNamespace("viridisLite", quietly = TRUE)) {
      cols <- viridisLite::viridis(256, option = palette_name)
      return(ggplot2::scale_fill_gradientn(colours = cols))
    }
    # Fallback approximations
    fallback <- list(
      viridis = c("#440154", "#31688e", "#35b779", "#fde725"),
      magma   = c("#000004", "#51127c", "#b73779", "#fc8961", "#fcfdbf"),
      plasma  = c("#0d0887", "#7e03a8", "#cc4778", "#f89540", "#f0f921"),
      inferno = c("#000004", "#420a68", "#932667", "#dd513a", "#fca50a", "#fcffa4"),
      cividis = c("#00224e", "#414d6b", "#7b7b78", "#b8a951", "#fdea45"),
      rocket  = c("#03051a", "#4c1d4e", "#a11a5b", "#e04f39", "#faebdd"),
      mako    = c("#0b0405", "#2a1858", "#245f8a", "#30b09a", "#def5e5"),
      turbo   = c("#30123b", "#4662d7", "#35abf8", "#1ae4b6", "#72fe5e",
                  "#c8ef34", "#faba39", "#f66b19", "#ca240e", "#7a0403")
    )
    anchors <- fallback[[palette_name]] %||% fallback[["viridis"]]
    message(sprintf("viridisLite not installed; using approximate '%s' palette.", palette_name))
    return(ggplot2::scale_fill_gradientn(colours = anchors))
  }

  # Single colour: white → colour
  if (is.character(cmap) && length(cmap) == 1L) {
    return(ggplot2::scale_fill_gradient(low = "white", high = cmap))
  }

  # Custom vector of colours
  if (is.character(cmap) && length(cmap) >= 2L) {
    return(ggplot2::scale_fill_gradientn(colours = cmap))
  }

  warning("Unrecognised `cmap`; using default viridis-like palette.")
  ggplot2::scale_fill_gradientn(
    colours = c("#440154", "#31688e", "#35b779", "#fde725")
  )
}


#' Build a palette function (value -> colour) matching the ggplot scale
#' @noRd
.build_palette_function <- function(cmap, domain) {
  # Remove NA values from domain
  domain <- domain[!is.na(domain)]
  
  if (length(domain) == 0) {
    stop("No valid values in domain for palette function.")
  }
  viridis_names <- c("viridis", "magma", "plasma", "inferno", "cividis",
                     "rocket", "mako", "turbo")

  if (is.null(cmap)) {
    cols <- c("#440154", "#31688e", "#35b779", "#fde725")
    return(scales::col_numeric(cols, domain = domain))
  }

  if (is.character(cmap) && length(cmap) == 1L &&
      tolower(cmap) %in% viridis_names) {

    palette_name <- tolower(cmap)

    if (requireNamespace("viridisLite", quietly = TRUE)) {
      cols <- viridisLite::viridis(256, option = palette_name)
    } else {
      fallback <- list(
        viridis = c("#440154", "#31688e", "#35b779", "#fde725"),
        magma   = c("#000004", "#51127c", "#b73779", "#fc8961", "#fcfdbf"),
        plasma  = c("#0d0887", "#7e03a8", "#cc4778", "#f89540", "#f0f921"),
        inferno = c("#000004", "#420a68", "#932667", "#dd513a", "#fca50a", "#fcffa4"),
        cividis = c("#00224e", "#414d6b", "#7b7b78", "#b8a951", "#fdea45"),
        rocket  = c("#03051a", "#4c1d4e", "#a11a5b", "#e04f39", "#faebdd"),
        mako    = c("#0b0405", "#2a1858", "#245f8a", "#30b09a", "#def5e5"),
        turbo   = c("#30123b", "#4662d7", "#35abf8", "#1ae4b6", "#72fe5e",
                    "#c8ef34", "#faba39", "#f66b19", "#ca240e", "#7a0403")
      )
      cols <- fallback[[palette_name]] %||% fallback[["viridis"]]
    }

    return(scales::col_numeric(cols, domain = domain))
  }

  if (is.character(cmap) && length(cmap) == 1L) {
    return(scales::col_numeric(c("white", cmap), domain = domain))
  }

  if (is.character(cmap) && length(cmap) >= 2L) {
    return(scales::col_numeric(cmap, domain = domain))
  }

  cols <- c("#440154", "#31688e", "#35b779", "#fde725")
  scales::col_numeric(cols, domain = domain)
}



#' @noRd
.get_text_color <- function(bg_color) {
  contrast_black <- colorspace::contrast_ratio(bg_color, "black")
  contrast_white <- colorspace::contrast_ratio(bg_color, "white")

  ifelse(contrast_black > contrast_white, "black", "white")
}