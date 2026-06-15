# ----------------------------- 2) Plot helpers -------------------------------
theme_pub <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.2, color = "grey88"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(color = "grey25"),
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey92", color = "grey70"),
      strip.text = element_text(face = "bold")
    )
}

save_gg <- function(p, name, width = 7, height = 5, dpi = 320) {
  pdf_file <- file.path(FIGDIR, paste0(name, ".pdf"))
  png_file <- file.path(FIGDIR, paste0(name, ".png"))
  ggsave(pdf_file, p, width = width, height = height, units = "in", useDingbats = FALSE)
  ggsave(png_file, p, width = width, height = height, units = "in", dpi = dpi)
  invisible(c(pdf_file, png_file))
}


# scales::label_number_si() was deprecated in scales 1.2.0 and is defunct in
# recent versions. This wrapper uses the modern API when available and falls
# back to comma labels on older/incomplete scales installations.
label_number_compact_safe <- function(...) {
  if (exists("cut_si", where = asNamespace("scales"), inherits = FALSE)) {
    scales::label_number(scale_cut = scales::cut_si(""), ...)
  } else {
    scales::label_comma(...)
  }
}

# Robust discrete palette helper. The grid error
# "the 'gpar' element 'fill' must not have length 0" is commonly triggered
# when a manual fill/color scale receives an empty named vector. These helpers
# always return either a non-empty named palette for the values actually present
# in the plotted data, or NULL. ggplot can safely add NULL layers.
safe_named_palette <- function(x, preferred = NULL, fallback_palette = scales::hue_pal()) {
  vals <- unique(as.character(x))
  vals <- vals[!is.na(vals) & nzchar(vals)]
  vals <- sort(vals)
  if (length(vals) == 0) return(NULL)
  
  if (!is.null(preferred) && length(preferred) > 0 && all(vals %in% names(preferred))) {
    pal <- preferred[vals]
  } else {
    pal <- setNames(fallback_palette(length(vals)), vals)
  }
  
  if (length(pal) == 0 || any(is.na(pal)) || any(!nzchar(as.character(pal)))) {
    pal <- setNames(fallback_palette(length(vals)), vals)
  }
  pal
}

manual_fill_or_null <- function(pal, drop = FALSE) {
  if (is.null(pal) || length(pal) == 0) return(NULL)
  ggplot2::scale_fill_manual(values = pal, drop = drop)
}

manual_color_or_null <- function(pal, drop = FALSE) {
  if (is.null(pal) || length(pal) == 0) return(NULL)
  ggplot2::scale_color_manual(values = pal, drop = drop)
}

save_base_plot <- function(name, expr, width = 7, height = 5, dpi = 320) {
  expr_sub <- substitute(expr)
  pdf_file <- file.path(FIGDIR, paste0(name, ".pdf"))
  png_file <- file.path(FIGDIR, paste0(name, ".png"))
  pdf(pdf_file, width = width, height = height, useDingbats = FALSE)
  eval(expr_sub, envir = parent.frame())
  dev.off()
  png(png_file, width = width, height = height, units = "in", res = dpi)
  eval(expr_sub, envir = parent.frame())
  dev.off()
  invisible(c(pdf_file, png_file))
}

save_pheatmap <- function(mat, name, width = 7, height = 7, ...) {
  pdf_file <- file.path(FIGDIR, paste0(name, ".pdf"))
  png_file <- file.path(FIGDIR, paste0(name, ".png"))
  pdf(pdf_file, width = width, height = height, useDingbats = FALSE)
  print(pheatmap::pheatmap(mat, silent = TRUE, ...))
  dev.off()
  png(png_file, width = width, height = height, units = "in", res = 320)
  print(pheatmap::pheatmap(mat, silent = TRUE, ...))
  dev.off()
  invisible(c(pdf_file, png_file))
}

message_step <- function(...) message("\n[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
