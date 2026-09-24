small_axis <- function(label = NULL, fontsize = 7, arrow_length = 10, label_offset = 1, fix_coord = TRUE, remove_axes = TRUE,
                       arrow_spec = grid::arrow(ends = "both", type = "closed", angle = 20, length = ggplot2::unit(arrow_length / 7, units)),
                       units = "mm", ...){
  coord <- if(fix_coord){
    ggplot2::coord_fixed(clip = "off", ...)
  }else{
    NULL
  }
  axis_theme <- if(remove_axes){
    ggplot2::theme(axis.line = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank())
  }else{
    NULL
  }
  lines <- ggplot2::annotation_custom(grid::polylineGrob(x = ggplot2::unit(c(0, 0, arrow_length), units), y = ggplot2::unit(c(arrow_length, 0, 0), units),
                                                         gp = grid::gpar(fill = "black"),
                                                         arrow = arrow_spec))
  text <- if(! is.null(label)){
    ggplot2::annotation_custom(grid::textGrob(label = label, gp = grid::gpar(fontsize = fontsize),
                                              x = ggplot2::unit(label_offset, units), y = ggplot2::unit(label_offset, units), hjust = 0, vjust = 0))
  }else{
    NULL
  }
  list(coord, axis_theme, lines, text)
}

signif_to_zero <- function(x, digits = 6){
  n_signif_digits <- digits - ceiling(log10(abs(x)))
  sign(x) * floor(abs(x) * 10^n_signif_digits) / 10^n_signif_digits
}

scale_color_de_gradient <- function(abs_max, mid_width = 0.1, ..., oob = scales::squish, limits = c(-1, 1) * abs_max, breaks = c(-1, 0, 1) * signif_to_zero(abs_max, 1)){
  colors <- c(scales::muted("blue"), "white", "white", scales::muted("red"))
  values <- c(0, 0.5 - mid_width/2, 0.5 + mid_width/2, 1)
  ggplot2::scale_color_gradientn(oob = oob, limits = limits, breaks = breaks, colors = colors, values = values, ...)
}