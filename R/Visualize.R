#' Visualize Variant Heat Map
#'
#' Sometimes, it may be beneficial to know approximately which
#' genomic regions have the most detected variants relative to
#' reference sequence. This function generates a simple figure
#' depicting the positions on the horizontal axis and the number
#' of variants on the vertical axis.
#'
#' @param heatmap The data frame containing the variant distribution of the
#' various nucleotides at each positions.
#'
#' @param nucbase The nucleotides to graph. Can be any combination of
#' "A", "T", "C", or "G".
#'
#' @examples
#' plot <- visVariantHeatMap(heatmap)
#' View(plot)
#'
#' @export
#' @import ggplot2
visVariantHeatMap <- function(heatmap, nucbases) {
  refColNames <- base::paste("ref", nucbases, sep = "")
  altColNames <- base::paste("alt", nucbases, sep = "")
  distribution <- data.frame(
    refs = rowSums(heatmap[refColNames]),
    alts = rowSums(heatmap[altColNames])
  )

  ggplot2::ggplot(data = distribution) +
    ggplot2::geom_line(aes(x = as.numeric(row.names(distribution)), y = refs), colour = "blue", alpha = 0.6) +
    ggplot2::geom_point(aes(x = as.numeric(row.names(distribution)), y = refs), colour = "blue", alpha = 0.6) +
    ggplot2::geom_line(aes(x = as.numeric(row.names(distribution)), y = alts), colour = "red", alpha = 0.6) +
    ggplot2::geom_point(aes(x = as.numeric(row.names(distribution)), y = alts), colour = "red", alpha = 0.6) +
    ggplot2::labs(x = "Position", y = "Frequency")
  return(plot)
}
