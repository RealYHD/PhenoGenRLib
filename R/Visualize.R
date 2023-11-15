library(ggplot2)


#' Visualize Variant Distribution
#' 
#' Sometimes, it may be beneficial to know approximately which
#' genomic regions have the most detected variants relative to 
#' reference sequence. This function generates a simple figure
#' depicting the positions on the horizontal axis and the number 
#' of variants on the vertical axis.
#' 
#' @param variants The data.frame containing all the variant
#' information.
#' 
#' @export 
visVariantDistribution <- function(variants) {
    greatestNVPos <- base::max(
        variants$POS + base::apply(variants["REF"], MARGIN = 1, base::nchar)
    )
    smallestNVPos <- min(variants$POS) 
    caseOccurrences <- base::vector(
        mode = "integer",
        length = greatestNVPos - smallestNVPos)
    for (row in i:base::nrow(variants)) {
        variant <- variants[row,]
        start <- variant[["POS"]] - smallestNVPos
        end <- start + base::nchar(variant[["REF"]]) - 1
        caseOccurrences[start:end] <- caseOccurrences[start:end] + 1
    }
    positions <- smallestNVPos:(greatestNVPos - 1)
    distribution <- base::data.frame(positions, caseOccurrences)
    
    ggplot2::ggplot(distribution, aes(x = positions, y = caseOccurrences)) +
        geom_bar(stat = "identity", fill = "red", width = 1) +
        geom_point(size = 2, colour = "red") +
        labs(x = "Position", y = "Occurrences", title = "Distribution of Occurrences") +
        theme(aspect.ratio = 1/2)
    
    return(NULL)
}