library(ggplot2)

variantDistributionEnsembl <- function(variants) {
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
    
    ggplot(distribution, aes(x = positions, y = caseOccurrences)) +
        geom_bar(stat = "identity", fill = "red", width = 1) +
        geom_point(size = 2, colour = "red") +
        labs(x = "Position", y = "Occurrences", title = "Distribution of Occurrences") +
        theme(aspect.ratio = 1/2)
}