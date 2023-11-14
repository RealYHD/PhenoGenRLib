library(stringr)
library(tibble)

varianceCCAnalysisEnsembl <- function(variants, rsids, totalCaseSamples, useChi = FALSE) {
    rsidsWithFreq <- rsids[!base::is.na(rsids["minor_allele_count"]),]
    uniqueRsids <- base::intersect(base::unique(variants[["refsnp_id"]]), base::unique(rsidsWithFreq[["refsnp_id"]]))
    results <- NULL
    pivotTables <- list()
    for (rsid in uniqueRsids) {
        specificVariants <- variants[which(variants["refsnp_id"]==rsid),]
        refAllele <- specificVariants[[1, "REF"]]
        allAlleles <- base::unique(c(
            refAllele,
            specificVariants[["ALT"]],
            unlist(stringr::str_split(rsids[rsids$refsnp_id==rsid,]$allele, "/"))
        ))
        numRef <- totalCaseSamples - nrow(specificVariants)
        caseFreqs <- base::lapply(allAlleles,
            function(x) {
                if (x == refAllele) return(numRef)
                return(length(which(specificVariants$ALT==x)))
            }
        )
        specificRsidStats <- rsidsWithFreq[rsidsWithFreq$refsnp_id==rsid,]
        totalControlSamples <- base::ceiling(specificRsidStats$minor_allele_count / specificRsidStats$minor_allele_freq)
        controlFreqs <- base::lapply(allAlleles, function(x) {
            if (x == specificRsidStats[["minor_allele"]]) {
                return(specificRsidStats[["minor_allele_count"]])
            } else {
                return(totalControlSamples - specificRsidStats[["minor_allele_count"]])
            }
        })
        freqTable <- data.frame(
            case = base::unlist(caseFreqs),
            control = base::unlist(controlFreqs),
            row.names = allAlleles
        )
        pivotTables <- base::append(pivotTables, freqTable)
        result <- if (!useChi) stats::fisher.test(base::as.matrix(freqTable)) else stats::chisq.test(base::as.matrix(freqTable))
        if (base::is.null(results)) {
            results <- data.frame(
                rsid = c(rsid),
                pVal = c(result$p.value),
                method = c(result$method),
                name = c(result$data.name)
            )
        } else {
            results <- base::rbind(
                results,
                c(
                    rsid,
                    result$p.value,
                    result$method,
                    result$data.name
                )
            )
        }
    }
    return(c(results, pivotTables))
}
