library(stringr)
library(tibble)

varianceCCAnalysisPheno <- function(variants, totalCaseSamples, phenotypeName, useChi = FALSE) {
    phenotypes <- base::unique(variants[[phenotypeName]])
    uniqueRsids <- base::unique(variants[["refsnp_id"]])
    results <- NULL
    for (rsid in uniqueRsids) {
        rsidVariants <- variants[which(variants["refsnp_id"]==rsid),]
        refAllele <- rsidVariants[[1, "REF"]]
        allAlleles <- base::unique(c(
            refAllele,
            rsidVariants[["ALT"]]
        ))
        freqMatrix <- NULL
        for (phenotype in phenotypes) {
            phenotypeVariants <- rsidVariants[rsidVariants[phenotypeName]==phenotype,]
            freqs <- base::lapply(allAlleles, 
                function(x) {
                    if (x == refAllele) return(totalCaseSamples - nrow(phenotypeVariants))
                    return(length(which(phenotypeVariants$ALT==x)))
                }
            )
            if (base::is.null(freqMatrix)) {
                freqMatrix <- freqs
            } else {
                freqMatrix <- base::cbind(freqMatrix, freqs)
            }
        }
        row.names(freqMatrix) <- allAlleles
        colnames(freqMatrix) <- phenotypes
        results <- multipleAssociationTests(
            generate2wayFromMxN(freqMatrix),
            useChi = useChi,
            groupName = rsid
        )
    }
    results <- if (base::is.nan(results)) result else base::rbind(results, result)
}

varianceCCAnalysisEnsembl <- function(variants, rsids, totalCaseSamples, useChi = FALSE) {
    rsidsWithFreq <- rsids[!base::is.na(rsids["minor_allele_count"]),]
    uniqueRsids <- base::intersect(base::unique(variants[["refsnp_id"]]), base::unique(rsidsWithFreq[["refsnp_id"]]))
    results <- NULL
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
        result <- multipleAssociationTests(
            generate2wayFromMxN(freqTable),
            useChi = useChi,
            groupName = rsid
        )
        results <- if (base::is.nan(results)) result else base::rbind(results, result)
    }
    return(results)
}

generate2wayFromMxN <- function(matrix) {
    matrix <- base::as.matrix(matrix)
    columnCombs <- utils::combn(base::ncol(matrix), m = 2)
    rowCombs <- utils::combn(base::nrow(matrix), m = 2)
    results <- list()
    
    for (col in 1:base::ncol(columnCombs)) {
        colComb <- columnCombs[,col]
        for (row in 1:base::ncol(rowCombs)) {
            rowComb <- rowCombs[,row]
            pairwise <- matrix[rowComb, colComb]
            results <- c(results, list(pairwise))
        }
    }
    return(results)
}

multipleAssociationTests <- function(matrices, useChi = FALSE, groupName = NULL) {
    results = NULL
    for (mat in matrices) {
        mat <- base::as.matrix(mat)
        # print(mat)
        result <- base::tryCatch(
            {
                result <- if (!useChi) stats::fisher.test(mat) else stats::chisq.test(mat)
                testResultSummary <- c(
                    if (base::is.null(groupName)) "No name" else groupName,
                    list(base::colnames(mat)),
                    list(base::row.names(mat)),
                    result$data.name,
                    result$p.value,
                    result$method
                )
                return(testResultSummary)
            },
            error = function(e) {
                testResultSummary <- c(
                    if (base::is.null(groupName)) "No name" else groupName,
                    list(base::colnames(mat)),
                    list(base::row.names(mat)),
                    base::toString(e),
                    NA,
                    "NA"
                )
                return(testResultSummary)
            }
        )
        results <- (if (base::is.null(results)) result else rbind(results, result))
    }
    base::colnames(results) <- c("Test Group", "Groups", "Categories", "Call", "P-Value", "Method")
    base::row.names(results) <- NULL
    return(results)
}
