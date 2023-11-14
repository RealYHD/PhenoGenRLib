library(stringr)
library(tibble)

#' Variance Case-Control Analysis where the Case and Control are User Provided
#' 
#' Similar to varianceCCAnalysisEnsembl, however, instead of using Ensembl,
#' the name of a column from the metadata used to load the VCFs is taken.
#' The function then finds all possible categorical values for the column
#' and builts pairwise tables based off these.
#' 
#' @param variants The data.frame containing the variant information.
#' @param totalCaseSamples The total number of case samples that goes into
#' the variants data.frame. In other words, this should equal the number of rows
#' found in the metadata file used to load the VCFs as each row in that file
#' should represent one case.
#' @param phenotypeName The name of the phenotype to categorize upon. E.g., if
#' each VCF had trinary column named "Type" and the values in the column are of
#' A, B, or C, then the categorical values this function would discover are A,
#' B, and C and therefore, would require the breaking down of the tables pairwise.
#' @param useChi Instead of using Fisher's exact test, the chi-square test will be employed.
#' Default is FALSE.
#' 
#' @return A data.frame containing all the p-values and the shape of the 2-way table
#' employed for each test.
#' 
#' @export 
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


#' Variance Case-Control Analysis where the Ensembl Database is the Control
#' 
#' This function makes it easy to test for potential associations between the
#' studied case, and the average population. To elaborate, this function takes
#' a data.frame of all the variants, and the rsIDs found to be linked with said
#' variants, and uses the pre-existing allele frequency information from the rsID
#' data.frame (fetched from Ensembl). Then, it moves over each mutation with population
#' information in public databases and runs a Fisher's exact test. The Fisher's exact
#' test is beneficial for low sample counts as it is an exact measure of correlation in
#' contrast to various other tests. It is more computational intensive however. Hence,
#' this function also allows the use of the chi-square test. The table formed is between
#' whether or not the position is mutated, and whether or not the mutation occurred in a
#' healthy individual (the control group), or the case study individual.
#' 
#' @param variants The data.frame containing all the variants, variant information, and rsIDs.
#' @param rsids The data.frame containing all the extra information regarding the known
#' variants.
#' @param totalCaseSamples The total number of case samples. This should be the number of rows
#' in the metadata file as each row represents one case, hence, one VCF.
#' @param useChi Instead of using Fisher's exact test, the chi-square test will be employed.
#' Default is FALSE.
#' 
#' @return A data.frame containing all the p-values and the shape of the 2-way table
#' employed for each test.
#'  
#' @export 
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
