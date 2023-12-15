#' Variance Case-Control Analysis where the Case and Control are User Provided
#'
#' Similar to varianceCCAnalysisEnsembl, however, instead of using Ensembl,
#' the name of a column from the metadata used to load the VCFs is taken.
#' The function then finds all possible categorical values for the column
#' and built pairwise tables based off these.
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
  if (!base::is.data.frame(variants)) {
    base::stop("Variants should be a data frame.")
  }
  if (!base::is.numeric(totalCaseSamples)) {
    base::stop("The total number of case samples should be numeric.")
  }
  if (! phenotypeName %in% variants) {
    base::stop("The name of the phenotype must be a column in the variants data frame.")
  }
  if (!base::is.logical(useChi)) {
    base::stop("Parameter determining whether or not to use Chi-Square test must be logical.")
  }

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
      generate2WayFromMxN(freqMatrix),
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
#' employed for each test. The structure is that of
#' @seealso PhenoGenRLib::MultipleAssociationTests
#'
#' @export
#'
#' @importFrom dplyr distinct
#' @importFrom stringr str_split
varianceCCAnalysisEnsembl <- function(variants, rsids, totalCaseSamples, useChi = FALSE) {
  if (!base::is.data.frame(variants)) {
    base::stop("Variants should be a data frame.")
  }
  if (!base::is.data.frame(rsids)) {
    base::stop("RSIDs should be a data frame.")
  }
  if (!base::is.numeric(totalCaseSamples)) {
    base::stop("The total number of case samples should be numeric.")
  }
  if (!base::is.logical(useChi)) {
    base::stop("Parameter determining whether or not to use Chi-Square test must be logical.")
  }


  rsidsWithFreq <- dplyr::distinct(rsids[!base::is.na(rsids["minor_allele_count"]),]) # Only keep variants with frequency counts
  uniqueRsids <- base::intersect( # Get all the unique variants that exist between the two groupings
    base::unique(variants[["refsnp_id"]]),
    base::unique(rsidsWithFreq[["refsnp_id"]])
  )
  results <- data.frame(matrix(ncol = 5, nrow = 0))

  for (rsid in uniqueRsids) { # for every unique variant
    specificVariants <- variants[which(variants["refsnp_id"]==rsid),] # get all same variants
    refAllele <- specificVariants[[1, "REF"]] # Fetch the reference allele for the set of variants
    controlAlleles <- unlist(stringr::str_split(rsidsWithFreq[rsidsWithFreq$refsnp_id==rsid,]$allele, "/"))
    allAlleles <- base::unique(c( # Find all unique Allele variants
      refAllele,
      specificVariants[["ALT"]],
      controlAlleles
    ))
    numRef <- totalCaseSamples - nrow(specificVariants) # Calculate number of cases with reference allele
    caseFreqs <- base::lapply(allAlleles, # Count number of alleles for each alternate allele
                              function(x) {
                                if (x == refAllele) return(numRef) # If the requested allele to count is ref, return ref count
                                return(length(which(specificVariants$ALT==x))) # Otherwise, count the alternate alleles that equal X
                              }
    )
    specificRsidStats <- rsidsWithFreq[rsidsWithFreq$refsnp_id==rsid,] # Get all variants with the same IDs that have frequency data
    totalControlSamples <- base::ceiling(specificRsidStats$minor_allele_count / specificRsidStats$minor_allele_freq)[[1]] # Should be the same across all alleles
    controlFreqs <- base::lapply(allAlleles, function(x) {
      for (mai in 1:nrow(specificRsidStats)) {
        ma <- specificRsidStats[mai, c("minor_allele", "minor_allele_count", "minor_allele_freq")]
        if (x == ma[["minor_allele"]]) {
          count <- ma[["minor_allele_count"]]
          return(count)
        }
      }
      if (x %in% controlAlleles) { # IF not a novel allele found in case dataset
        return(totalControlSamples - sum(specificRsidStats[["minor_allele_count"]])) # Then it must be the major allele
      }
      return(0) # Assume no catalogued SNVs of this variety
    })
    freqTable <- base::data.frame(
      case = base::unlist(caseFreqs),
      control = base::unlist(controlFreqs)
    )
    row.names(freqTable) <- allAlleles

    result <- multipleAssociationTests(
      generate2WayFromMxN(freqTable),
      useChi = useChi,
      groupName = rsid
    )
    results <- base::rbind(results, result, stringsAsFactors = FALSE)
  }
  return(results)
}

#' Generates Two-Way Tables from MxN Data Frame
#'
#' This function takes a MxN data and solves for all pairwise combinations
#' of rows and columns and creates lists of the 2x2 matrices that result
#' from these combinations.
#'
#' @param mxn The MxN data frame to solve for pairwise permutations of.
#'
#' @return Returns a list of 2x2 two-way tables in the form of data frames
#' while retaining the row and column names.
#'
#' @export
generate2WayFromMxN <- function(mxn) {
  if (!base::is.data.frame(mxn)) {
    base::stop("The input should be a data.frame.")
  }

  columnCombs <- utils::combn(base::ncol(mxn), m = 2)
  rowCombs <- utils::combn(base::nrow(mxn), m = 2)
  results <- list()

  for (col in 1:base::ncol(columnCombs)) {
    colComb <- columnCombs[,col]
    for (row in 1:base::ncol(rowCombs)) {
      rowComb <- rowCombs[,row]
      pairwise <- mxn[rowComb, colComb]
      results[[length(results) + 1]] <- data.frame(pairwise)
    }
  }
  return(results)
}

#' Performs Multiple Association Tests on Arrays of Matrices
#'
#' For a list of 2x2 data frames, this function calculates the association
#' between the the individual tables categorical variables by application of
#' either Fisher's exact test or Chi-Square test.
#'
#' @param mxns A list of 2x2 data frames.
#' @param useChi Whether or not to use Chi-Square tests instead of Fisher's exact test
#' @param groupName The name of this test group. Annotates the tests ran with this name
#' in the results.
#'
#' @return Returns a data.frame with columns describing the 2-way table that was used for
#' the test:
#' \itemize{
#'   \item test_group - The group the test belongs too. Often used to group based on position.
#'   \item groups - The names of the columns of the 2-way table tested on.
#'   \item categories - The names of rows of the 2-way table tested on.
#'   \item p_value - The resulting p-value
#'   \item additional - Any additional information regarding the automated test.
#' }
#'
#' @examples
#' testResults <- multipleAssociationTests(
#'   mxns = list(data.frame(A = c(1, 2, 3), B = c(1, 2, 3), row.names = c("C", "D"))),
#'   groupName = "First Group"
#' )
#' View(testResults)
#'
#' @export
multipleAssociationTests <- function(mxns, useChi = FALSE, groupName = "Untitled") {
  if (!base::is.list(mxns)) {
    base::stop("The input should be a list of data frames.")
  }
  if (!base::is.logical(useChi)) {
    base::stop("Parameter determining whether or not to use Chi-Square test must be logical.")
  }

  results <- base::data.frame(base::matrix(ncol = 5, nrow = 0))
  for (mxn in mxns) {
    mat <- base::as.matrix(mxn)

    base::tryCatch(
      {
        result <- if (!useChi) stats::fisher.test(mat) else stats::chisq.test(mat)
        testResultSummary <- base::data.frame(
          test_group = groupName,
          group = base::toString(base::colnames(mxn)),
          categories = base::toString(base::row.names(mxn)),
          p_value = base::as.numeric(result$p.value),
          additional = "N/A"
        )
        # Append result as a row to results
        results <- rbind(results, testResultSummary, stringsAsFactors = FALSE)
      },
      error = function(e) {
        testResultSummary <- base::data.frame(
          test_group = groupName,
          group = base::toString(base::colnames(mxn)),
          categories = base::toString(base::colnames(mxn)),
          p_value = NA,
          additional = base::toString(e)
        )
        results <- base::rbind(results, testResultSummary, stringsAsFactors = FALSE)
      }
    )
  }
  base::colnames(results) <- c(
    "test_group",
    "groups",
    "categories",
    "p_value",
    "additional"
  )
  return(results)
}


#' Generates a vector of NV Occurrences
#'
#' Takes a data frame of known nucleotides variants and counts their occurrences
#' on different positions. Groups ocurrences based on the variant allele.
#'
generatePositionHeatmap <- function(variants) {
  greatestNVPos <- base::max(
    variants$POS + base::max(
      base::apply(variants["REF"], MARGIN = 1, base::nchar), # Accounts for MNVs
      base::apply(variants["ALT"], MARGIN = 1, base::nchar)
    )
  ) # Find the largest variant position
  smallestNVPos <- min(variants$POS) # find the smallest variant position
  varianceRegionLength <- greatestNVPos - smallestNVPos + 1
  occurrences <- data.frame(
    altA = base::vector(
      mode = "integer",
      length = varianceRegionLength
    ),
    altT = base::vector(
      mode = "integer",
      length = varianceRegionLength
    ),
    altC = base::vector(
      mode = "integer",
      length = varianceRegionLength
    ),
    altG = base::vector(
      mode = "integer",
      length = varianceRegionLength
    ),
    refA = base::vector(
      mode = "integer",
      length = varianceRegionLength
    ),
    refT = base::vector(
      mode = "integer",
      length = varianceRegionLength
    ),
    refC = base::vector(
      mode = "integer",
      length = varianceRegionLength
    ),
    refG = base::vector(
      mode = "integer",
      length = varianceRegionLength
    ),
    row.names = smallestNVPos:greatestNVPos
  )
  for (row in 1:base::nrow(variants)) {
    variant <- variants[row,]
    start <- variant[["POS"]] - smallestNVPos
    end <- start + base::nchar(variant[["REF"]]) - 1

    for (relPos in 1:base::nchar(variant[["REF"]])) {
      baseType <- base::substring(variant[["REF"]], relPos, relPos)
      pos <- start + relPos - 1 # because relPos starts at 1.
      colName <- paste("ref", baseType, sep = "")
      occurrences[pos, colName] <- occurrences[pos, colName] + 1
    }

    # Account for multiple variants
    for (variantSeq in base::unlist(stringr::str_split(variant[["ALT"]], pattern = ","))) {
      for (relPos in 1:base::nchar(variantSeq)) {
        baseType <- base::substring(variantSeq, relPos, relPos)
        pos <- start + relPos - 1
        colName <- paste("alt", baseType, sep = "")
        occurrences[pos, colName] <- occurrences[pos, colName] + 1
      }
    }
  }
  return(occurrences)
}

