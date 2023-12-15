#' Link Existing VCF files with Metadata
#'
#' This function's purpose is to join a series of VCFs together and
#' associate sample, or other metadata with each of the variants inside
#' each of the previously mentioned VCFs. To do this, one column in the
#' metadata CSV file must contain the file names of each VCF, and on the
#' same row, contain the associated metadata with that VCF. Then,
#' it looks inside `vcfDir` for the VCF filenames as they appear in the
#' CSV.
#'
#' @param metadata A path to the metadata file in the CSV format or a data frame
#' . This file must at the very least, contain one column which lists the
#' various VCF's filenames. For each filename, associated metadata should appear
#' on the same row.
#' @param vcfDir A path to the directory containing the VCFs. The previously
#' mentioned filenames will be appended in a smart manner to from the full path.
#' @param vcfColName The name of the column containing the VCF files.
#' @param progress A function to call each iteration to track progress. The
#' function should accept two numeric parameters, one representing the current
#' iteration, and the second representing the total number of iterations.
#'
#' @return Returns a data.frame like object where each row is a variant and
#' the associated metadata is attached column-wise to the right.
#'
#' @examples
#'
#' variants <- PhenoGenRLib::linkVariantsWithMetadata(
#'   metadata = system.file("extdata/huntingtons_datasheet_shortened.csv",
#'     package = "PhenoGenRLib"),
#'   vcfDir = system.file("extdata/", package = "PhenoGenRLib"),
#'   vcfColName = "vcfs"
#' )
#' # Use View(variants) to see results!
#'
#' @export
#' @importFrom readr read_csv
#' @importFrom bedr read.vcf
#' @importFrom stringr str_trim
linkVariantsWithMetadata <- function(metadata, vcfColName, vcfDir = NULL, progress = NULL) {
  if (is.data.frame(metadata)) {
    # Nothing
  } else if (is.character(metadata)) {
    metadata <- readr::read_csv(metadata, show_col_types = FALSE)
  } else {
    base::stop("Metadata must be either data frame or path to CSV.")
  }

  variantData <- NULL
  for(i in 1:nrow(metadata)) {
    if (is.null(progress)) {
      # do nothing
    } else {
      progress(i, nrow(metadata))
    }
    metadataRow <- metadata[i,]
    vcfFile <- base::unlist(metadataRow[[vcfColName]])
    if (is.null(vcfDir) || base::nchar(vcfDir) == 0) {
      vcfFile <- base::normalizePath(vcfFile)
    } else {
      vcfFile <- base::normalizePath(base::file.path(vcfDir, vcfFile))
    }
    if (file.exists(vcfFile)) {
      # Nothing
    } else {
      stop("Could not find some of the VCFs.")
    }
    vcfs <- bedr::read.vcf(vcfFile, split.info = TRUE, split.samples = TRUE)
    header <- vcfs[[1]]
    vcfs <- vcfs[2:length(vcfs)]
    for (vcf in vcfs) {
      # Build the extension in terms of columns
      metadataExt <- base::do.call(rbind, replicate(nrow(vcf), metadataRow, simplify = FALSE))
      extendedVcf <- cbind(vcf, metadataExt)
      if (is.null(variantData)) {
        variantData <- extendedVcf
      } else {
        variantData <- rbind(variantData, extendedVcf)
      }
    }
  }
  variantData["local_id"] <- apply(
    variantData[c("REF", "POS")],
    MARGIN = 1,
    FUN = function (data) {
      ref <- data[["REF"]]
      pos <- data[["POS"]]
      return(paste(stringr::str_trim(ref), stringr::str_trim(pos), sep = ""))
  })
  return(variantData)
}

#' Maps RefSeq IDs to the variants
#'
#' This function takes a data.frame like object containing the variant
#' information and attached metadata. It then searches the ensembl database
#' to find rsIDs for the well-known variants that appear. The rsIDs are then
#' added column-wise to each of the variants in the original data.frame. This
#' updated data.frame is what is returned.
#'
#' @param chromCol The column containing the name of the chromosome. It may be
#' beneficial to use the metadata file to add this as often times, the VCF's chromosome
#' is a ID rather than just a number.
#' @param variants The variant data.frame to map rsIDs for.
#' @param offset VCFs report positions relative to the reference sequence. Of course,
#' if this reference sequence is not the entire chromosome, the reported coordinate will
#' not match with any known rsIDs. Therefore, an offset may need to be introduced. E.g.,
#' if the reference sequence started at 3074681 (ex. NC_000004.12:3074681-3243960), we may
#' may choose the offset to be 3074680. Default is 0 (no offset).
#' @param hostGenVersion The Host genome version. The two supported options are 38 (GRCh38),
#' and 37 (GRCh37). Default is 38.
#' @param batchSize The number coordinates to send to Ensembl servers at once. Smaller batches
#' may reduce chances of timing out. Default is 100.
#' @param progress A function to call each iteration to track progress. The
#' function should accept two numeric parameters, one representing the current
#' iteration, and the second representing the total number of iterations.
#'
#' @return Returns variant data.frame containing all pre-existing information (including metadata and
#' other built-in VCF annotations), as well as the associated rsIDs, chromosomal and known chromosomal
#' start and end positions. Also returns a data.frame of all the rsIDs and associated information such
#' allele frequencies in the databases. The two values are returned in a data.frame where `nvs` is the
#' accessor for the variants, and `rsids` is the accessor for the rsID information.
#'
#' @examples
#'
#' mappedVariants <- mapRsidsForVariants(
#' chromCol = "chromosome",
#' variants = huntingtonsVariants,
#' offset = 3074680
#' )
#' # Use View(mappedVariants$nvs) to see the variants and their mappings
#' # Use View(mappedVariants$rsids) to see just the RSID information
#' @export
#'
mapRsidsForVariants <- function(chromCol, variants, offset = 0, hostGenVersion = 38, batchSize = 100, progress = NULL) {
  nvCoords <- coordinatesFromVariants(
    variants = variants,
    chromCol = chromCol,
    offset = offset
  )
  mappingData <- variants
  mappingData["chrom_start"] <- nvCoords$POS
  mappingData["chrom_end"] <- nvCoords$END

  batchIndex <- 1
  rsids <- NULL
  while(batchIndex <= base::nrow(nvCoords)) {
    if (is.null(progress)) {
      # Do nothing
    } else {
      progress(batchIndex, nrow(nvCoords))
    }
    endIndex <- base::min(batchIndex + batchSize - 1, base::nrow(nvCoords))
    nvCoordBatch <- nvCoords[batchIndex:endIndex,]
    rsidBatch <- mapRsidsForVariantPositions(coordinates = nvCoordBatch, hostGenVersion = hostGenVersion)
    rsids <- if (base::is.null(rsids)) rsidBatch else base::rbind(rsids, rsidBatch)
    batchIndex <- endIndex + 1
  }

  mappedVariants <- base::merge(
    mappingData,
    rsids[base::c("chrom_start", "chrom_end", "refsnp_id")],
    by = base::c("chrom_start", "chrom_end"), all.x = TRUE
  )
  return(base::list(nvs = mappedVariants, rsids = rsids))
}

#' Convert From Coordinates to Variants
#'
#' Helper function for rapidly performing vectorized variant positional calculations.
#'
coordinatesFromVariants <- function(variants, chromCol, offset = 0) {
  queryData <- variants[base::c(chromCol, "POS")]
  queryData$POS <- queryData$POS + offset
  queryData["END"] <- queryData$POS + base::apply(variants["REF"], MARGIN = 1, nchar) - 1
  return(queryData)
}

#' Maps RSIDs based on Variant Positions
#'
#' Helper function for fetching RSIDs based off set of specific coordinates
#'
#' @importFrom biomaRt useEnsembl getBM
mapRsidsForVariantPositions <- function(coordinates, hostGenVersion = 38) {
  coords <- base::apply(coordinates, 1, paste, collapse = ":")
  snpMart <- biomaRt::useEnsembl(
    biomart = "snps",
    dataset = "hsapiens_snp",
    host = (if (hostGenVersion == 38) "https://www.ensembl.org" else "https://grch37.ensembl.org"))
  rsids <- biomaRt::getBM(
    attributes = base::c(
      "refsnp_id",
      "allele",
      "minor_allele",
      "minor_allele_count",
      "minor_allele_freq",
      "chr_name",
      "chrom_start",
      "chrom_end"
    ),
    filters = base::c("chromosomal_region"),
    values = coords,
    mart = snpMart,
    uniqueRows = TRUE
  )
  return(rsids)
}

# [END]
