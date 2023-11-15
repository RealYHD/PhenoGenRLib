library(readr)
library(tibble)
library(bedr)
library(biomaRt)


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
#' @param metadataFile A path to the metadata file in the CSV format. This
#' file must at the very least, contain one column which lists the various
#' VCF's filenames. For each filename, associated metadata should appear
#' on the same row.
#' @param vcfDir A path to the directory containing the VCFs. The previously
#' mentioned filenames will be appended in a smart manner to from the full path.
#' @param vcfColName The name of the column containing the VCF files.
#' 
#' @return Returns a data.frame like object where each row is a variant and
#' the associated metadata is attached column-wise to the right. 
#' 
#' @export 
linkVariantsWithMetadata <- function(metadataFile, vcfDir, vcfColName) {
    metadata <- readr::read_csv(metadataFile)
    variantData <- NULL
    for(i in 1:nrow(metadata)) {
        metadataRow <- metadata[i,]
        vcfFile <- base::normalizePath(base::file.path(vcfDir, metadataRow[vcfColName]))
        vcfs <- bedr::read.vcf(vcfFile, split.info= TRUE, split.samples = TRUE)
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
#' 
#' @return Returns variant data.frame containing all pre-existing information (including metadata and
#' other built-in VCF annotations), as well as the associated rsIDs, chromosomal and known chromosomal
#' start and end positions. Also returns a data.frame of all the rsIDs and associated information such
#' allele frequencies in the databases. The two values are returned in a data.frame where `nvs` is the 
#' accessor for the variants, and `rsids` is the accessor for the rsID information. 
#' 
#' @export 
mapRsidsForVariants <- function(chromCol, variants, offset = 0, hostGenVersion = 38, batchSize = 100) {
    nvCoords <- coordinatesFromVariants(variants = variants, offset = offset)
    mappingData <- variants
    mappingData["chrom_start"] <- nvCoords$POS
    mappingData["chrom_end"] <- nvCoords$END

    batchIndex <- 1
    rsids <- NULL
    while(batchIndex <= nrow(nvCoords)) {
        endIndex <- min(batchIndex + batchSize - 1, nrow(nvCoords))
        nvCoordBatch <- nvCoords[batchIndex:endIndex,]
        rsidBatch <- mapRsidsForVariantPositions(coordinates = nvCoordBatch, hostGenVersion = hostGenVersion)
        rsids <- if (base::is.null(rsids)) rsidBatch else base::rbind(rsids, rsidBatch)
        batchIndex <- endIndex + 1
    }

    mappedVariants <- base::merge(mappingData, rsids[c("chrom_start", "chrom_end", "refsnp_id")], by = c("chrom_start", "chrom_end"), all.x = TRUE)
    return(base::data.frame(nvs = mappedVariants, rsids = rsids))
}

coordinatesFromVariants <- function(variants, offset = 0) {
    queryData <- variants[c(chromCol, "POS")]
    queryData$POS <- queryData$POS + offset
    queryData["END"] <- queryData$POS + apply(variants["REF"], MARGIN = 1, nchar) - 1
    return(queryData)
}

mapRsidsForVariantPositions <- function(coordinates, hostGenVersion = 38) {
    coords <- apply(coordinates, 1, paste, collapse = ":")
    snpMart <- biomaRt::useEnsembl(biomart = "snps", dataset = "hsapiens_snp", host = (if (hostGenVersion == 38) "https://www.ensembl.org" else "https://grch37.ensembl.org"))
    rsids <- biomaRt::getBM(
        attributes = c("refsnp_id", "allele", "minor_allele", "minor_allele_count", "minor_allele_freq", "chr_name", "chrom_start", "chrom_end"),
        filters = c("chromosomal_region"),
        values = coords,
        mart = snpMart,
        uniqueRows = TRUE
    )
    return(rsids)
}
