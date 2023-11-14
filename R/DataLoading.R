library(readr)
library(tibble)
library(bedr)
library(biomaRt)

linkVariantsWithMetadata <- function(metadataFile, vcfDir, vcfColName) {
    metadata <- readr::read_csv(metadataFile)
    variantData <- NULL
    for(i in 1:nrow(metadata)) {
        metadataRow <- metadata[i,]
        vcfFile <- normalizePath(file.path(vcfDir, metadataRow[vcfColName]))
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
    return(c(mappedVariants, rsids))
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
