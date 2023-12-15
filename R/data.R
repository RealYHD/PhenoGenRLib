#' Demo Huntingtons Patient Variants with Metadata and Mapped to Known RSIDs
#'
#' Similar to `huntingtonVariants`. This dataset is derived from the previous
#' by additionally calling the RSID mapping function on `huntingtonVariants`.
#'
#' This dataset is a list of two different data frames, one of which under the
#' the name `mappedHuntingtonVariants$nvs` is largely the same with the
#' `huntingtonVariants` data frame, however, contains extra information on the
#' associated RSIDs. The second data frame is a list of all RSIDs pulled and
#' contains more information for each of the RSIDs. The latter can be found
#' under the name `mappedHuntingtonVariants$rsids`.
#'
#' @seealso [PhenoGenRLib::mapRsidsForVariants()]
#'
#' @source NCBI BioProject PRJNA299309
#'
#'
"mappedHuntingtonsVariants"

#' Demo Huntingtons Patient Variants Linked with Metadata
#'
#' FASTQ files obtained from BioProject with accession PRJNA299309 was loaded
#' into Bowtie2, for alignment, and Freebayes for variant calling. Produced
#' Variant Call Format (VCF) files were loaded with the linking function
#' built-in to the library.
#'
#' This is a data frame produced from reading multiple VCFs and linking their
#' associated metadata. On top of the the standard tags typically seen in a VCF
#' this file also contains the filename itself, an additional chromosome tag
#' with human readable numbers indicating the chromosome of origin for the
#' sample, and a `dummy_pheno`, representing whether or not this patient has
#' phenotype A, or B.
#'
#' @seealso [PhenoGenRLib::linkVariantsWithMetadata()]
#'
#' @source NCBI BioProject PRJNA299309
#'
#'
"huntingtonsVariants"

#' Demo Huntingtons Patient Variants Heatmap
#'
#' A Dataframe containing 8 columns, two columns representing the reference,
#' and the alternate alleles of A, T, C and G. Heatmap generated from the
#' data set `huntingtonsVariants` by use of the heatmap function.
#'
#' @seealso [PhenoGenRLib::generatePositionHeatmap()]
#'
#' @source NCBI BioProject PRJNA299309
#'
"huntingtonsVariantsHeatmap"
