library(PhenoGenRLib)

test_that("loaded VCFs do have metadata columns", {
  variants <- PhenoGenRLib::linkVariantsWithMetadata(
    metadata = system.file("extdata", "huntingtons_datasheet_shortened.csv", package = "PhenoGenRLib"),
    vcfDir = system.file("extdata", package = "PhenoGenRLib"),
    vcfColName = "vcfs"
  )
  expect_contains(colnames(variants), c("dummy_pheno", "chromosome", "vcfs"))
})

test_that("loaded VCFs have correct data in metadata columns", {
  variants <- PhenoGenRLib::linkVariantsWithMetadata(
    metadata = system.file("extdata", "huntingtons_datasheet_shortened.csv", package = "PhenoGenRLib"),
    vcfDir = system.file("extdata", package = "PhenoGenRLib"),
    vcfColName = "vcfs"
  )
  expect_equal(variants[1,"dummy_pheno"], "A")
})

test_that("rsid mappings return two distinct values", {
  mapped <- PhenoGenRLib::mapRsidsForVariants("chromosome", variants = UnmappedVariants, offset = 3074680)
  expect_length(mapped, 2)
})

test_that("rsid mappings result in correct data columns", {
  mapped <- PhenoGenRLib::mapRsidsForVariants("chromosome", variants = UnmappedVariants, offset = 3074680)
  mappedNvs <- mapped$nvs
  rsids <- mapped$rsids

  expect_contains(
    colnames(mappedNvs),
    c("refsnp_id", "chrom_start", "chrom_end")
  )
  expect_contains(
    colnames(rsids),
    c(
      "refsnp_id",
      "allele",
      "minor_allele",
      "minor_allele_count",
      "minor_allele_freq",
      "chr_name",
      "chrom_start",
      "chrom_end")
    )
})

test_that("rsid mappings have at least 1 mapping", {
  mapped <- PhenoGenRLib::mapRsidsForVariants("chromosome", variants = UnmappedVariants, offset = 3074680)
  mappedNvs <- mapped$nvs
  expect_gte(length(mappedNvs[which(!is.na(mappedNvs$refsnp_id)),]), 1)
})