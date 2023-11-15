library(PhenoGenRLib)

test_that("case control with ensembl database as control produces expected table", {
  results <- PhenoGenRLib::varianceCCAnalysisEnsembl(
    MappedVariants,
    RsIDs,
    totalCaseSamples = 8
  )
})
