library(PhenoGenRLib)

testthat::test_that("pairwise table generation produces correct permutations of 3x3", {
  results <- PhenoGenRLib:::generate2WayFromMxN(
    matrix(
      data = c(
        1,2,3,
        4,5,6,
        7,8,9
      ),
      nrow = 3,
      ncol = 3,
      byrow = TRUE
    )
  )
  testthat::expect_type(results, "list")
})

test_that("case control with ensembl database as control produces expected table", {
  results <- PhenoGenRLib::varianceCCAnalysisEnsembl(
    MappedVariants,
    RsIDs,
    totalCaseSamples = 8
  )
  expect_type(results, data.frame)
})
