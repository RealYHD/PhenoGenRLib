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
  testthat::expect_equal(
    results,
    list(
      matrix(data = c(1,2,4,5), nrow = 2, ncol = 2, byrow = TRUE),
      matrix(data = c(1,2,7,8), nrow = 2, ncol = 2, byrow = TRUE),
      matrix(data = c(4,5,7,8), nrow = 2, ncol = 2, byrow = TRUE),
      matrix(data = c(1,3,4,6), nrow = 2, ncol = 2, byrow = TRUE),
      matrix(data = c(1,3,7,9), nrow = 2, ncol = 2, byrow = TRUE),
      matrix(data = c(4,6,7,9), nrow = 2, ncol = 2, byrow = TRUE),
      matrix(data = c(2,3,5,6), nrow = 2, ncol = 2, byrow = TRUE),
      matrix(data = c(2,3,8,9), nrow = 2, ncol = 2, byrow = TRUE),
      matrix(data = c(5,6,8,9), nrow = 2, ncol = 2, byrow = TRUE)
    )
  )
})

test_that("case control with ensembl database as control produces expected table", {
  results <- PhenoGenRLib::varianceCCAnalysisEnsembl(
    MappedVariants,
    RsIDs,
    totalCaseSamples = 8
  )
  expect_type(results, "data.frame")
})
