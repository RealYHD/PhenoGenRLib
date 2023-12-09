library(PhenoGenRLib)

testthat::test_that("pairwise table generation produces correct permutations of 3x3", {
  results <- PhenoGenRLib::generate2WayFromMxN(
    data.frame(
      C1 = c(1,2,3),
      C2 = c(4,5,6),
      C3 = c(7,8,9),
      row.names = c("R1", "R2", "R3")
    )
  )
  testthat::expect_type(results, "list")
  testthat::expect_setequal(
    object = results,
    expected = list(
      data.frame(C1 = c(1,2), C2 = c(4,5), row.names = c("R1", "R2")),
      data.frame(C1 = c(1,3), C2 = c(4,6), row.names = c("R1", "R3")),
      data.frame(C1 = c(2,3), C2= c(5,6), row.names = c("R2", "R3")),
      data.frame(C1 = c(1,2), C3 = c(7,8), row.names = c("R1", "R2")),
      data.frame(C1 = c(1,3), C3 = c(7,9), row.names = c("R1", "R3")),
      data.frame(C1 = c(2,3), C3 = c(8,9), row.names = c("R2", "R3")),
      data.frame(C2 = c(4,5), C3 = c(7,8), row.names = c("R1", "R2")),
      data.frame(C2 = c(4,6), C3 = c(7,9), row.names = c("R1", "R3")),
      data.frame(C2 = c(5,6), C3 = c(8,9), row.names = c("R2", "R3"))
    )
  )
})

testthat::test_that("Verify the multiple tests for correlation function returns correct p-values", {
  results <- PhenoGenRLib::multipleAssociationTests(
    list(
      data.frame(C1 = c(1,2), C2 = c(4,28), row.names = c("R1", "R2")),
      data.frame(C1 = c(4,2), C2 = c(4,28), row.names = c("R1", "R3"))
    )
  )
  testthat::expect_s3_class(results, "data.frame")
  testthat::expect_equal(
    object = results,
    expected = data.frame(
      test_group = c("Untitled", "Untitled"),
      groups = c("C1, C2", "C1, C2"),
      categories = c("R1, R2", "R1, R3"),
      p_value = c(0.37967, 0.01164),
      method = c(
        "Fisher's Exact Test for Count Data",
        "Fisher's Exact Test for Count Data"
      )
    ),
    tolerance = 1.0e-4
  )
})

test_that("case control with ensembl database as control results not changed", {
  results <- PhenoGenRLib::varianceCCAnalysisEnsembl(
    MappedVariants,
    RsIDs,
    totalCaseSamples = 8
  )
  testthat::expect_s3_class(results, "data.frame")
  testthat::expect_snapshot(results)
})
