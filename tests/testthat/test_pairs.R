test_that("pairs_copula_data works", {
  u <- replicate(3, runif(10))
  pairs_copula_data(u)
})
