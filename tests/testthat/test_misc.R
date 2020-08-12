test_that("pairs_copula_data works", {
  u <- replicate(3, runif(10))
  expect_silent(pairs_copula_data(u))
})

test_that("summary_df is printed correctly", {
  x <- data.frame(a = 1:20, B = "A")
  output <- capture_output_lines(x2 <- print.summary_df(x))
  expect_identical(x, x2)
  expect_length(output, 13)
})

