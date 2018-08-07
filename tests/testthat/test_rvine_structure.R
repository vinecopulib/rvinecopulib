context("Class 'rvine_structure'")

mat <- matrix(c(1, 2, 3, 4, 1, 2, 3, 0, 1, 2, 0, 0, 1, 0, 0, 0), 4, 4)

test_that("constructor and as/is generics work", {
    expect_silent(as.rvine_structure(mat))
    expect_silent(as.rvine_matrix(as.rvine_structure(mat)))
    expect_silent(struct <- rvine_structure(order = 1:4, 
                                            struct_array = list(c(1, 1, 1), 
                                                                c(2, 2), 
                                                                3)))
    expect_true(is.rvine_structure(struct))
    expect_true(is.rvine_matrix(as.rvine_matrix(struct)))
})

test_that("print/dim generics work", {
    struct <- rvine_structure(order = 1:4, 
                              struct_array = list(c(1, 1, 1), 
                                                  c(2, 2), 
                                                  3))
    expect_output(print(as.rvine_matrix(struct)))
    expect_output(print(struct))
    expect_equivalent(dim(as.rvine_matrix(struct)), c(4, 3))
    expect_equivalent(dim(struct), c(4, 3))
})
