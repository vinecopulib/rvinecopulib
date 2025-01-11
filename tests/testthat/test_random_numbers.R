context("Random numbers")

set.seed(1)

u <- rbicop(2, "gumbel", 90, 3)
u_expected <- matrix(c(0.779635026119649, 0.957448964705691,
                       0.191619620469997, 0.00465105664209572), 2, 2)
expect_equal(u, u_expected)
