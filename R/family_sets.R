family_set_onepar <- c(
    "gaussian", "clayton", "gumbel", "frank", "joe"
)

family_set_bb <- c(
    "bb1", "bb6", "bb7", "bb8"
)

family_set_twopar <- c(
    "student", family_set_bb
)

family_set_parametric <- c(
    "indep", family_set_onepar, family_set_twopar
)

family_set_nonparametric <- c(
    "indep", "tll0"
)

family_set_all <- unique(
    c(family_set_parametric, family_set_nonparametric)
)