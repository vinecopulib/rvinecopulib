family_set_archimedean <- c(
    "clayton", "gumbel", "frank", "joe"
)

family_set_elliptical <- c(
    "gaussian", "t"
)

family_set_bb <- c(
    "bb1", "bb6", "bb7", "bb8"
)

family_set_onepar <- c(
    "gaussian", family_set_archimedean
)

family_set_twopar <- c(
    "t", family_set_bb
)

family_set_parametric <- c(
    "indep", family_set_onepar, family_set_twopar
)

family_set_nonparametric <- c(
    "indep", "tll"
)

family_set_itau <- c(
    "t", family_set_onepar
)

family_set_all <- unique(
    c(family_set_parametric, family_set_nonparametric)
)

family_set_defs <- c(
    "archimedean", "elliptical", "bbs", "oneparametric", "twoparametric", 
    "parametric", "nonparametric", "itau", "all"
)

family_set_all_defs <- c(
    family_set_all, family_set_defs
)