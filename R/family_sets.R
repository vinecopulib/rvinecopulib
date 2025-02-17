family_set_archimedean <- c(
  "clayton", "gumbel", "frank", "joe"
)

family_set_elliptical <- c(
  "gaussian", "t"
)

family_set_extreme_value <- c(
  "tawn", "gumbel"
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

family_set_threepar <- c(
  "tawn"
)

family_set_parametric <- c(
  "indep", family_set_onepar, family_set_twopar, family_set_threepar
)

family_set_nonparametric <- c(
  "indep", "tll"
)

family_set_itau <- c(
  "indep", family_set_onepar, "t"
)

family_set_rotationless <- c(
  "frank", family_set_elliptical, family_set_nonparametric
)

family_set_all <- unique(
  c(family_set_parametric, family_set_nonparametric)
)

family_set_defs <- c(
  "archimedean", "elliptical", "ev", "bbs", "oneparametric", "twoparametric",
  "threeparametric", "parametric", "nonparametric", "itau", "all"
)

family_set_all_defs <- c(
  family_set_all, family_set_defs
)
