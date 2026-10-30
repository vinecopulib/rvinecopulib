## ---------------------------------------------------------------------------
## Replication script for
##
##   Nagler T, Vatter T. "rvinecopulib: A Fast and Comprehensive Implementation
##   of Vine Copula Models in C++ and R." Journal of Statistical Software.
##
## Run with:  Rscript code.R
##
## Set RVCL_FAST=1 for a reduced grid that completes in a few minutes.  The
## full run reproduces every number in the paper; the timings of Section 7.3
## are only meaningful on an idle machine with uniform cores (see README.md).
## An HTML report of the run is produced by knitr::spin("code.R").
## ---------------------------------------------------------------------------

FAST <- nzchar(Sys.getenv("RVCL_FAST"))

options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
set.seed(20260826)

library("rvinecopulib")

## ===========================================================================
## Section 4.1: bivariate copulas
## ===========================================================================

cop <- bicop_dist("bb1", rotation = 90, parameters = c(3, 2))
cop

par_to_ktau(cop)
tail_dep(cop)
blomqvist_beta(cop)

## ===========================================================================
## Section 4.2: vine structures
## ===========================================================================

structure <- dvine_structure(order = c(3, 1, 4, 2))
as_rvine_matrix(structure)

## ===========================================================================
## Section 4.4: the two data layouts for discrete variables
##
## The expanded layout carries 2d columns, the compact layout d + k where k is
## the number of discrete variables.  Both must give the same density.
## ===========================================================================

pair_copulas <- list(
  list(bicop_dist("gaussian", 0, 0.5), bicop_dist("clayton", 0, 1.5)),
  list(bicop_dist("frank", 0, 2))
)
model <- vinecop_dist(
  pair_copulas,
  dvine_structure(1:3),
  var_types = c("c", "d", "c")
)
x <- cbind(rnorm(500), rpois(500, 1), rnorm(500))
u_expanded <- cbind(
  pnorm(x[, 1]),
  ppois(x[, 2], 1),
  pnorm(x[, 3]),
  pnorm(x[, 1]),
  ppois(x[, 2] - 1, 1),
  pnorm(x[, 3])
)
u_compact <- u_expanded[, c(1, 2, 3, 5)]
all.equal(dvinecop(u_expanded, model), dvinecop(u_compact, model))

## ===========================================================================
## Section 4.7: derivatives and uncertainty
## ===========================================================================

## The stepwise gradient is (near) zero at a sequentially fitted model; the
## joint gradient generally is not.
S <- matrix(0.4, 4, 4)
diag(S) <- 1
u <- pseudo_obs(matrix(rnorm(800 * 4), 800, 4) %*% chol(S))
fit <- vinecop(u, family_set = "parametric", keep_data = TRUE)

max(abs(colSums(scores(u, fit, step_wise = TRUE))))
max(abs(colSums(scores(u, fit, step_wise = FALSE))))

## The sandwich is the only form offered; step_wise selects the objective.
round(sqrt(diag(vcov(fit))), 3)
round(sqrt(diag(vcov(fit, step_wise = FALSE))), 3)
round(confint(fit, level = 0.95), 3)

## ===========================================================================
## Section 5.1: a user-defined zero-inflated log-normal margin family
## ===========================================================================

fit_zilnorm <- function(x, weights = numeric(), type = "zi") {
  if (length(weights) == 0) {
    weights <- rep(1, length(x))
  }
  pos <- x > 0
  p0 <- sum(weights[!pos]) / sum(weights)
  ml <- sum(weights[pos] * log(x[pos])) / sum(weights[pos])
  sl <- sqrt(sum(weights[pos] * (log(x[pos]) - ml)^2) / sum(weights[pos]))
  margin_dist(
    d = function(q) {
      ifelse(q > 0, (1 - p0) * dlnorm(q, ml, sl), ifelse(q == 0, p0, 0))
    },
    p = function(q) {
      ifelse(q > 0, p0 + (1 - p0) * plnorm(q, ml, sl), ifelse(q < 0, 0, p0))
    },
    q = function(p) {
      ifelse(p <= p0, 0, qlnorm(pmax((p - p0) / (1 - p0), 0), ml, sl))
    },
    family = "zilnorm",
    type = "zi",
    support = c(0, Inf),
    npars = 3,
    loglik = sum(weights[!pos]) *
      log(p0) +
      sum(weights[pos] * (log(1 - p0) + dlnorm(x[pos], ml, sl, log = TRUE)))
  )
}
zilnorm <- margin_family(fit_zilnorm, family_name = "zilnorm", types = "zi")

## a check on simulated data; the paper fits it to the insurance claims below
z <- ifelse(runif(400) < 0.6, 0, rlnorm(400, 1, 0.5))
str(margin_info(fit_margin(zilnorm, z))[c("npars", "loglik")])

## ===========================================================================
## Section 5.2: a user-supplied tree-selection criterion
##
## Chatterjee's xi, symmetrized because a tree edge is undirected; it rejects
## weights rather than ignoring them.
## ===========================================================================

xi_one <- function(x, y) {
  n <- length(x)
  r <- rank(y[order(x)], ties.method = "max")
  1 - 3 * sum(abs(diff(r))) / (n^2 - 1)
}
xi <- function(data, weights) {
  if (length(weights) > 0) {
    stop("xi does not support weights")
  }
  max(xi_one(data[, 1], data[, 2]), xi_one(data[, 2], data[, 1]))
}
fit_xi <- vinecop(u, tree_crit = xi, family_set = "parametric")
fit_xi

## ===========================================================================
## Section 4.7: coverage of the Wald intervals under known vs estimated margins
##
## The two-stage approximation, not an implementation error, explains the
## shortfall reported in the paper.
## ===========================================================================

coverage_study <- function(R = 300, n = 1000, d = 3, rho = 0.5) {
  truth <- c(rep(rho, d - 1), rho / (1 + rho))
  out <- list(known = matrix(NA, R, d), ranks = matrix(NA, R, d))
  S <- matrix(rho, d, d)
  diag(S) <- 1
  for (r in seq_len(R)) {
    x <- matrix(rnorm(n * d), n, d) %*% chol(S)
    for (nm in c("known", "ranks")) {
      uu <- if (nm == "known") pnorm(x) else pseudo_obs(x)
      f <- vinecop(
        uu,
        family_set = "gaussian",
        structure = dvine_structure(1:d),
        keep_data = TRUE
      )
      est <- unlist(lapply(f$pair_copulas, function(t) {
        lapply(t, function(p) p$parameters)
      }))
      se <- sqrt(diag(vcov(f)))
      out[[nm]][r, ] <- (truth >= est - 1.96 * se) & (truth <= est + 1.96 * se)
    }
  }
  ## per-parameter coverage, the pooled average, and the standard error of the
  ## paired difference -- the three quantities reported in Section 4.7
  diff <- rowMeans(out$known) - rowMeans(out$ranks)
  list(
    per_parameter = rbind(
      known = colMeans(out$known),
      ranks = colMeans(out$ranks)
    ),
    pooled = c(known = mean(out$known), ranks = mean(out$ranks)),
    difference = c(estimate = mean(diff), se = sd(diff) / sqrt(R)),
    mc_se_of_one_proportion = sqrt(0.95 * 0.05 / R)
  )
}

set.seed(11)
coverage_study(R = if (FAST) 30 else 1000)

## ===========================================================================
## Section 3.1: how much of the h-function work the structure index avoids
##
## Reimplements the backend index in R and counts the h-functions some later
## tree actually consumes, against computing both at every edge.
## ===========================================================================

needed_frac <- function(st) {
  d <- st$d
  sa <- st$struct_array
  minarr <- vector("list", length(sa))
  for (t in seq_along(sa)) {
    minarr[[t]] <- if (t == 1) {
      sa[[1]]
    } else {
      pmin(sa[[t]], minarr[[t - 1]][seq_along(sa[[t]])])
    }
  }
  need1 <- lapply(sa, function(r) logical(length(r)))
  need2 <- lapply(sa, function(r) logical(length(r)))
  for (t in 2:(d - 1)) {
    for (e in seq_len(d - t)) {
      need2[[t - 1]][e] <- TRUE
      m <- minarr[[t]][e]
      if (identical(m, sa[[t]][e])) {
        need2[[t - 1]][m] <- TRUE
      } else {
        need1[[t - 1]][m] <- TRUE
      }
    }
  }
  avail <- sum(sapply(sa[1:(d - 2)], length)) * 2
  (sum(unlist(need1)) + sum(unlist(need2))) / avail
}

set.seed(3)
for (d in c(5, 10, 20, 50)) {
  rand <- replicate(if (FAST) 20 else 200, needed_frac(rvine_structure_sim(d)))
  cat(sprintf(
    "d=%-3d random %.3f   C-vine %.3f   D-vine %.3f\n",
    d,
    mean(rand),
    needed_frac(cvine_structure(1:d)),
    needed_frac(dvine_structure(1:d))
  ))
}

## ===========================================================================
## Section 6.1: financial returns
## ===========================================================================

library("qrmdata")
library("xts")
data("SP500_const", package = "qrmdata")

prices <- SP500_const["2010-01-01/2015-12-31"]
prices <- prices[, colSums(is.na(prices)) == 0]
set.seed(20260826)
tickers <- sort(sample(colnames(prices), 50))
returns <- diff(log(prices[, tickers]))[-1, ]
u_ret <- pseudo_obs(as.matrix(returns))

settings <- list(
  list(
    lab = "parametric, trunc_lvl = NA",
    fs = "parametric",
    tc = "tau",
    tl = NA
  ),
  list(lab = "all families, trunc_lvl = NA", fs = "all", tc = "tau", tl = NA),
  list(
    lab = "parametric, mcor, trunc_lvl = NA",
    fs = "parametric",
    tc = "mcor",
    tl = NA
  ),
  list(
    lab = "parametric, trunc_lvl = 10",
    fs = "parametric",
    tc = "tau",
    tl = 10
  )
)
if (!FAST) {
  ret_times <- numeric(length(settings))
  ret_fits <- lapply(seq_along(settings), function(i) {
    r <- settings[[i]]
    t0 <- Sys.time()
    f <- vinecop(
      u_ret,
      family_set = r$fs,
      tree_crit = r$tc,
      trunc_lvl = r$tl,
      selcrit = "mbicv",
      cores = 4
    )
    ret_times[i] <<- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    f
  })
  ## seconds per fit on four cores; the first is quoted in Section 6.1
  print(round(ret_times))

  ## Table: fitted vine copulas
  print(do.call(
    rbind,
    Map(
      function(r, f) {
        data.frame(
          setting = r$lab,
          trees = f$structure$trunc_lvl,
          non_indep = sum(unlist(get_all_families(f)) != "indep"),
          npars = round(f$npars, 1),
          loglik = round(f$loglik),
          bic = round(BIC(f)),
          mbicv = round(mBICV(f))
        )
      },
      settings,
      ret_fits
    )
  ))

  ## the nonparametric family is never selected at this dimension
  print("tll" %in% unlist(get_all_families(ret_fits[[2]])))

  ## the Section 4.3 call is the second setting
  print(head(summary(ret_fits[[2]]), 4))

  ## families of the selected model, by tree
  smr <- summary(ret_fits[[1]])
  print(sort(table(smr$family), decreasing = TRUE))
  print(table(smr$tree, smr$family))

  ## truncation path of the full fit (Figure: truncation.pdf)
  full <- ret_fits[[1]]
  path <- sapply(1:49, function(k) {
    mBICV(truncate_model(full, k), newdata = u_ret)
  })
  print(round(path))
  print(c(selected = full$structure$trunc_lvl, path_min = which.min(path)))

  ## Whole-model comparison, out of sample.  Information criteria compare
  ## parametric candidates; they are not a basis for choosing between a
  ## parametric family and a nonparametric estimator, so we split the sample.
  train <- u_ret[1:1000, ]
  test <- u_ret[1001:1509, ]
  heldout <- function(m) mean(dvinecop(test, m, log = TRUE))
  timed <- function(expr) {
    secs <- system.time(f <- expr)[["elapsed"]]
    attr(f, "secs") <- secs
    f
  }
  fits_ho <- list(
    parametric_itau = timed(vinecop(
      train,
      family_set = "itau",
      par_method = "itau",
      selcrit = "mbicv",
      trunc_lvl = NA,
      cores = 8
    )),
    parametric_mle = timed(vinecop(
      train,
      family_set = "parametric",
      selcrit = "mbicv",
      trunc_lvl = NA,
      cores = 8
    ))
  )
  ## nonparametric vines at the parametric level and at two shallower ones
  lvls <- unique(c(fits_ho$parametric_itau$structure$trunc_lvl, 3, 1))
  for (k in lvls) {
    fits_ho[[paste0("nonparametric_", k)]] <- timed(vinecop(
      train,
      family_set = "tll",
      trunc_lvl = k,
      cores = 8
    ))
  }
  print(data.frame(
    model = names(fits_ho),
    trees = sapply(fits_ho, function(m) m$structure$trunc_lvl),
    npars = round(sapply(fits_ho, function(m) m$npars), 1),
    fit_secs = round(sapply(fits_ho, attr, "secs"), 1),
    heldout = round(sapply(fits_ho, heldout), 2)
  ))
}

## ===========================================================================
## Section 6.2: insurance claims
## ===========================================================================

library("insuranceData")
data("dataCar", package = "insuranceData")

set.seed(20260826)
idx <- sample(nrow(dataCar), if (FAST) 3000 else 10000)
dat <- data.frame(
  cost = zero_inflated(dataCar$claimcst0[idx]),
  veh_value = dataCar$veh_value[idx],
  exposure = dataCar$exposure[idx],
  veh_age = ordered(dataCar$veh_age[idx]),
  agecat = ordered(dataCar$agecat[idx])
)
predictors <- c("veh_value", "exposure", "veh_age", "agecat")

## claim cost is positive exactly when the claim count is, so only one of the
## two can enter the model
all((dataCar$claimcst0 > 0) == (dataCar$numclaims > 0))

## the user-defined family of Section 5.1, fitted directly
str(margin_info(fit_margin(zilnorm, dat$cost))[c("npars", "loglik")])

fit_ins <- vine(
  dat,
  margins_controls = list(
    family_set = list(
      cost = zilnorm,
      veh_value = c("kde1d", "norm", "sstd", "std", "logis", "cauchy"),
      exposure = c("kde1d", "beta", "norm", "unif"),
      veh_age = "kde1d",
      agecat = "kde1d"
    ),
    selcrit = "bic"
  ),
  copula_controls = list(
    family_set = "parametric",
    selcrit = "bic",
    conditioning_set = predictors,
    keep_data = TRUE
  ),
  keep_data = TRUE,
  cores = 4
)

## Table: selected margins
do.call(
  rbind,
  lapply(seq_along(dat), function(i) {
    m <- margin_info(fit_ins$margins[[i]])
    data.frame(
      variable = names(dat)[i],
      type = m$type,
      family = m$family_name,
      npars = round(m$npars, 2),
      loglik = round(m$loglik, 1)
    )
  })
)

summary(fit_ins$copula)

## the modified BIC would not truncate this model
round(
  sapply(1:4, function(k) {
    mBICV(truncate_model(fit_ins$copula, k), newdata = fit_ins$copula$data)
  }),
  1
)

## Wald intervals condition on the fitted margins; see ?parameter_uncertainty
round(confint(fit_ins$copula), 3)

## Rosenblatt residuals: uniform margins and no remaining dependence
set.seed(1)
res <- rosenblatt(dat, fit_ins)
pv <- function(x) suppressWarnings(ks.test(x, "punif"))$p.value
round(apply(res, 2, pv), 2)
round(max(abs(cor(res, method = "kendall")[upper.tri(diag(5))])), 3)

## Conditional simulation on the original data scale.
profile <- data.frame(
  veh_value = 1.5,
  exposure = 0.75,
  veh_age = ordered(2, levels = levels(dat$veh_age)),
  agecat = ordered(2, levels = levels(dat$agecat))
)
n_sim <- if (FAST) 1e5 else 1e6
set.seed(20260826)
sim_cond <- rvine(
  n_sim,
  fit_ins,
  x_cond = profile,
  conditioning_set = predictors
)
sim_marg <- rvine(n_sim, fit_ins)

c(
  conditional = mean(sim_cond[, "cost"] > 0),
  unconditional = mean(sim_marg[, "cost"] > 0)
)
c(
  conditional = mean(sim_cond[, "cost"]),
  unconditional = mean(sim_marg[, "cost"])
)
c(
  conditional = median(sim_cond[sim_cond[, "cost"] > 0, "cost"]),
  unconditional = median(sim_marg[sim_marg[, "cost"] > 0, "cost"])
)

## the conditioning values are reproduced exactly
all.equal(unname(sim_cond[, "veh_value"]), rep(1.5, n_sim))

## ===========================================================================
## Section 7: validation and benchmarks
## ===========================================================================

library("VineCopula")
library("kdecopula")
library("bench")

## --- Bivariate parity against VineCopula --------------------------------
## The full driver is in parity.R, which writes the table reproduced in the
## manuscript.  Run it with
##
##   Rscript parity.R
##
## Its headline result is that everything except Kendall's tau agrees to
## between 1e-16 and 1e-11.  Frank's tau differs by 7.8e-04; a high-precision
## reference settles which implementation is responsible:
theta <- 3
debye1 <- integrate(
  function(t) ifelse(t < 1e-12, 1, t / expm1(t)),
  0,
  theta,
  rel.tol = 1e-12,
  subdivisions = 2000
)$value /
  theta
reference <- 1 - 4 / theta + 4 * debye1 / theta
c(
  reference = reference,
  rvinecopulib = par_to_ktau(bicop_dist("frank", 0, theta)),
  VineCopula = BiCopPar2Tau(5, theta)
)

## Joe's tau has a removable singularity at theta = 2; the limit is 2 - pi^2/6.
c(value = par_to_ktau(bicop_dist("joe", 0, 2)), limit = 2 - pi^2 / 6)

## --- Runtime, parametric track ------------------------------------------

## Every ratio in Section 7 is quoted only after checking that the two
## implementations reached a comparable optimum on the same data: a speed
## difference between models that differ is not a speed difference.  The fitted
## models are captured out of the timed expressions themselves via the assign()
## calls below, so this reports on exactly the fits that were measured rather
## than on a separate pair fitted for the purpose.
cap <- new.env()

## RVineStructureSelect() leaks its socket cluster on every call with
## cores > 1: each call leaves `cores` worker processes running.  Left alone
## they accumulate across repetitions and across cells, so a later measurement
## competes with the debris of every earlier one -- and only the multi-core
## cells are affected, which is exactly where the comparison lives.  Reaping
## between repetitions keeps every timing starting from the same machine state.
## The reaping itself is outside the measured region.
reap_cluster <- function() {
  cl <- tryCatch(parallel::getDefaultCluster(), error = function(e) NULL)
  if (!is.null(cl)) {
    try(parallel::stopCluster(cl), silent = TRUE)
    try(parallel::setDefaultCluster(NULL), silent = TRUE)
  }
  invisible(NULL)
}

## We report the fastest of the repetitions rather than the median.  Anything
## else running on the machine can only ever make a repetition slower, so the
## minimum is the best available estimate of the uncontended cost, and it also
## discards the first-call warm-up.  Medians of three repetitions moved by up
## to 37% between two runs of this script; minima are stable to a few percent.
##
## `f` is a function, not an expression: a promise would be forced once and
## every later repetition would measure nothing.
time_best <- function(f, reps) {
  t <- numeric(reps)
  for (i in seq_len(reps)) {
    t[i] <- system.time(f())[["elapsed"]]
    reap_cluster()
  }
  min(t)
}

agreement <- function(u) {
  f <- get("f_rv", envir = cap)
  v <- get("f_vc", envir = cap)
  ll_rv <- f$loglik
  ll_vc <- RVineLogLik(u, v, verbose = FALSE, calculate.V = FALSE)$loglik
  c(
    d = ncol(u),
    rv_loglik = ll_rv,
    vc_loglik = ll_vc,
    gap = ll_rv - ll_vc,
    rel_gap = abs(ll_rv - ll_vc) / abs(ll_rv),
    rv_npc = sum(summary(f)$family != "indep"),
    vc_npc = sum(v$family != 0)
  )
}

RV_FAM <- c(
  "indep",
  "gaussian",
  "t",
  "clayton",
  "gumbel",
  "frank",
  "joe",
  "bb1",
  "bb6",
  "bb7",
  "bb8"
)
VC_FAM <- c(0:10, 13, 14, 16:20, 23, 24, 26:30, 33, 34, 36:40)

runtime_cell <- function(d, n, iterations) {
  u <- u_ret[seq_len(n), seq_len(d), drop = FALSE]
  rv <- function(cores) {
    function() {
      f <- vinecop(
        u,
        family_set = RV_FAM,
        selcrit = "aic",
        tree_crit = "tau",
        cores = cores
      )
      if (cores == 1) {
        assign("f_rv", f, envir = cap)
      }
      invisible(f)
    }
  }
  vc <- function(cores) {
    function() {
      f <- RVineStructureSelect(
        u,
        familyset = VC_FAM,
        selectioncrit = "AIC",
        treecrit = "tau",
        indeptest = FALSE,
        progress = FALSE,
        cores = cores
      )
      if (cores == 1) {
        assign("f_vc", f, envir = cap)
      }
      invisible(f)
    }
  }
  list(
    time = c(
      rv1 = time_best(rv(1), iterations),
      rv4 = time_best(rv(4), iterations),
      vc1 = time_best(vc(1), iterations),
      vc4 = time_best(vc(4), iterations)
    ),
    fit = agreement(u)
  )
}
if (!FAST) {
  mle <- list(
    runtime_cell(5, 1000, 7),
    runtime_cell(10, 1000, 5),
    runtime_cell(20, 1000, 5)
  )
  ## fastest seconds
  print(do.call(rbind, lapply(mle, `[[`, "time")))
  ## the same fits, checked against each other
  print(do.call(rbind, lapply(mle, `[[`, "fit")))
}

## --- Runtime: tau-inversion fitting, where threading matters most --------
RV_ITAU <- c("indep", "gaussian", "t", "clayton", "gumbel", "frank", "joe")
VC_ITAU <- c(0, 1, 2, 3, 4, 5, 6, 13, 14, 16, 23, 24, 26, 33, 34, 36)

itau_cell <- function(d, n = 1509, iterations = 5) {
  u <- u_ret[seq_len(n), seq_len(d), drop = FALSE]
  rv <- function(cores) {
    function() {
      f <- vinecop(
        u,
        family_set = RV_ITAU,
        par_method = "itau",
        selcrit = "aic",
        cores = cores
      )
      if (cores == 1) {
        assign("f_rv", f, envir = cap)
      }
      invisible(f)
    }
  }
  vc <- function(cores) {
    function() {
      f <- RVineStructureSelect(
        u,
        familyset = VC_ITAU,
        method = "itau",
        selectioncrit = "AIC",
        treecrit = "tau",
        indeptest = FALSE,
        progress = FALSE,
        cores = cores
      )
      if (cores == 1) {
        assign("f_vc", f, envir = cap)
      }
      invisible(f)
    }
  }
  list(
    time = c(
      rv1 = time_best(rv(1), iterations),
      rv4 = time_best(rv(4), iterations),
      rv8 = time_best(rv(8), iterations),
      vc1 = time_best(vc(1), iterations),
      vc4 = time_best(vc(4), iterations),
      vc8 = time_best(vc(8), iterations)
    ),
    fit = agreement(u)
  )
}
if (!FAST) {
  itau <- lapply(c(5, 10, 20, 50), itau_cell)
  print(do.call(rbind, lapply(itau, `[[`, "time")))
  print(do.call(rbind, lapply(itau, `[[`, "fit")))
}

## --- Runtime: evaluation on a fitted model -------------------------------
eval_cell <- function(d, n = 1509) {
  u <- u_ret[seq_len(n), seq_len(d), drop = FALSE]
  f <- vinecop(
    u,
    family_set = RV_ITAU,
    par_method = "itau",
    selcrit = "aic",
    cores = 8
  )
  v <- RVineStructureSelect(
    u,
    familyset = VC_ITAU,
    method = "itau",
    selectioncrit = "AIC",
    treecrit = "tau",
    indeptest = FALSE,
    progress = FALSE,
    cores = 8
  )
  ## the two fits above are setup, not measurement; drop the leaked workers
  ## before timing anything
  reap_cluster()
  b <- bench::mark(
    rv_pdf1 = dvinecop(u, f, cores = 1),
    rv_pdf8 = dvinecop(u, f, cores = 8),
    vc_pdf = RVinePDF(u, v),
    rv_ll1 = sum(dvinecop(u, f, cores = 1, log = TRUE)),
    rv_ll8 = sum(dvinecop(u, f, cores = 8, log = TRUE)),
    vc_ll = RVineLogLik(u, v)$loglik,
    rv_ros1 = rosenblatt(u, f, cores = 1),
    rv_ros8 = rosenblatt(u, f, cores = 8),
    vc_ros = RVinePIT(u, v),
    rv_inv1 = inverse_rosenblatt(u, f, cores = 1),
    rv_inv8 = inverse_rosenblatt(u, f, cores = 8),
    vc_inv = RVineSim(n, v, U = u),
    rv_sim1 = rvinecop(n, f, cores = 1),
    rv_sim8 = rvinecop(n, f, cores = 8),
    vc_sim = RVineSim(n, v),
    iterations = 10,
    check = FALSE,
    filter_gc = FALSE
  )
  setNames(as.numeric(b$min), as.character(b$expression))
}
if (!FAST) {
  print(t(sapply(c(5, 10, 20, 50), eval_cell)))
}

## --- Nonparametric track -------------------------------------------------
u2 <- u_ret[1:1000, 1:2]
tll_fit <- bicop(u2, family_set = "tll", nonpar_method = "quadratic")
kde_fit <- kdecop(u2, method = "TLL2")
gr <- as.matrix(expand.grid(
  seq(0.05, 0.95, length.out = 100),
  seq(0.05, 0.95, length.out = 100)
))
bench::mark(
  rv_pdf = dbicop(gr, tll_fit),
  kd_pdf = dkdecop(gr, kde_fit),
  rv_h = hbicop(gr, 1, tll_fit),
  kd_h = hkdecop(gr, kde_fit, cond.var = 2),
  rv_hinv = hbicop(gr, 1, tll_fit, inverse = TRUE),
  kd_hinv = hkdecop(gr, kde_fit, cond.var = 2, inverse = TRUE),
  iterations = 5,
  check = FALSE,
  filter_gc = FALSE
)[, c("expression", "min")]

## Accuracy against a known truth, so the timings above cannot be explained by
## one implementation doing less work than the other.
tll_accuracy <- function(tau, n, reps) {
  truth <- bicop_dist("clayton", 0, ktau_to_par("clayton", tau))
  g <- as.matrix(expand.grid(
    seq(0.01, 0.99, length.out = 60),
    seq(0.01, 0.99, length.out = 60)
  ))
  d0 <- dbicop(g, truth)
  err <- replicate(reps, {
    x <- rbicop(n, truth)
    c(
      rv = mean(abs(
        dbicop(g, bicop(x, family_set = "tll", nonpar_method = "quadratic")) -
          d0
      )),
      kd = mean(abs(dkdecop(g, kdecop(x, method = "TLL2")) - d0))
    )
  })
  ratio <- err["rv", ] / err["kd", ]
  c(
    rv = mean(err["rv", ]),
    kd = mean(err["kd", ]),
    excess_pct = 100 * (mean(ratio) - 1),
    se_pct = 100 * sd(ratio) / sqrt(reps)
  )
}

set.seed(4242)
tll_grid <- expand.grid(n = c(500, 1000, 2000), tau = c(0.3, 0.6))
cbind(
  tll_grid,
  t(mapply(
    function(tau, n) tll_accuracy(tau, n, if (FAST) 10 else 200),
    tll_grid$tau,
    tll_grid$n
  ))
)

## ===========================================================================
sessionInfo()
