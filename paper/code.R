## ---------------------------------------------------------------------------
## Replication script for
##
##   Nagler T, Vatter T. "rvinecopulib: A Fast and Comprehensive Implementation
##   of Vine Copula Models in C++ and R." Journal of Statistical Software.
##
## Run with:  Rscript code.R
##
## Set FAST <- TRUE (or export RVCL_FAST=1) for the reduced grid, which
## completes in a few minutes.  The full run reproduces every number in the
## paper and is intended to stay under one hour on a recent desktop.
## ---------------------------------------------------------------------------

FAST <- nzchar(Sys.getenv("RVCL_FAST"))

options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
set.seed(20260826)

library("rvinecopulib")

## ===========================================================================
## Section 4.1: bivariate copulas
## ===========================================================================

## A fixed pair copula, printed as it appears in the paper.
cop <- bicop_dist("bb1", rotation = 90, parameters = c(3, 2))
cop

## Dependence summaries.
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

n <- 500
x <- cbind(rnorm(n), rpois(n, 1), rnorm(n))
model <- vinecop_dist(
  list(
    list(bicop_dist("gaussian", 0, 0.5), bicop_dist("clayton", 0, 1.5)),
    list(bicop_dist("frank", 0, 2))
  ),
  dvine_structure(1:3),
  var_types = c("c", "d", "c")
)

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
## Section 4.7 and 5.1: derivatives, uncertainty, and the margin protocol
## ===========================================================================

## The stepwise gradient is (near) zero at a sequentially fitted model; the
## joint gradient generally is not.
S <- matrix(0.4, 4, 4)
diag(S) <- 1
u <- pseudo_obs(matrix(rnorm(800 * 4), 800, 4) %*% chol(S))
fit <- vinecop(u, family_set = "parametric", keep_data = TRUE)

max(abs(colSums(scores(u, fit, step_wise = TRUE))))
max(abs(colSums(scores(u, fit, step_wise = FALSE))))

## Covariance and Wald intervals.
## the sandwich is the only form offered; step_wise selects the objective
sqrt(diag(vcov(fit)))
sqrt(diag(vcov(fit, step_wise = FALSE)))
confint(fit, level = 0.95)

## A user-defined zero-inflated log-normal margin family.
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

z <- ifelse(runif(400) < 0.6, 0, rlnorm(400, 1, 0.5))
margin_info(fit_margin(zilnorm, z))

## ===========================================================================
## Section 5.2: a user-supplied tree-selection criterion
##
## Chatterjee's xi, which detects functional relationships that rank
## correlations miss.
## ===========================================================================

xi <- function(data, weights) {
  n <- nrow(data)
  r <- rank(data[order(data[, 1]), 2], ties.method = "max")
  1 - 3 * sum(abs(diff(r))) / (n^2 - 1)
}
vinecop(u, tree_crit = xi, family_set = "parametric")

## ===========================================================================
## Section 4.7: coverage of the Wald intervals under known vs estimated margins
##
## Establishes that the two-stage approximation, not an implementation error,
## explains the shortfall reported in the paper.
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
      se <- sqrt(diag(vcov(f, type = "sandwich")))
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
## Section 6: illustrations
## ===========================================================================

## --- Case study A: S&P 500 returns --------------------------------------
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
    lab = "parametric, untruncated",
    fs = "parametric",
    tc = "tau",
    tl = Inf
  ),
  list(lab = "parametric, mBICV", fs = "parametric", tc = "tau", tl = NA),
  list(lab = "all families, mBICV", fs = "all", tc = "tau", tl = NA),
  list(lab = "parametric, mBICV, mcor", fs = "parametric", tc = "mcor", tl = NA)
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
  ## the untruncated fit time quoted in Section 6.1
  round(ret_times[1])
  do.call(
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
  )

  ## the nonparametric family is never selected at this dimension
  "tll" %in% unlist(get_all_families(ret_fits[[3]]))

  ## truncation path (Figure: truncation.pdf)
  full <- ret_fits[[1]]
  sapply(1:15, function(k) mBICV(truncate_model(full, k), newdata = u_ret))

  ## Whole-model comparison, out of sample.  Information criteria compare
  ## parametric candidates; they are not a basis for choosing between a
  ## parametric family and a nonparametric estimator, so we split the sample.
  train <- u_ret[1:1000, ]
  test <- u_ret[1001:1509, ]
  heldout <- function(m) mean(log(dvinecop(test, m)))
  fits_ho <- list(
    parametric_itau = vinecop(
      train,
      family_set = "itau",
      par_method = "itau",
      selcrit = "mbicv",
      trunc_lvl = NA,
      cores = 8
    ),
    parametric_mle = vinecop(
      train,
      family_set = "parametric",
      selcrit = "mbicv",
      trunc_lvl = NA,
      cores = 8
    ),
    nonparametric_9 = vinecop(
      train,
      family_set = "tll",
      trunc_lvl = 9,
      cores = 8
    ),
    nonparametric_3 = vinecop(
      train,
      family_set = "tll",
      trunc_lvl = 3,
      cores = 8
    ),
    nonparametric_1 = vinecop(
      train,
      family_set = "tll",
      trunc_lvl = 1,
      cores = 8
    )
  )
  data.frame(
    model = names(fits_ho),
    npars = round(sapply(fits_ho, function(m) m$npars), 1),
    heldout = round(sapply(fits_ho, heldout), 2)
  )
}


## --- Case study B: insurance claims -------------------------------------
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

## intervals condition on the fitted margins; see ?parameter_uncertainty
confint(fit_ins$copula)

## Conditional simulation on the original data scale.
profile <- data.frame(
  veh_value = 1.5,
  exposure = 0.75,
  veh_age = ordered(2, levels = levels(dat$veh_age)),
  agecat = ordered(2, levels = levels(dat$agecat))
)
sim_cond <- rvine(
  20000,
  fit_ins,
  x_cond = profile,
  conditioning_set = predictors
)
sim_marg <- rvine(20000, fit_ins)

c(
  conditional = mean(sim_cond[, "cost"] > 0),
  unconditional = mean(sim_marg[, "cost"] > 0)
)
c(
  conditional = mean(sim_cond[, "cost"]),
  unconditional = mean(sim_marg[, "cost"])
)

## the conditioning values are reproduced exactly
all.equal(unname(sim_cond[, "veh_value"]), rep(1.5, 20000))

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
  b <- bench::mark(
    rv1 = vinecop(
      u,
      family_set = RV_FAM,
      selcrit = "aic",
      tree_crit = "tau",
      cores = 1
    ),
    rv4 = vinecop(
      u,
      family_set = RV_FAM,
      selcrit = "aic",
      tree_crit = "tau",
      cores = 4
    ),
    vc1 = RVineStructureSelect(
      u,
      familyset = VC_FAM,
      selectioncrit = "AIC",
      treecrit = "tau",
      indeptest = FALSE,
      progress = FALSE,
      cores = 1
    ),
    vc4 = RVineStructureSelect(
      u,
      familyset = VC_FAM,
      selectioncrit = "AIC",
      treecrit = "tau",
      indeptest = FALSE,
      progress = FALSE,
      cores = 4
    ),
    iterations = iterations,
    check = FALSE,
    filter_gc = FALSE
  )
  setNames(as.numeric(b$median), as.character(b$expression))
}
if (!FAST) {
  rbind(
    runtime_cell(5, 1000, 5),
    runtime_cell(10, 1000, 5),
    runtime_cell(20, 1000, 3)
  )
}

## --- Runtime: tau-inversion fitting, where threading matters most --------
RV_ITAU <- c("indep", "gaussian", "t", "clayton", "gumbel", "frank", "joe")
VC_ITAU <- c(0, 1, 2, 3, 4, 5, 6, 13, 14, 16, 23, 24, 26, 33, 34, 36)

itau_cell <- function(d, n = 1509, iterations = 3) {
  u <- u_ret[seq_len(n), seq_len(d), drop = FALSE]
  b <- bench::mark(
    rv1 = vinecop(
      u,
      family_set = RV_ITAU,
      par_method = "itau",
      selcrit = "aic",
      cores = 1
    ),
    rv4 = vinecop(
      u,
      family_set = RV_ITAU,
      par_method = "itau",
      selcrit = "aic",
      cores = 4
    ),
    rv8 = vinecop(
      u,
      family_set = RV_ITAU,
      par_method = "itau",
      selcrit = "aic",
      cores = 8
    ),
    vc1 = RVineStructureSelect(
      u,
      familyset = VC_ITAU,
      method = "itau",
      selectioncrit = "AIC",
      treecrit = "tau",
      indeptest = FALSE,
      progress = FALSE,
      cores = 1
    ),
    vc4 = RVineStructureSelect(
      u,
      familyset = VC_ITAU,
      method = "itau",
      selectioncrit = "AIC",
      treecrit = "tau",
      indeptest = FALSE,
      progress = FALSE,
      cores = 4
    ),
    vc8 = RVineStructureSelect(
      u,
      familyset = VC_ITAU,
      method = "itau",
      selectioncrit = "AIC",
      treecrit = "tau",
      indeptest = FALSE,
      progress = FALSE,
      cores = 8
    ),
    iterations = iterations,
    check = FALSE,
    filter_gc = FALSE
  )
  setNames(as.numeric(b$median), as.character(b$expression))
}
if (!FAST) {
  t(sapply(c(5, 10, 20, 50), itau_cell))
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
  b <- bench::mark(
    rv_pdf1 = dvinecop(u, f, cores = 1),
    rv_pdf8 = dvinecop(u, f, cores = 8),
    vc_pdf = RVinePDF(u, v),
    rv_ll1 = sum(log(dvinecop(u, f, cores = 1))),
    rv_ll8 = sum(log(dvinecop(u, f, cores = 8))),
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
  setNames(as.numeric(b$median), as.character(b$expression))
}
if (!FAST) {
  t(sapply(c(5, 10, 20, 50), eval_cell))
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
)[, c("expression", "median")]

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
