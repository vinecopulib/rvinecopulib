## ---------------------------------------------------------------------------
## Figures for the manuscript.  Refits the two case-study models of code.R.
##
##   Rscript make_figures.R
##
## Writes figures/margins.pdf, figures/contours.pdf, figures/condsim.pdf and
## figures/truncation.pdf, each at the text width of the jss class.
## ---------------------------------------------------------------------------

suppressMessages({
  library("rvinecopulib")
  library("insuranceData")
  library("qrmdata")
  library("xts")
})
dir.create("figures", showWarnings = FALSE)
tw <- 6.1 # text width in inches

## --- Case study B: margins, contours and conditional simulation ----------
data("dataCar", package = "insuranceData")
set.seed(20260826)
idx <- sample(nrow(dataCar), 10000)
dat <- data.frame(
  cost = zero_inflated(dataCar$claimcst0[idx]),
  veh_value = dataCar$veh_value[idx],
  exposure = dataCar$exposure[idx],
  veh_age = ordered(dataCar$veh_age[idx]),
  agecat = ordered(dataCar$agecat[idx])
)

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
predictors <- c("veh_value", "exposure", "veh_age", "agecat")

fit <- vine(
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

pdf("figures/margins.pdf", width = tw, height = 2.3, pointsize = 9)
par(
  mfrow = c(1, 3),
  mar = c(3.4, 3.4, 2.0, 0.6),
  mgp = c(2.1, 0.7, 0),
  bty = "n"
)
pos <- dat$cost[dat$cost > 0]
p0 <- 1 - mean(dat$cost > 0)
hist(
  log(pos),
  breaks = 40,
  freq = FALSE,
  col = "grey88",
  border = "white",
  main = "cost (zero-inflated)",
  xlab = "log claim amount, positive part",
  ylab = "density"
)
xs <- seq(min(log(pos)), max(log(pos)), length.out = 300)
lines(xs, dmargin(exp(xs), fit$margins[[1]]) * exp(xs) / (1 - p0), lwd = 2)
legend("topright", bty = "n", legend = sprintf("atom at 0: %.3f", p0))

hist(
  dat$veh_value,
  breaks = 60,
  freq = FALSE,
  col = "grey88",
  border = "white",
  xlim = c(0, 6),
  main = "veh_value (continuous)",
  xlab = "vehicle value",
  ylab = "density"
)
xs <- seq(0.01, 6, length.out = 400)
lines(xs, dmargin(xs, fit$margins[[2]]), lwd = 2)

## ordered factors are fitted on the internal 0-based coding
lv <- seq_along(levels(dat$agecat)) - 1L
pk <- pmargin(lv, fit$margins[[5]]) - pmargin(lv - 1, fit$margins[[5]])
bp <- barplot(
  prop.table(table(as.integer(dat$agecat))),
  col = "grey88",
  border = "white",
  main = "agecat (discrete)",
  xlab = "driver age category",
  ylab = "probability"
)
points(bp, pk, pch = 19, cex = 1.1)
segments(bp, 0, bp, pk, lwd = 1.5)
dev.off()

pdf("figures/contours.pdf", width = tw, height = 1.9, pointsize = 9)
contour(fit$copula, tree = 1)
dev.off()

profile <- data.frame(
  veh_value = 1.5,
  exposure = 0.75,
  veh_age = ordered(2, levels = levels(dat$veh_age)),
  agecat = ordered(2, levels = levels(dat$agecat))
)
set.seed(20260826)
sc <- rvine(1e6, fit, x_cond = profile, conditioning_set = predictors)
su <- rvine(1e6, fit)

pdf("figures/condsim.pdf", width = tw, height = 2.4, pointsize = 9)
par(
  mfrow = c(1, 2),
  mar = c(3.6, 3.6, 2.0, 0.6),
  mgp = c(2.3, 0.7, 0),
  bty = "n"
)
p <- c(
  unconditional = mean(su[, "cost"] > 0),
  conditional = mean(sc[, "cost"] > 0)
)
b <- barplot(
  p,
  ylim = c(0, 0.13),
  col = c("grey80", "grey35"),
  border = "white",
  ylab = "P(claim)",
  main = "probability of a claim"
)
text(b, p + 0.008, sprintf("%.3f", p))
du <- density(log(su[su[, "cost"] > 0, "cost"]))
dc <- density(log(sc[sc[, "cost"] > 0, "cost"]))
plot(
  du,
  lwd = 2,
  col = "grey55",
  main = "claim amount given a claim",
  xlab = "log claim amount",
  ylab = "density",
  ylim = range(0, du$y, dc$y)
)
lines(dc, lwd = 2)
legend(
  "topright",
  bty = "n",
  lwd = 2,
  col = c("grey55", "black"),
  legend = c("unconditional", "given profile")
)
dev.off()

## --- Case study A: truncation path ---------------------------------------
data("SP500_const", package = "qrmdata")
prices <- SP500_const["2010-01-01/2015-12-31"]
prices <- prices[, colSums(is.na(prices)) == 0]
set.seed(20260826)
tickers <- sort(sample(colnames(prices), 50))
u_ret <- pseudo_obs(as.matrix(diff(log(prices[, tickers]))[-1, ]))

full <- vinecop(
  u_ret,
  family_set = "parametric",
  tree_crit = "tau",
  trunc_lvl = Inf,
  selcrit = "mbicv",
  cores = 4
)
auto <- vinecop(
  u_ret,
  family_set = "parametric",
  tree_crit = "tau",
  trunc_lvl = NA,
  selcrit = "mbicv",
  cores = 4
)
sel <- auto$structure$trunc_lvl
ks <- seq_len(ncol(u_ret) - 1)
ll <- sapply(ks, function(k) truncate_model(full, k)$loglik)
mb <- sapply(ks, function(k) mBICV(truncate_model(full, k), newdata = u_ret))
kmin <- which.min(mb)

pdf("figures/truncation.pdf", width = tw, height = 2.4, pointsize = 9)
par(
  mfrow = c(1, 2),
  mar = c(3.6, 3.9, 2.0, 0.6),
  mgp = c(2.4, 0.7, 0),
  bty = "n"
)
plot(
  ks,
  ll,
  type = "l",
  lwd = 2,
  xlab = "truncation level",
  ylab = "log-likelihood",
  main = "fit"
)
points(sel, ll[sel], pch = 19)
## the level-1 value lies far below the rest; show the path from level 2
plot(
  ks[-1],
  mb[-1],
  type = "l",
  lwd = 2,
  xlab = "truncation level",
  ylab = "mBICV",
  main = "criterion (levels 2 to 49)"
)
points(sel, mb[sel], pch = 19)
points(kmin, mb[kmin], pch = 1, cex = 1.4)
legend(
  "topright",
  bty = "n",
  pch = c(19, 1),
  legend = c("selected by trunc_lvl = NA", "minimum over the path")
)
dev.off()

cat(
  "figures written to figures/; selected level",
  sel,
  "path minimum",
  kmin,
  "\n"
)
sessionInfo()
