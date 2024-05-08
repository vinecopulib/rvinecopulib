library(rvinecopulib)
library(truncnorm)

n <- 4000

u <- rbicop(n, "clayton", 0, 3)
z <- rnorm(n)

k <- 5

plot(u)

x <- cbind(
  qbinom(u[, 1], k, pnorm(z)),
  # qbinom(u[, 2], k, pnorm(0.5 * z))
  u[, 2]
)

ut1 <- cbind(
  pbinom(x[, 1], k, pnorm(z)),
  # pbinom(x[, 2], k, pnorm(0.5 * z))
  u[, 2]
)

# plot(ut1)


ut2 <- cbind(
  pbinom(x[, 1] - 1, k, pnorm(z)),
  # pbinom(x[, 2] - 1, k, pnorm(0.5 * z))
  u[, 2]
)

w <- runif(n)
ut <- ut1
ut[, 1] <- w * ut1[, 1] + (1 - w) * ut2[, 1]
w <- runif(n)
ut[, 2] <- w * ut1[, 2] + (1 - w) * ut2[, 2]


uu <- rvinecopulib:::find_latent_sample_cpp(cbind(ut1, ut2), 0.1, 10)
plot(qnorm(uu))
points(qnorm(u), col = 3)

bicop(cbind(ut1, ut2), var_types = c("d", "d"), family = "tll")


plot(ut)

contour(bicop(ut, family = "tll"))
contour(bicop(u, family = "tll"))
contour(bicop(uu, family = "tll"))

ugrid <- pnorm(seq(-3.25, 3.25, l = 30))
ugrid <- as.matrix(expand.grid(ugrid, ugrid)[, 2:1])


fit <- fit_k <- bicop(ut, family = "tll", mult = 0.5, nonpar = "constant")
contour(fit)

b <- 0.01

for (k in 1:10) {
  chat <- dbicop(ugrid, fit_k)
  for (m in 1:nrow(ugrid)) {
    ix <- which(
      ugrid[m, 1] <= ut1[, 1] + b & ugrid[m, 2] <= ut1[, 2] + b &
        ugrid[m, 1] >= ut2[, 1] - b &  ugrid[m, 2] >= ut2[, 2] - b
    )
    n_sel <- length(ix)
    if (n_sel == 0) {
      chat[m] <- 0
      next
    }
    x_mc <- cbind(
      # rtruncnorm(5 * n_sel, qnorm(ut2[ix, 1]), qnorm(ut1[ix, 1]), qnorm(ugrid[m, 1]), b),
      # rtruncnorm(5 * n_sel, qnorm(ut2[ix, 2]), qnorm(ut1[ix, 2]), qnorm(ugrid[m, 2]), b),
      rnorm(5 * n_sel, qnorm(ugrid[m, 1]), b),
      rnorm(5 * n_sel, qnorm(ugrid[m, 2]), b)
    )
    ps <- sapply(ix, \(i) {
      mean(dbicop(cbind(runif(10, ut2[i, 1], ut1[i, 1]),
                        runif(10, ut2[i, 2], ut1[i, 2])),
                  fit_k))
    })
    mean(dbicop(pnorm(x_mc), fit_k) / rep(ps, each = 5))
    chat[m] <- mean(dbicop(pnorm(x_mc), fit_k))
  }
  fit_k$parameters <- matrix(chat, 30, 30)
  contour(fit_k)
}
#
# plot(chat, dbicop(ugrid, fit))
# fit_k <- fit
# fit_k$parameters <- matrix(chat, 30, 30)
#
# plot(chat, dbicop(ugrid, fit))
#
# lattice::wireframe(chat ~ ugrid[, 1] + ugrid[, 2])
#
# v <- ut
#
# for (k in 1:5) {
#
#   fit <- bicop(v, family = "tll", nonpar_method = "constant", mult = 0.1)
#
#   v1 <- runif(n, ut2[, 1], ut1[, 1])
#   v2 <- runif(n, hbicop(cbind(v1, ut2[, 2]), 1, fit), hbicop(cbind(v1, ut1[, 2]), 1, fit))
#   v[, 1] <- v1
#   v[, 2] <- hbicop(cbind(v1, v2), 1, fit, inverse = TRUE)
#   v <- pseudo_obs(v)
#
#   contour(bicop(v, family = "tll"))
#   Sys.sleep(0.1)
# }


xt1 <- qnorm(ut1)
xt2 <- qnorm(ut2)
xt <- qnorm(ut)
x <- qnorm(u)

plot(xt)
points(x, col = 2)

b <- 0.1

for (k in 1:k) {
  for (i in 1:nrow(x)) {
    ix <- which(xt[, 1] <= xt1[i, 1] + b & xt[, 1] >= xt2[i, 1] - b &
                  xt[, 2] <= xt1[i, 2] + b & xt[, 2] >= xt2[i, 2] - b)
    if (length(ix) == 0)
      next
    ix <- sample(ix, 1)
    # xt[i, 1] <- truncnorm::rtruncnorm(1, xt2[i, 1], xt1[i, 1], xt[ix, 1], b)
    # xt[i, 2] <- truncnorm::rtruncnorm(1, xt2[i, 2], xt1[i, 2], xt[ix, 2], b)
    xt[i, 1] <- rnorm(1, xt[ix, 1], b)
    xt[i, 2] <- rnorm(1, xt[ix, 2], b)
  }
  plot(xt)
  points(x, col = 2)

  Sys.sleep(0.1)
  xt <- qnorm(pseudo_obs(xt))
  contour(bicop(pseudo_obs(xt), family = "tll"))
}
# contour(bicop(ut, family = "tll"))

# contour(bicop(u, family = "tll"))





