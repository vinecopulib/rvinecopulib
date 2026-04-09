# install.packages("../rvinecopulib_0.7.3.1.0.999.tar.gz")

library(rvinecopulib)
library(copula)
set.seed(0)

rho_true <- 0.6
cop <- normalCopula(param = rho_true)
U <- rCopula(2000, cop)
X <- cbind(qnbinom(U[, 1], size = 2, mu = 300),
           qnbinom(U[, 2], size = 1, mu = 50))

# Compute pseudo-observations
FX <- cbind(pnbinom(X[, 1], size = 2, mu = 300),
            pnbinom(X[, 2], size = 1, mu = 50))
FXm <- cbind(pnbinom(X[, 1] - 1, size = 2, mu = 300),
             pnbinom(X[, 2] - 1, size = 1, mu = 50))
U_disc <- cbind(FX, FXm)

data.frame(
  version = packageVersion("rvinecopulib")
)


vinecop_fit <- rvinecopulib::vinecop(
  data = U_disc,
  var_types = c("d", "d"),
  family_set = "gaussian",
  par_method = "mle"
)
summary(vinecop_fit)

vinecop_fit2 <- rvinecopulib::vinecop(
  data = U_disc,
  var_types = c("d", "d"),
  family_set = "gaussian",
  par_method = "mle",
  structure = vinecop_fit$structure
)
summary(vinecop_fit2)

bicop_fit <- bicop(
  data = U_disc,
  var_types = c("d", "d"),
  family_set = "gaussian",
  par_method = "mle"
)

data.frame(
  vine1 = vinecop_fit$pair_copula[[1]][[1]]$parameters,
  vine2 = vinecop_fit2$pair_copula[[1]][[1]]$parameters,
  bicop = bicop_fit$parameters
)

