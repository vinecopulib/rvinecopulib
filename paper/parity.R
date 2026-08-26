## ---------------------------------------------------------------------------
## Bivariate parity against VineCopula.
##
##   Rscript parity.R
##
## Compares the density, distribution function, both h-functions, both inverse
## h-functions, Kendall's tau, Blomqvist's beta and the tail dependence
## coefficients for every parametric family the two packages share, in every
## admissible rotation, on a fixed grid of the unit square.  Produces the
## summary reported in Table 7 of the manuscript.
##
## Note the argument convention: hbicop(u, 1, .) conditions on the FIRST
## variable and corresponds to VineCopula's BiCopHfunc1.
## ---------------------------------------------------------------------------

suppressMessages({
  library("rvinecopulib")
  library("VineCopula")
})

grid <- as.matrix(expand.grid(seq(0.05, 0.95, by = 0.05),
                              seq(0.05, 0.95, by = 0.05)))

## family -> (VineCopula base code, representative parameters)
spec <- list(
  gaussian = list(1, 0.5),        t      = list(2, c(0.5, 5)),
  clayton  = list(3, 2.0),        gumbel = list(4, 2.0),
  frank    = list(5, 3.0),        joe    = list(6, 2.0),
  bb1      = list(7, c(1, 1.5)),  bb6    = list(8, c(1.5, 1.5)),
  bb7      = list(9, c(1.5, 1)),  bb8    = list(10, c(3, 0.8))
)

## VineCopula codes rotations by an offset; rotationless families keep the base
rot_code <- function(base, rot) {
  if (rot == 0 || base %in% c(1, 2, 5)) return(base)
  c("90" = base + 20, "180" = base + 10, "270" = base + 30)[[as.character(rot)]]
}
## and negates the parameter for the 90 and 270 degree rotations
neg <- function(rot, par) if (rot %in% c(90, 270)) -par else par

out <- NULL
for (fam in names(spec)) {
  base <- spec[[fam]][[1]]
  par <- spec[[fam]][[2]]
  rotations <- if (base %in% c(1, 2, 5)) 0 else c(0, 90, 180, 270)
  for (rot in rotations) {
    rv <- bicop_dist(fam, rot, par)
    p <- neg(rot, par)
    vc <- BiCop(rot_code(base, rot), p[1], if (length(p) > 1) p[2] else 0)
    td <- BiCopPar2TailDep(vc$family, vc$par, vc$par2)
    out <- rbind(out, data.frame(
      family = fam, rotation = rot,
      pdf  = max(abs(dbicop(grid, rv)   - BiCopPDF(grid[, 1], grid[, 2], vc))),
      cdf  = max(abs(pbicop(grid, rv)   - BiCopCDF(grid[, 1], grid[, 2], vc))),
      h1   = max(abs(hbicop(grid, 1, rv) - BiCopHfunc1(grid[, 1], grid[, 2], vc))),
      h2   = max(abs(hbicop(grid, 2, rv) - BiCopHfunc2(grid[, 1], grid[, 2], vc))),
      hi1  = max(abs(hbicop(grid, 1, rv, inverse = TRUE) -
                       BiCopHinv1(grid[, 1], grid[, 2], vc))),
      hi2  = max(abs(hbicop(grid, 2, rv, inverse = TRUE) -
                       BiCopHinv2(grid[, 1], grid[, 2], vc))),
      tau  = abs(par_to_ktau(rv) - BiCopPar2Tau(vc$family, vc$par, vc$par2)),
      beta = abs(blomqvist_beta(rv) - BiCopPar2Beta(vc$family, vc$par, vc$par2)),
      taildep = max(abs(diag(tail_dep(rv)) - c(td$lower, td$upper)))
    ))
  }
}

cat("configurations compared:", nrow(out), "\n\n")
cat("largest absolute deviation by quantity:\n")
print(signif(apply(out[, -(1:2)], 2, max), 3))

cat("\nKendall's tau by family:\n")
print(signif(tapply(out$tau, out$family, max), 3))

## Frank is the outlier.  A high-precision reference settles which
## implementation is responsible.
theta <- 3
debye1 <- integrate(function(t) ifelse(t < 1e-12, 1, t / expm1(t)),
                    0, theta, rel.tol = 1e-12, subdivisions = 2000)$value / theta
reference <- 1 - 4 / theta + 4 * debye1 / theta
cat("\nFrank's tau at theta = 3:\n")
print(c(reference    = reference,
        rvinecopulib = par_to_ktau(bicop_dist("frank", 0, theta)),
        VineCopula   = BiCopPar2Tau(5, theta)))

## Joe's tau has a removable singularity at theta = 2; the limit is 2 - pi^2/6.
cat("\nJoe's tau at theta = 2:\n")
print(c(value = par_to_ktau(bicop_dist("joe", 0, 2)), limit = 2 - pi^2 / 6))

saveRDS(out, "parity.rds")
sessionInfo()
