# JSS manuscript: rvinecopulib

Working draft of

> Nagler T, Vatter T. *rvinecopulib: A Fast and Comprehensive Implementation of
> Vine Copula Models in C++ and R.* To be submitted to the Journal of
> Statistical Software.

## Files

| File | Purpose |
|---|---|
| `rvinecopulib.tex` | The manuscript (`jss` document class). |
| `rvinecopulib.bib` | Bibliography. Every entry was verified against Crossref or the publisher. |
| `code.R` | Replication script. `RVCL_FAST=1 Rscript code.R` runs a reduced grid. |
| `parity.R` | Bivariate parity against VineCopula; reproduces the validation table. |
| `make_figures.R` | Regenerates `figures/*.pdf`. |
| `build.sh` | Builds the PDF, adding R's texmf tree to `TEXINPUTS`/`BSTINPUTS`. |
| `jsslogo.pdf` | Not tracked. `build.sh` generates a placeholder if it is missing; replace it with the logo from the JSS author kit before submitting. |

## Building

```sh
./build.sh
```

`jss.cls` and `jss.bst` ship with R under `$(R RHOME)/share/texmf`; `build.sh`
puts that on the TeX search path. The dozens of "Overfull \hbox ... while
\output is active" warnings come from the page-number box of that `jss.cls`
and not from the manuscript.

## Status

Draft under revision after a full review (October 2026). Before submission:

* Fill the acknowledgments.
* Release VineCopula 2.6.2 to CRAN (the multicore timings need it; CRAN has
  2.6.1) and rvinecopulib with the backend fix below.
* Re-run `code.R` in full on the benchmark machine and paste its outputs;
  produce `code.html` with `knitr::spin("code.R")`, which JSS asks for.
* Replace the placeholder logo.

## Replication

`code.R` reproduces every number, table and code example in the manuscript
except Figures 1 and 2 (diagrams) and the capability table (compiled from
package documentation). It ends with `sessionInfo()`. `make_figures.R` draws
the data figures from the same fits, and `parity.R` the bivariate parity
table. Every simulation is seeded.

The non-timing results were computed under R 4.3.3 on an Intel Core Ultra 9
285H laptop (Ubuntu 24.04). They are identical to a run under R 4.6.1 on the
benchmark machine apart from optimizer noise in the seventh significant digit,
which the manuscript never prints.

The timings of Section 7.3 were measured on an idle AMD Ryzen Threadripper
2990WX (32 cores, 125 GB, Ubuntu 22.04, reference BLAS, R 4.6.1). They are
minima over repetitions: contention can only add time, so the minimum
estimates the uncontended cost and discards warm-up; medians of three moved by
up to 37% between runs on the same machine, minima by a few percent. Do not
quote ratios measured on a laptop with mixed performance and efficiency cores:
single-core cells moved by up to 2x between runs on the laptop above. The full
run takes roughly an hour on the Threadripper; `RVCL_FAST=1` skips the timing
grid and shrinks the simulation grids.

Comparator versions: VineCopula 2.6.2 (contains the three parallelization
fixes described in the manuscript; install from GitHub until it reaches
CRAN), kdecopula 0.9.3, kde1d 1.1.1, univariateML 1.5.0, insuranceData 1.0,
qrmdata, bench 1.1.4.

## Package changes carried on this branch

* `inst/include/vinecopulib/vinecop/implementation/tools_select.ipp`:
  `get_num_non_indeps_of_tree()` counted independence copulas, inverting the
  prior of the modified BIC during automatic truncation and thresholding
  (vinecopulib#805). The vendored copy is fixed here; it must be re-vendored
  once upstream releases the fix. Editing a vendored header does not trigger
  recompilation: delete `src/*.o` or install with `--preclean`.

## Known limitations that the manuscript states

* `kde1d` reports no effective parameter count for a zero-inflated fit, so a
  zero-inflated `kde1d` margin can only be compared with other candidates by
  log-likelihood; the insurance example uses a user-defined family for that
  variable.
* vinecopulib 1.0.1 (October 2026) changes the TLL distribution function,
  effective degrees of freedom and discrete TLL fits. Re-vendoring it will move
  the nonparametric rows of the held-out table and the kdecopula accuracy
  comparison.
* Only the sandwich covariance is offered, and only for copula objects; see
  Section 4.7 of the manuscript.

The drafting log that this file used to hold (benchmark fairness checks,
corrected claims, review passes) is in the git history of `paper/README.md`.
