# JSS manuscript: rvinecopulib

Working draft of

> Nagler T, Vatter T. *rvinecopulib: A Fast and Comprehensive Implementation of
> Vine Copula Models in C++ and R.* Submitted to the Journal of Statistical
> Software.

## Files

| File | Purpose |
|---|---|
| `rvinecopulib.tex` | The manuscript (`jss` document class). |
| `rvinecopulib.bib` | Bibliography. Every entry was verified against Crossref or the publisher. |
| `code.R` | Replication script. `RVCL_FAST=1 Rscript code.R` runs the reduced grid. |
| `make_figures.R` | Regenerates `figures/*.pdf`. |
| `parity.R` | Bivariate parity against VineCopula; reproduces Table 7. |
| `build.sh` | Builds the PDF, adding R's texmf tree to `TEXINPUTS`/`BSTINPUTS`. |
| `jsslogo.pdf` | Not tracked. `build.sh` generates a 1pt placeholder if it is missing; replace it with the real logo from the JSS author kit before submitting. It is the sole source of the one remaining overfull-box warning. |

## Building

```sh
./build.sh
```

`jss.cls` and `jss.bst` ship with R under `$(R RHOME)/share/texmf`; `build.sh`
puts that on the TeX search path so no system-wide installation is needed.

## Status

**Complete draft.** 37 pages, compiles with no errors, no undefined
references, no overfull boxes other than the placeholder logo. Every section is
written; both case studies use real data and real fitted numbers; the benchmark
and validation numbers are measured, not estimated.

`R CMD check` on the modified package: **Status: 1 NOTE** (installed size
82 Mb, unavoidable for a vendored header-only C++ library).

`code.R` runs end to end (`RVCL_FAST=1 Rscript code.R`); `make_figures.R`
regenerates all three figures.

### Remaining placeholders

Seven, all mechanical:

| Marker | Count | What it needs |
|---|---|---|
| `[OUTPUT: from code.R]` | 6 | Paste printed output into the `CodeOutput` blocks. |
| `[ACKNOWLEDGMENTS]` | 1 | Author input. |

### Not yet done

* `code.R` has been run end to end in `FAST` mode only; the full path has not
  been timed (you said the one-hour budget does not matter, so this is
  informational).
* Table 8's maximum-likelihood rows stop at `d = 20` — the `d = 50` MLE cell
  would take VineCopula a very long time and adds nothing the τ-inversion row
  does not already show.
* Acknowledgments are still a placeholder.

### Package changes made while drafting

Writing the paper surfaced three defects in the package. All are fixed, with
tests.

1. **`vcov()` / `confint()` did not exist.** The manuscript claimed them before
   they did. Implemented in `R/vcov.R` with methods for `bicop_dist`,
   `vinecop_dist` and `vine_dist`, plus `tests/testthat/test_vcov.R`.

2. **The sandwich estimator was wrong** (found by adversarial review, and the
   most serious of the three). Under the default `step_wise = TRUE`,
   `hessian()` does not return a Hessian: it returns the transpose of the
   Jacobian of the step-wise estimating equations, which is *exactly block
   upper-triangular* and not symmetric. The first implementation symmetrized it
   before inverting, discarding the cross-tree propagation terms, and used
   `A^-1 B A^-1` instead of `A^-1 B A^-T`.

   Verified: `max|H - Jᵀ| = 2.8e-09` while `max|H - J| = 0.28`; the lower
   triangle is exactly zero. Monte Carlo (1200 replications, 4-d Gaussian
   D-vine, n = 400) against the true sampling SD:

   | | tree 1 | | | tree 2 | | tree 3 |
   |---|---|---|---|---|---|---|
   | corrected | −0.7% | +0.9% | +1.8% | +1.4% | −2.4% | +0.1% |
   | old (symmetrized) | +3.9% | +7.0% | +5.3% | −1.1% | −4.2% | −0.4% |

   Coverage of nominal 95% intervals moved from 0.951/0.964/0.961 (over-covering)
   to 0.941/0.948/0.955. `type = "model"` is now refused for `step_wise = TRUE`,
   since the information equality does not hold for the sequential estimator.

3. **`par_to_ktau()` returned `NaN` for the Joe copula at θ = 2.** Joe's
   Kendall tau has a removable singularity there; the implementation evaluated
   the 0/0 quotient directly. It is reachable in practice because
   `ktau_to_par("joe", 2 - pi^2/6)` returns exactly 2, so the round trip
   failed. Found by comparing against VineCopula while building Table 3.

   > **This fix is in a vendored file**
   > (`inst/include/vinecopulib/bicop/implementation/joe.ipp`) and will be lost
   > at the next re-vendor. **It needs to go upstream to vinecopulib.**

   Note also that editing a vendored header does *not* trigger recompilation —
   `src/*.o` must be deleted by hand. This is the same problem upstream fixed
   in vinecopulib #753, which the vendored copy predates.

4. **`margin_info.univariateML()` broke the serialization round trip.** It was
   the only one of the four `univariateML` methods missing the
   `check_univariateML()` guard, and it calls `stats::logLik()`, whose method
   lives in that namespace. Reading a saved model in a fresh session therefore
   failed. One-line fix in `R/margins.R`, regression test in `test_margins.R`.

### What the benchmarks actually showed

The headline was drafted as "orders of magnitude faster". Measurement supports
that only for some operations, and the manuscript now separates them.

**Maximum-likelihood fitting** vs VineCopula 2.6.1, matched settings,
log-likelihoods agreeing to 1e-6: 1.5× (d=5), 1.6× (d=10), 1.7× (d=20) on one
core. Modest, and the paper says so.

**τ-inversion fitting** — where parameter estimation stops dominating and the
parallelized parts (edge weights, per-tree pair-copula fits) take over:

| d | 1 core | 4 cores | 8 cores |
|---|---|---|---|
| 5 | 2.0× | 6.6× | 4.1× |
| 10 | 2.6× | 8.3× | 8.9× |
| 20 | 1.1× | 6.9× | 9.4× |
| 50 | 1.9× | 8.5× | **23.4×** |

VineCopula's `cores` uses socket clusters, so its τ-inversion fits get *slower*
with more cores (d=50: 52 s → 67 s → 124 s) while ours keep scaling
(28 s → 8 s → 5 s). The paper says explicitly that the 23× is not an
algorithmic claim — the single-core column is, and there it is ~2×.

**Evaluation on a fitted model** (speed-up, 1 core / 8 cores):

| Operation | d=5 | d=10 | d=20 | d=50 |
|---|---|---|---|---|
| Density | 2.6 / 6.6 | 2.0 / 8.6 | 2.4 / 6.1 | 2.0 / 6.1 |
| Log-likelihood | 2.2 / 4.5 | 2.1 / 9.1 | 2.5 / 14.0 | 2.3 / 14.5 |
| Rosenblatt | 2.2 / 3.4 | 1.9 / 5.6 | 1.9 / 9.8 | 1.3 / 6.2 |
| Inverse Rosenblatt | 1.6 / 6.7 | 1.7 / 8.8 | 1.7 / 11.0 | 2.3 / 8.0 |
| Simulation | 1.6 / 6.3 | 1.8 / 9.3 | 1.2 / 11.2 | 1.2 / 5.7 |

**Nonparametric evaluation** vs kdecopula: density 36×, h-function 337×,
inverse h-function ~5300×. But the two use different bandwidth rules, so they
implement the same *method* rather than the same estimator, and kdecopula is
18–41% more accurate at τ=0.6. The paper states both, and also makes the point
that speed is not the substance there: kdecopula estimates one bivariate
copula, with no vine, no discrete data and no way to compete against parametric
families in selection.

**Not measured:** memory. `bench`'s `mem_alloc` only tracks R-level allocation,
which flatters us enormously (315 KB vs 382 MB at d=5) because our work happens
in C++. Peak RSS would be the fair measure; the abstract no longer promises
memory benchmarks.

### Second review pass

A second adversarial review of the rewritten sections raised 36 findings and
confirmed **26**, including two critical. All are fixed. The most serious:

* **Fabricated output** (critical) — see below.
* **"The replication material reproduces every number, figure and table"**
  (critical) was false. Table 2 is compiled from package documentation and
  Figures 1–2 are diagrams; neither is script output. The sentence now names
  the three scripts and states what is *not* in them.
* **A sign error in the Joe series expansion I had written** — see below.
* **The TLL accuracy numbers came from 15 replications and were wrong at the
  top end.** Re-measured over 200 replications with Monte Carlo standard
  errors: at τ=0.6 our error exceeds kdecopula's by 17.4% / 23.9% / 28.3%
  (SE ≈ 1.3), not "18% rising to 41%". At τ=0.3 we are marginally *better*, not
  identical. `code.R` now runs that grid.
* **"Four orders of magnitude"** and **"11× to 14×"** were both arithmetic that
  did not survive checking. Replaced with the measured values (the τ-inversion
  log-likelihood gap is 0.05%–0.6%, always in our favour — now stated, with the
  direction).
* **`coverage_study()` returned only proportions**, not the pooled average or
  the paired standard error the paper quotes. It now returns all three.

The rest:

* **A sign error in the Joe series expansion I had written.** The correct
  expansion is `1 - psi'(2) + d[psi'(2)/2 + psi''(2)/4]`; I had a minus on the
  second term. Verified numerically: the correct sign gives ~1e-10 error
  against the closed form, mine gave ~2e-6. The regression test now checks the
  expansion against the closed form at 1e-8, which sits between the two.
* **The Joe round-trip claim was false.** `ktau_to_par("joe", 2 - pi^2/6)`
  returns 1.99999999974, not exactly 2. The NaN was still reachable (θ = 2 is
  an ordinary interior point) but the stated mechanism was wrong. Reworded.
* **The fitted Joe copula has *lower* tail dependence, not upper.** It is
  rotated 180° (survival Joe), τ_lower = 0.59. The corrected reading is more
  interesting: low-exposure policies almost never claim, more reliably than
  high-exposure policies always do.
* **"Multi-core advantage grows with dimension, reaching 11× to 14× at d=20"**
  was wrong twice — the d=20 column spans 6.1 to 14.0, and the advantage falls
  at d=50 for four of five operations. Rewritten with the honest pattern.
* **Atom mass 0.932 in the text vs 0.929 in the figure.** The figure was right
  (the model is fitted on the subsample).
* **§7.4 said "a tenth to a half" of h-functions are skipped**, contradicting
  §3.1, where a D-vine at d=50 skips only 3.9%. Rewritten.
* **The "same log-likelihood" claim hid absolute gaps** of 33 nats at d=20 and
  169 at d=50 under τ-inversion, always in our favour. Now stated in both
  relative and absolute terms, with the caveat that the τ-inversion rows
  compare two implementations of a procedure rather than two computations of
  one number.
* **The socket-cluster explanation was overstated** — VineCopula's MLE fits *do*
  speed up with cores. Narrowed to the τ-inversion case.
* **`code.R` ran 300 coverage replications where the paper reports 1000**, and
  ran only one of the three TLL sample sizes. Both fixed.
* **`code.R` referenced a `parity.R` that did not exist.** Written; it
  reproduces Table 7 exactly.

### Fabricated output — found and removed

The most serious finding of the second pass. When I transcribed
`summary(fit$copula)` into the manuscript I hand-wrote `zi,c` and `zi,d` in the
`var_types` column. **The software cannot produce that.** `simplify_var_types()`
(`R/vine.R:633`) collapses `"zi"` to `"d"` before the copula sees it, and
`correct_var_types()` (`R/tools.R:375`) accepts only `"c"`/`"d"` at the copula
layer. The real output says `d,c`.

The block is now verbatim from the fitted model, and §6.2 gained a paragraph
explaining *why* the column says `d` — the third type lives at the
vine-distribution level, where it selects the left limit; the copula layer only
needs to know whether a variable has atoms.

### Two more claims that were false

* **"which the modified BIC would have truncated away had we asked it to"**
  (§6.2) was wrong twice. mBICV *prefers* keeping trees 2–3, by 34.6 units
  (verified: −6298.80 / −6324.43 / −6333.44 for k = 1/2/3), and
  `conditioning_set` refuses truncation anyway
  (`R/vinecop.R:245`). Replaced with the accurate statement, including that the
  two features cannot currently be combined.
* **The truncation figure's caption promised an open circle marking the
  criterion minimum. No such marker was ever drawn.** Worse, the underlying
  tension was never explained: Table 6 shows the untruncated model with the
  better mBICV (−52247 vs −51810) while `trunc_lvl = NA` returned 10 trees,
  which reads as a broken selector. §6.1 now explains that the search is
  *sequential* — it stops the first time a tree fails to improve the criterion
  and never revisits (`tools_select.ipp:318`,
  `pass.mbicv_trunc >= pass.mbicv`) — so the two rows answer different
  questions. Caption corrected to match the figure.

### Claims corrected after measurement

* **"The masks roughly halve h-function evaluations"** was wrong. Measured: the
  fraction consumed is exactly 1/2 for C-vines, ≈0.57–0.61 for random R-vines,
  and 0.67 rising to 0.96 for D-vines as d grows. The manuscript now reports
  the measurements and notes that D-vines in high dimensions gain almost
  nothing.
* **Family counts.** rvinecopulib has 11 parametric families, not 12;
  VineCopula has 12, one more than us. The table now says so.
* **Three Table 2 cells** understated VineCopula (observation weights, a
  user-supplied tree criterion, and standard errors via `RVineStdError()`).
  Corrected.
* **`tll` derivatives.** The draft said argument derivatives were exact. They
  are not implemented at all.
* **BB8 tail dependence** is zero over the interior of the parameter space.
* **Parametric vs nonparametric was being decided by mBICV**, which is not what
  information criteria are for — they compare parametric candidates with
  commensurable parameter counts. §6.1 now compares *whole models* out of
  sample instead. Result: parametric-itau 16.12, parametric-mle 15.85,
  nonparametric 8.29 (9 trees) / 13.33 (3) / 11.44 (1), independence 0.
  The nonparametric vine has 16,818 effective d.f. against 1000 training
  observations, so it overfits — the expected outcome at this d and n, not a
  defect.
* **`itau` beats `mle`** on these data: 23× faster to fit (1.6 s vs 37 s) *and*
  better out of sample (16.12 vs 15.85). The two-parameter families MLE can
  reach buy in-sample fit that does not survive the split.
* **Automatic truncation stopping at 10 trees** was written up as a discrepancy
  in an earlier draft. It is not — the criterion has flattened by then, and
  where to stop in the flat region is a modeling choice. Reworded.
* **Coverage study** was over-read from 300 replications. Re-run at 1000 with
  the Monte Carlo standard error reported; only the aggregate difference
  (0.027, SE 0.004) is now claimed.

### Verified while drafting
* Every code example in Sections 4 and 5 runs against `rvinecopulib` 1.0.0.1.0
  and is reproduced in `code.R`.
* `vcov()` and `confint()` were **implemented during drafting** (`R/vcov.R`) —
  the manuscript had claimed them before they existed. Wald intervals cover the
  truth; see below.
* Coverage study: nominal 95% intervals attain 0.930 / 0.957 / 0.953 with known
  margins and 0.920 / 0.923 / 0.940 with rank-estimated margins (300
  replications, `n = 1000`, three-dimensional Gaussian D-vine). The shortfall is
  the two-stage margin effect, not an implementation error; it is reported in
  §4.7.


### Known blockers inherited from the release

These are not manuscript problems but they gate submission; see the publication
plan for detail.

1. The vendored `vinecopulib` headers are behind upstream `main`, including
   correctness fixes that change every TLL and discrete fit.
2. Version 1.0.0.1.0 is not on CRAN. JSS requires it.
3. `VineCopula` on CRAN is 2.6.1; the derivative parity comparison in §7.2
   needs 2.6.2, which currently exists only on GitHub.

## Software used for this draft

R 4.3.3; rvinecopulib 1.0.0.1.0 (vinecopulib 1.0.0); VineCopula 2.6.1;
kdecopula 0.9.3; kde1d 1.1.1; wdm 0.2.6; univariateML 1.5.0;
insuranceData 1.0; bench 1.1.4.
