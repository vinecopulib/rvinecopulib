#' Fit vine copula models
#' 
#' Automated fitting for vine copula models
#' 
#' @inheritParams bicop_fit
#' @param matrix an R-vine matrix specifying the structure matrix, or `NA` for
#'   automatic structure selection (default).
#' @param trunc_lvl the truncation level of the vine copula.
#' @param tree_crit the criterion for tree selection, one of `"tau"`, `"rho"`,
#'    `"hoeffd"` for Kendall's \eqn{tau}, Spearman's \eqn{rho}, and Hoeffding's
#'    \eqn{D}, respectively.
#' @param threshold for thresholded vine copulas.
#' 
#' @details
#' In addition, to the families in `bicop_dist()`, the following convenience
#' defintion can be used (and combined):
#' \describe{
#' \item{"all"}{all families}.
#' \item{"parametric"}{parametric families.}
#' \item{"nonparametric"}{nonparametric families.}
#' \item{"archimedean"}{archimedean families.}
#' \item{"elliptical"}{elliptical families.}
#' \item{"bbs"}{BB families.}
#' \item{"oneparametric "}{one parameter families.}
#' \item{"twoparametric "}{two parameter families.}
#' }
#' Partial matching is activated. For example, you can write `"nonpar"` instead 
#' of the full name `"nonparametric"`.
#'
#' @return An object inherting from `bicop_fit` and `bicop_dist`.
#'
#' @examples
#' u <- sapply(1:3, function(i) runif(50))
#' fit <- vinecop_select(u, "par")
#' str(fit, 3)
#' 
#' @export
vinecop_select <- function(data, family_set = "all", matrix = NA, method = "mle", 
                        mult = 1, selcrit = "bic", presel = TRUE, 
                        trunc_lvl = Inf, tree_crit = "tau", threshold = 0, 
                        keep_data = TRUE) {
    
    ## family_set can only use standard family names in cpp
    family_set <- family_set_all_defs[pmatch(family_set, family_set_all_defs)]
    family_set <- expand_family_set(family_set)
    
    ## pre-process input
    data <- if_vec_to_matrix(data)
    if (any(is.na(matrix)))
        matrix <- as.matrix(0)
    
    ## fit and select copula model
    vinecop <- vinecop_select_cpp(
        data = data, 
        matrix = matrix,
        family_set = family_set,
        method = method,
        mult = mult,
        selection_criterion = selcrit,
        preselect_families = presel,
        truncation_level = ifelse(  # Inf cannot be passed to C++
            is.finite(trunc_lvl),
            trunc_lvl, 
            .Machine$integer.max
        ),
        tree_criterion = tree_crit,
        threshold = threshold
    )
    
    ## make all pair-copulas bicop_fit objects
    vinecop$pair_copulas <- lapply(
        vinecop$pair_copulas, 
        function(tree) lapply(tree, as.bicop_fit)
    )
    
    ## add information about the fit
    if (keep_data) {
        vinecop$data <- data
    }
    vinecop$controls <- list(
        family_set = family_set,
        method = method,
        mult = mult,
        selcrit = selcrit,
        presel = presel,
        trunc_lvl = trunc_lvl,
        tree_crit = tree_crit,
        threshold = threshold
    )
    
    structure(vinecop, class = c("vinecop_fit", "vinecop_dist"))
}
