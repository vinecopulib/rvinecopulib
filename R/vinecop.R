#' Vine copula models
#' 
#' Automated fitting or creation of custom vine copula models
#' 
#' @aliases vinecop_dist
#' @inheritParams bicop
#' @param family_set a character vector of families; see [bicop()] for 
#' additional options.
#' @param matrix a quadratic matrix specifying the structure matrix (see 
#'   [check_rvine_matrix()]); for [vinecop_dist()], the dimension must be 
#'   `length(pair_copulas)-1`; for [vinecop()], `matrix = NA` performs
#'   automatic structure selection.
#' @param trunc_lvl the truncation level of the vine copula; `Inf` means no
#'   truncation, `NA` indicates that the truncation level should be selected
#'   automatically by GIC.
#' @param tree_crit the criterion for tree selection, one of `"tau"`, `"rho"`,
#'    `"hoeffd"` for Kendall's \eqn{tau}, Spearman's \eqn{rho}, and Hoeffding's
#'    \eqn{D}, respectively.
#' @param threshold for thresholded vine copulas; `NA` indicates that the 
#'   threshold should be selected automatically by GIC.
#' @param show_trace logical; whether a trace of the fitting progress should be 
#'    printed.
#' @param cores number of cores to use; if more than 1, estimation of pair 
#'    copulas within a tree is done in parallel.
#' 
#' @details
#' `vinecop_dist()` creates a vine copula by specifying a nested list of 
#' `bicop_dist` objects and a quadratic structure matrix. 
#' 
#' `vinecop()` provides automated fitting for vine copula models. 
#' The function inherits the parameters of [bicop()]. 
#' Optionally, a quadratic `matrix` can be used as
#' input to prespecify the vine structure. `tree_crit` describes the
#' criterion for tree selection, one of `"tau"`, `"rho"`, `"hoeffd"` for
#' Kendall's tau, Spearman's rho, and Hoeffding's D, respectively. Additionally, 
#' `threshold` allows to threshold the `tree_crit` and `trunc_lvl` to truncate 
#' the vine copula, with `threshold_sel` and `trunc_lvl_sel` to automatically 
#' select both parameters.
#'
#' @return Objects inheriting from `vinecop_dist` for [vinecop_dist()], and
#' `vinecop` and `vinecop_dist` for [vinecop()].
#'
#' @examples
#' # specify pair-copulas
#' bicop <- bicop_dist("bb1", 90, c(3, 2))
#' pcs <- list(
#'     list(bicop, bicop),  # pair-copulas in first tree 
#'     list(bicop)          # pair-copulas in second tree 
#'  )
#'  
#' # specify R-vine matrix
#' mat <- matrix(c(1, 2, 3, 1, 2, 0, 1, 0, 0), 3, 3) 
#' 
#' # set up vine copula model
#' vc <- vinecop_dist(pcs, mat)
#' 
#' # show model
#' summary(vc)
#' 
#' # get some data
#' u <- sapply(1:3, function(i) runif(50))
#' 
#' # estimate a vine copula model
#' fit <- vinecop(u, "par")
#' fit
#' summary(fit)
#' str(fit, 3)
#' 
#' @export
vinecop <- function(data, family_set = "all", matrix = NA, 
                    par_method = "mle", nonpar_method = "constant",
                    mult = 1, selcrit = "bic", presel = TRUE, 
                    trunc_lvl = Inf, tree_crit = "tau", threshold = 0, 
                    keep_data = TRUE, show_trace = FALSE, cores = 1) {
    # check if families known (w/ partial matching) and expand convenience defs
    family_set <- process_family_set(family_set)
    
    ## pre-process input
    data <- if_vec_to_matrix(data)
    if (any(is.na(matrix)))
        matrix <- as.matrix(0)
    
    ## fit and select copula model
    vinecop <- vinecop_select_cpp(
        data = data, 
        matrix = matrix,
        family_set = family_set,
        par_method = par_method,
        nonpar_method = nonpar_method,
        mult = mult,
        selection_criterion = selcrit,
        preselect_families = presel,
        truncation_level = ifelse(  # Inf cannot be passed to C++
            is.finite(trunc_lvl),
            trunc_lvl, 
            .Machine$integer.max
        ),
        tree_criterion = tree_crit,
        threshold = threshold,
        select_truncation_level = is.na(trunc_lvl),
        select_threshold = is.na(threshold),
        show_trace = show_trace,
        num_threads = cores
    )
    
    ## make all pair-copulas bicop objects
    vinecop$pair_copulas <- lapply(
        vinecop$pair_copulas, 
        function(tree) lapply(tree, as.bicop)
    )

    ## add information about the fit
    vinecop$names <- colnames(data)
    if (keep_data) {
        vinecop$data <- data
    }
    vinecop$controls <- list(
        family_set = family_set,
        par_method = par_method,
        nonpar_method = nonpar_method,
        mult = mult,
        selcrit = selcrit,
        presel = presel,
        trunc_lvl = trunc_lvl,
        tree_crit = tree_crit,
        threshold = threshold
    )
    vinecop$nobs <- nrow(data)
    
    structure(vinecop, class = c("vinecop", "vinecop_dist"))
}

#' @param pair_copulas A nested list of 'bicop_dist' objects, where 
#'    \code{pair_copulas[[t]][[e]]} corresponds to the pair-copula at edge `e` in
#'    tree `t`.
#' @rdname vinecop
#' @export
vinecop_dist <- function(pair_copulas, matrix) {
    # create object
    vinecop <- structure(
        list(pair_copulas = pair_copulas, matrix = matrix),
        class = "vinecop_dist"
    )
    
    # sanity checks
    stopifnot(is.list(pair_copulas))
    pc_lst <- unlist(pair_copulas, recursive = FALSE)
    if (!all(sapply(pc_lst, function(x) inherits(x, "bicop_dist")))) {
        stop("some objects in pair_copulas aren't of class 'bicop_dist'")
    }
    vinecop_check_cpp(vinecop)
    check_rvine_matrix(matrix)
    vinecop$npars <- sum(sapply(pc_lst, function(x) x[["npars"]]))
    
    vinecop
}
