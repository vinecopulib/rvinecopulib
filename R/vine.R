#' Vine copula models
#' 
#' Automated fitting or creation of custom vine copula models
#' 
#' @aliases vine_dist
#' @param data a matrix or data.frame.
#' @param marg_controls a list with arguments to be passed to 
#' [kde1d::kde1d()]. Currently, there can be 
#'   * `mult` numeric; all bandwidhts for marginal kernel density estimation
#'   are multiplied with \code{mult_1d}. Defaults to `log(1 + d)` where `d` is
#'   the number of variables after applying [cctools::expand_as_numeric()].
#'   * `xmin` numeric vector of length d; see [kde1d::kde1d()].
#'   * `xmax` numeric vector of length d; see [kde1d::kde1d()].
#'   * `bw` numeric vector of length d; see [kde1d::kde1d()].
#' @param cop_controls a list with arguments to be passed to [vinecop()].
#' 
#' @details
#' `vine_dist()` creates a vine copula by specifying the marg, a nested list 
#' of `bicop_dist` objects and a quadratic structure matrix. 
#' 
#' `vine()` provides automated fitting for vine copula models. 
#' `marg_controls` is a list with the same parameters as 
#' [kde1d::kde1d()] (except for `x`). `cop_controls` is a list 
#' with the same parameters as [vinecop()] (except for `data`). 
#'
#' @return Objects inheriting from `vine_dist` for [vine_dist()], and
#' `vine` and `vine_dist` for [vine()].
#'
#' @examples
#' # TODO
#' @importFrom kde1d kde1d dkde1d pkde1d qkde1d
#' @importFrom cctools cont_conv expand_vec
#' @export
vine <- function(data, 
                 marg_controls = list(mult = NULL, 
                                      xmin = -Inf, 
                                      xmax = Inf, 
                                      bw = NULL), 
                 cop_controls = list(family_set = "all", 
                                     matrix = NA, 
                                     par_method = "mle", 
                                     nonpar_method = "constant",
                                     mult = 1, 
                                     selcrit = "bic", 
                                     psi0 = 0.9, 
                                     presel = TRUE, 
                                     trunc_lvl = Inf, 
                                     tree_crit = "tau", 
                                     threshold = 0, 
                                     keep_data = TRUE, 
                                     show_trace = FALSE, 
                                     cores = 1)) {
    
    data_cc <- cont_conv(data)
    if (NCOL(data_cc) == 1)
        stop("data must be multivariate.")
    d <- ncol(data_cc)
    
    ## check that the correct arguments are there
    marg_controls_names <- c("mult", "xmin", "xmax", "bw")
    if (!is.list(marg_controls) | !setequal(names(marg_controls), 
                                            marg_controls_names)) {
        msg <- "marg controls should be a list with elements 'mult', 'xmin', 
        'xmax', and 'bw'."
        stop(msg)
    }
    
    ## expand the required arguments and compute default mult if needed
    marg_controls_names <- names(marg_controls)
    marg_controls <- sapply(seq_along(marg_controls), function(j) {
        par <- marg_controls[[j]]
        if (names(marg_controls)[j] == "mult") {
            if (is.null(par)) 
                par <- log(1 + ncol(data_cc))
        } else {
            par <- expand_vec(par, data)
        }
        return(par)
    })
    names(marg_controls) <- marg_controls_names
    
    ## estimation of the marginals
    vine <- list()
    i_disc <- attr(data_cc, "i_disc")
    theta <- attr(data_cc, "theta")
    vine$marg <- lapply(1:d, function(k) 
        kde1d(data_cc[, k], 
              xmin = marg_controls$xmin[k], 
              xmax = marg_controls$xmax[k],
              bw = marg_controls$bw[k],
              mult = marg_controls$mult,
              bw_min = ifelse(k %in% i_disc, 0.5 - theta, 0)))
    vine$marg_controls <- marg_controls
    
    ## estimation of the R-vine copula (only if d > 1)
    if (d > 1) {
        
        ## transform to copula data
        cop_controls$data <- sapply(1:d, function(k) pkde1d(data_cc[, k],
                                                            vine$marg[[k]]))
        
        ## to avoid saving copula data
        keep_data <- cop_controls$keep_data
        cop_controls$keep_data <- FALSE
        
        ## estimate the copula
        vine$cop  <- do.call(vinecop, cop_controls)
        
        ## to potentially save the data on the standard scale
        cop_controls$keep_data <- keep_data
    }
    
    ## add information about the fit
    vine$names <- colnames(data)
    if (keep_data) {
        vine$data <- data_cc
    }
    vine$cop_controls <- cop_controls[-which(names(cop_controls) == "data")]
    vine$nobs <- nrow(data)
    
    ## create and return object
    structure(vine, class = c("vine", "vine_dist"))
}

#' @param marg A list with with each element containing the specification of a 
#' marginal [stats:Distributions](distribution). Each marginal specification 
#' should be a list with containing at least the name and optionally the 
#' parameters, e.g. `list(list(name = "norm"), list(name = "norm", mu = 1), list(name = "beta", shape1 = 1, shape2 = 1))`.
#' Note that parameters that have no default values have to be provided. 
#' Furthermore, if `marg` has length one, it will be recycled for every component.
#' @param pair_copulas A nested list of 'bicop_dist' objects, where 
#'    \code{pair_copulas[[t]][[e]]} corresponds to the pair-copula at edge `e` in
#'    tree `t`.
#' @rdname vine
#' @export
vine_dist <- function(marg, pair_copulas, matrix) {
    
    # sanity checks for the marg
    if (!(length(marg) %in% c(1, ncol(matrix))))
        stop("marg should have length 1 or ncol(matrix)")
    stopifnot(is.list(marg))
    if (depth(marg) == 1) {
        check_marg <- check_distr(marg)
    } else {
        check_marg <- lapply(marg, check_distr)
    }
    is_ok <- sapply(check_marg, isTRUE)
    if (!all(is_ok)) {
        msg <- "Some objects in marg aren't properly defined.\n"
        msg <- c(msg, paste0("margin ", seq_along(check_marg)[!is_ok], " : ",
                             unlist(check_marg[!is_ok]), ".", sep = "\n"))
        stop(msg)
    }
    
    # create the vinecop object
    vinecop <- vinecop(pair_copulas, matrix)
    
    # create object
    structure(list(marg = marg, vinecop = vinecop), class = "vine_dist")
}
