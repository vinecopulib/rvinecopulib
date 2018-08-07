#' Inverse Rosenblatt transform
#' 
#' The inverse Rosenblatt transform computes conditional quantiles and can be 
#' used simulate from a stochastic model, see *Details*.
#' 
#' @param u matrix of numbers between 0 and 1.
#' @param model a model object; classes currently supported are 
#'    `bicop_dist()`, `vinecop_dist()`, and `vine_dist()`.
#' @param cores if `>1`, computation is parallelized over `cores` batches (rows
#'    of `u`).
#' 
#' @details 
#' The Rosenblatt transform (Rosenblatt, 1952) \eqn{U = T(V)} of a random vector
#' \eqn{V = (V_1,\ldots,V_d) ~ F} is defined as
#' \deqn{ 
#'   U_1 = V_1, U_2 = F(V_2|V_1), \ldots, U_d =F(V_d|V_1,\ldots,V_{d-1}), 
#' } 
#' where \eqn{F(v_k|v_1,\ldots,v_{k-1})} is the conditional distribution of 
#' \eqn{V_k} given \eqn{V_1 \ldots, V_{k-1}, k = 2,\ldots,d}. The vector \eqn{U}
#' are then independent standard uniform variables. The inverse operation 
#' \deqn{ 
#'   V_1 = U_1, V_2 = F^{-1}(U_2|U_1), \ldots, V_d =F^{-1}(U_d|U_1,\ldots,U_{d-1}), 
#' } 
#' can can be used to simulate from a distribution. For any copula \eqn{F}, if 
#' \eqn{U} is a vector of independent random variables, \eqn{V = T^{-1}(U)} has 
#' distribution \eqn{F}.
#' 
#' @examples 
#' # simulate data
#' x <- replicate(3, rnorm(200))
#' u <- replicate(3, runif(200))
#' 
#' # estimate a vine distribution model
#' fit <- vine(x, copula_controls = list(family_set = "par"))
#' 
#' # inverse rosenblatt transform for vine distribution
#' pairs(inverse_rosenblatt(u, fit))
#' 
#' # inverse rosenblatt transform for vine copula
#' pairs(inverse_rosenblatt(u, fit$copula))
#' 
#' # inverse rosenblatt transform for a bivariate copula
#' plot(inverse_rosenblatt(u[, 1:2], bicop_dist("clayton", 0, 3)))
#' @export
inverse_rosenblatt <- function(u, model, cores = 1) {
    assert_that(
        is.matrix(u) | is.data.frame(u), 
        inherits(model, c("bicop_dist", "vinecop_dist", "vine_dist")),
        is.number(cores)
    )
    col_names <- colnames(u)
    
    if (inherits(model, "bicop_dist")) {
        assert_that(ncol(u) == 2)
        u <- cbind(bicop_hinv2_cpp(u, model), u[, 2])
    } else if (inherits(model, "vinecop_dist")) {
        u <- vinecop_inverse_rosenblatt_cpp(u, model, cores)
    } else {
        u <- vinecop_inverse_rosenblatt_cpp(u, model$copula, cores)
        # prepare marginals if only one is specified
        if (!inherits(vine, "vine") & depth(model$margins) == 1) 
            model$margins <- replicate(dim(model)[1], model$margins, simplify = FALSE)
        u <- dpq_marg(u, model, "q")
    }
    colnames(u) <- col_names
    
    u
}

