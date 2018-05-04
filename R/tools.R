#' Internal: Turn vector input into a matrix with two columns
#'
#' @param u input data
#'
#' @return either a matrix with two columns, or an error if u is neither a
#' matrix, data.frame, or a length two vector
#'
#' @noRd
if_vec_to_matrix <- function(u) {
    assert_that(is.numeric(u) | is.data.frame(u))
    if (NCOL(u) == 1)
        u <- matrix(u, 1, length(u))
    if (!is.matrix(u))
        u <- as.matrix(u)
    
    u
}

#' Internal: Convert arguments to `bicop_dist` object.
#' @param family the family as passed in function call.
#' @param rotation the rotation as passed in function call.
#' @param parameters the parameters as passed in function call.
#' @return A `bicop_dist` object.
#' @noRd
args2bicop <- function(family, rotation, parameters) {
    if (all(inherits(family, "bicop_dist"))) {
        return(family)
    } else {
        if (missing(rotation))
            rotation <- 0
        if (missing(parameters))
            parameters <- numeric(0)
        assert_that(is.string(family), is.number(rotation), is.numeric(parameters))
        return(bicop_dist(family, rotation, parameters))
    }
}

process_family_set <- function(family_set) {
    family_set <- check_and_match_family_set(family_set)
    expand_family_set(family_set)
}

#' Internal: Expand shortcuts in the familyset.
#' @noRd
expand_family_set <- function(family_set) {
    unique(unlist(lapply(family_set, expand_family)))
}

expand_family <- function(family) {
    switch(
        family,
        "archimedean"   = family_set_archimedean,
        "ellipiltical"  = family_set_elliptical,
        "bbs"           = family_set_bb,
        "oneparametric" = family_set_onepar,
        "twoparametric" = family_set_twopar,
        "parametric"    = family_set_parametric,
        "nonparametric" = family_set_nonparametric,
        "itau"          = family_set_itau,
        "all"           = family_set_all,
        family  # default is no expansion
    )
}

#' Internal: Checks whether all families are known (including partial matching).
#' @noRd
check_and_match_family_set <- function(family_set) {
    matched_fams <- family_set_all_defs[pmatch(family_set, family_set_all_defs)]
    if (any(is.na(matched_fams))) {
        stop("unknown families in family_set: ",
            paste0('"', family_set[is.na(matched_fams)], '"', collapse = ", "))
    }
    matched_fams
}

#' @importFrom stats runif
prep_uniform_data <- function(n, d, U) {
    if (is.null(U)) {
        U <- matrix(runif(n * d), n, d)
    } else {
        assert_that(is.matrix(U), nrow(U) == n)
        if (d == 2) {
            assert_that(ncol(U) == 2)
        } else {
            assert_that(ncol(U) == eval(d))
        }
    }
    U
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
        print(plots[[1]])
        
    } else {
        # Set up the page
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(layout = 
                                              grid::grid.layout(nrow(layout), 
                                                                ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], 
                  vp = grid::viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
        }
    }
}

# Get the depth of a list
depth <- function(this) ifelse(is.list(this), 1L + max(sapply(this, depth)), 0L)

supported_distributions <- c("beta", "cauchy", "chisq", "exp", "f", "gamma", 
                             "lnorm", "norm", "t", "unif", "weibull")

#' @importFrom stats pbeta qbeta qbeta dcauchy pcauchy qcauchy dchisq pchisq
#' @importFrom stats qchisq dexp pexp qexp df pf qf dgamma pgamma qgamma
#' @importFrom stats dlnorm plnorm qlnorm dt pt qt dunif punif qunif
#' @importFrom stats dweibull pweibull qweibull
check_distr <- function(distr) {
    
    ## if provided with a kde1d object, then there is nothing to do
    if (inherits(distr, "kde1d"))
        return(TRUE)
    
    ## basic sanity checks
    if (!is.list(distr))
        return("a distribution should be a kde1d object or a list")
    if (!any(is.element(names(distr), "name")))
        return("a distribution should be a kde1d object or a list with a name element")
    nn <- distr[["name"]]
    if (!is.element(nn, supported_distributions))
        return("the provided name does not belong to supported distributions")
    
    ## check that the provided parameters are consistent with the distribution
    qfun <- get(paste0("q", nn))
    par <- distr[names(distr)!= "name"]
    par$p <- 0.5
    e <- tryCatch(do.call(qfun, par), error = function(e) e)
    if (any(class(e) == "error"))
        return(e$message)
    
    return(TRUE)
}

get_npars_distr <- function(distr) {
    switch(distr$name,
           beta = 2,
           cauchy = 2,
           chisq = ifelse("ncp" %in% names(distr), 2, 1),
           exp = 1,
           f = 3,
           gamma = 2,
           lnorm = 2,
           norm = 2,
           t = ifelse("ncp" %in% names(distr), 3, 2),
           unif = 2,
           weibull = ifelse("scale" %in% names(distr), 2, 1))
}

#' @noRd
#' @importFrom assertthat assert_that on_failure<-
#' @importFrom assertthat is.number is.string is.flag is.scalar
in_set <- function(el, set) {
    all(el %in% set)
}


on_failure(in_set) <- function(call, env) {
    paste0(deparse(call$el), 
           " must be one of {", 
           paste0(eval(call$set, env), collapse = ", "), 
           "}.")
}
