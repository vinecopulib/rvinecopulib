% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vinecop_methods.R
\name{vinecop_predict_and_fitted}
\alias{vinecop_predict_and_fitted}
\alias{predict.vinecop}
\alias{fitted.vinecop}
\title{Predictions and fitted values for a vine copula model}
\usage{
\method{predict}{vinecop}(object, newdata, what = "pdf", n_mc = 10^4, cores = 1, ...)

\method{fitted}{vinecop}(object, what = "pdf", n_mc = 10^4, cores = 1, ...)
}
\arguments{
\item{object}{a \code{vinecop} object.}

\item{newdata}{points where the fit shall be evaluated.}

\item{what}{what to predict, either \code{"pdf"} or \code{"cdf"}.}

\item{n_mc}{number of samples used for quasi Monte Carlo integration when
\code{what = "cdf"}.}

\item{cores}{number of cores to use; if larger than one, computations are
done in parallel on \code{cores} batches.}

\item{...}{unused.}
}
\value{
\code{fitted()} and \code{predict()} have return values similar to \code{\link[=dvinecop]{dvinecop()}}
and \code{\link[=pvinecop]{pvinecop()}}.
}
\description{
Predictions of the density and distribution function
for a vine copula model.
}
\details{
\code{fitted()} can only be called if the model was fit with the
\code{keep_data = TRUE} option.
\subsection{Discrete variables}{

When at least one variable is discrete, two types of
"observations" are required in \code{newdata}: the first \eqn{n \; x \; d} block
contains realizations of \eqn{F_{X_j}(X_j)}. The second \eqn{n \; x \; d}
block contains realizations of \eqn{F_{X_j}(X_j^-)}. The minus indicates a
left-sided limit of the cdf. For, e.g., an integer-valued variable, it holds
\eqn{F_{X_j}(X_j^-) = F_{X_j}(X_j - 1)}. For continuous variables the left
limit and the cdf itself coincide. Respective columns can be omitted in the
second block.
}
}
\examples{
u <- sapply(1:5, function(i) runif(50))
fit <- vinecop(u, family = "par", keep_data = TRUE)
all.equal(predict(fit, u), fitted(fit), check.environment = FALSE)
}
