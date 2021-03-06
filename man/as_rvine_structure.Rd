% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as_rvine_structure.R
\name{as_rvine_structure}
\alias{as_rvine_structure}
\alias{as_rvine_matrix}
\alias{as_rvine_structure.rvine_structure}
\alias{as_rvine_matrix.rvine_structure}
\alias{as_rvine_structure.list}
\alias{as_rvine_matrix.list}
\alias{as_rvine_structure.rvine_matrix}
\alias{as_rvine_matrix.rvine_matrix}
\alias{as_rvine_structure.matrix}
\alias{as_rvine_matrix.matrix}
\title{Coerce various kind of objects to R-vine structures and matrices}
\usage{
as_rvine_structure(x, ...)

as_rvine_matrix(x, ...)

\method{as_rvine_structure}{rvine_structure}(x, ..., validate = FALSE)

\method{as_rvine_matrix}{rvine_structure}(x, ..., validate = FALSE)

\method{as_rvine_structure}{list}(x, ..., is_natural_order = FALSE)

\method{as_rvine_matrix}{list}(x, ..., is_natural_order = FALSE)

\method{as_rvine_structure}{rvine_matrix}(x, ..., validate = FALSE)

\method{as_rvine_matrix}{rvine_matrix}(x, ..., validate = FALSE)

\method{as_rvine_structure}{matrix}(x, ..., validate = TRUE)

\method{as_rvine_matrix}{matrix}(x, ..., validate = TRUE)
}
\arguments{
\item{x}{An object of class \code{rvine_structure}, \code{rvine_matrix}, \code{matrix} or
\code{list} that can be coerced into an R-vine structure or R-vine matrix
(see \emph{Details}).}

\item{...}{Other arguments passed on to individual methods.}

\item{validate}{When `TRUE``, verifies that the input is a valid
rvine-structure (see \emph{Details}). You may want to suppress this when you
know that you already have a valid structure and you want to save some time,
or to explicitly enable it if you have a structure that you want to re-check.}

\item{is_natural_order}{A flag indicating whether the \code{struct_array} element
of \code{x} is assumed to be provided in natural order already (a structure is in
natural order if the anti-diagonal is 1, .., d from bottom left to top
right).}
}
\value{
Either an object of class \code{rvine_structure} or of class
\code{rvine_matrix} (see \code{\link[=rvine_structure]{rvine_structure()}} or \code{\link[=rvine_matrix]{rvine_matrix()}}).
}
\description{
\code{as_rvine_structure} and \code{as_rvine_matrix} are new S3 generics allowing to
coerce objects into R-vine structures and matrices (see \code{\link[=rvine_structure]{rvine_structure()}}
and \code{\link[=rvine_matrix]{rvine_matrix()}}).
}
\details{
The coercion to \code{rvine_structure} and \code{rvine_matrix} can be applied to
different kind of objects Currently, \code{rvine_structure}, \code{rvine_matrix},
\code{matrix} and \code{list} are supported.

For \code{as_rvine_structure}:
\itemize{
\item \code{rvine_structure} : the main use case is to re-check an object via
\code{validate = TRUE}.
\item \code{rvine_matrix} and \code{matrix} : allow to coerce matrices into R-vine
structures (see \code{\link[=rvine_structure]{rvine_structure()}} for more details). The main difference
between \code{rvine_matrix} and \code{matrix} is the nature of the validity
checks.
\item \code{list} : must contain named elements \code{order} and \code{struct_array} to be
coerced into an R-vine structure (see \code{\link[=rvine_structure]{rvine_structure()}} for more details).
}

For \code{as_rvine_matrix}:
\itemize{
\item \code{rvine_structure} : allow to coerce an \code{rvine_structure} into an
R-vine matrix (useful e.g. for printing).
\item \code{rvine_matrix}: similar to \code{as_rvine_structure} for \code{rvine_structure},
the main use case is to re-check an object via \code{validate = TRUE}.
\item \code{matrix} : allow to coerce matrices into R-vine
matrices (mainly by checking that the matrix defines a valid
R-vine, see \code{\link[=rvine_matrix]{rvine_matrix()}} for more details).
\item \code{list} : must contain named elements \code{order} and \code{struct_array} to be
coerced into an R-vine matrix (see \code{\link[=rvine_structure]{rvine_structure()}} for more details).
}
}
\examples{
# R-vine structures can be constructed from the order vector and struct_array
rvine_structure(order = 1:4, struct_array = list(
  c(4, 4, 4),
  c(3, 3),
  2
))

# ... or a similar list can be coerced into an R-vine structure
as_rvine_structure(list(order = 1:4, struct_array = list(
  c(4, 4, 4),
  c(3, 3),
  2
)))

# similarly, standard matrices can be coerced into R-vine structures
mat <- matrix(c(4, 3, 2, 1, 4, 3, 2, 0, 4, 3, 0, 0, 4, 0, 0, 0), 4, 4)
as_rvine_structure(mat)

# or truncate and construct the structure
mat[3, 1] <- 0
as_rvine_structure(mat)

# throws an error
mat[3, 1] <- 5
try(as_rvine_structure(mat))
}
\seealso{
rvine_structure rvine_matrix
}
