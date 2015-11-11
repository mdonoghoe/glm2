\name{glm2-package}
\alias{glm2-package}
\docType{package}
\title{
Fitting Generalized Linear Models
}
\description{
Fits generalized linear models using the same model specification as \code{glm} in the \bold{stats} package, but with a modified default fitting method. The method provides greater stability for models that may fail to converge using \code{glm}.
}
\details{
\tabular{ll}{
Package: \tab glm2\cr
Type: \tab Package\cr
Version: \tab 1.1.2\cr
Date: \tab 2014-05-07\cr
License: \tab GPL (>=2)\cr
LazyLoad: \tab yes\cr
}
There are two functions in the package, \code{glm2} and \code{glm.fit2}. The \code{glm2} function fits generalized linear models using the same model specification as \code{glm} in the \bold{stats} package. It is identical to \code{glm} except for minor modifications to change the default fitting method. The \code{glm.fit2} function provides the default fitting method for \code{glm2}. It is identical to \code{glm.fit} in the \bold{stats} package, except for modifications to the computational method that provide more stable convergence. Normally only \code{glm2} would be called directly, although like \code{glm.fit}, \code{glm.fit2} can be called directly. It can also be passed to \code{glm} as an alternative fitting method, using the \code{method} argument. See Marschner (2011) for a discussion of the need for a modified fitting method.
}
\author{
Ian Marschner (using code from \code{glm} and \code{glm.fit} in the \bold{stats} package)

Maintainer: Ian Marschner <ian.marschner@mq.edu.au>
}
\references{
Marschner, I.C. (2011) glm2: Fitting generalized linear models with convergence problems. The R Journal, Vol. 3/2, pp.12-15.
}
\keyword{package}
