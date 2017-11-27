\name{glm.fit2}
\alias{glm.fit2}
\title{
Generalized Linear Models Fitting Method
}
\description{
\code{glm.fit2} is identical to \code{glm.fit} in the \bold{stats} package, except for a modification to the computational method that provides improved convergence properties. It is the default fitting method for \code{glm2} and can also be used as an alternative fitting method for \code{glm}, instead of the default method \code{glm.fit}. 
}
\usage{
glm.fit2(x, y, weights = rep(1, nobs), start = NULL, 
    etastart = NULL, mustart = NULL, offset = rep(0, nobs), 
    family = gaussian(), control = list(), intercept = TRUE,
    singular.ok = TRUE)
}
\arguments{
  \item{x}{
as for \code{\link{glm.fit}}
}
  \item{y}{
as for \code{\link{glm.fit}}
}
  \item{weights}{
as for \code{\link{glm.fit}}
}
  \item{start}{
as for \code{\link{glm.fit}}
}
  \item{etastart}{
as for \code{\link{glm.fit}}
}
  \item{mustart}{
as for \code{\link{glm.fit}}
}
  \item{offset}{
as for \code{\link{glm.fit}}
}
  \item{family}{
as for \code{\link{glm.fit}}
}
  \item{control}{
as for \code{\link{glm.fit}}
}
  \item{intercept}{
as for \code{\link{glm.fit}}
}
  \item{singular.ok}{
as for \code{\link{glm.fit}}
}
}
\details{
\code{glm.fit2} is a modified version of \code{glm.fit} in the \bold{stats} package. The computational method in \code{glm.fit2} uses a stricter form of step-halving to deal with numerical instability in the iteratively reweighted least squares algorithm. Whereas \code{glm.fit} uses step-halving to correct divergence and parameter space violations, \code{glm.fit2} additionally uses step-halving to force the model deviance to decrease at each iteration, which improves the convergence properties. Like \code{glm.fit}, \code{glm.fit2} usually would not be called directly. Instead, it is called by \code{glm2} as the default fitting method. Like \code{glm.fit}, it is possible to call \code{glm.fit2} directly if the response vector and design matrix have already been calculated, in which case it may be more efficient than calling \code{glm2}. It can also be passed to \code{glm} in the \bold{stats} package, using the \code{method} argument, providing an alternative to the default fitting method \code{glm.fit}. See Marschner (2011) for a discussion of the need for a modified fitting method.
}
\value{
The value returned by \code{glm.fit2} has exactly the same documentation as the value returned by \code{glm.fit}.
}
\author{
\code{glm.fit2} uses the code from \code{glm.fit}, whose authors are listed in the help documentation for the \bold{stats} package. Modifications to this code were made by Ian Marschner.
}
\references{
Marschner, I.C. (2011) glm2: Fitting generalized linear models with convergence problems. The R Journal, Vol. 3/2, pp.12-15.
}
\seealso{
\code{\link{glm.fit}}
}
\examples{
library(glm2)

#--- logistic regression null model ---------------#
# (intercept estimate = log(0.75/0.25) = 1.098612)

y <- c(1,1,1,0)
x <- rep(1,4)

#--- divergence of glm.fit to infinite estimate ---#
fit1 <- glm.fit(x,y, family=binomial(link="logit"),start=-1.81)
fit2 <- glm.fit2(x,y, family=binomial(link="logit"),start=-1.81)
print.noquote(c(fit1$coef,fit2$coef))

}

