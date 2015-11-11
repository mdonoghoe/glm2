\name{crabs}
\alias{crabs}
\docType{data}
\title{
Horseshoe Crab Data
}
\description{
This data set is derived from Agresti (2007, Table 3.2, pp.76-77). It gives 4 variables for each of 173 female horseshoe crabs. Also provided are two random samples of the data with replacement, which are useful for illustrating the convergence properties of \code{glm} and \code{glm2}.    
}
\usage{data(crabs)}
\format{
  A data frame with 173 observations on the following 6 variables:
  \describe{
    \item{\code{Satellites}}{number of male partners in addition to the female's primary partner}
    \item{\code{Width}}{width of the female in centimeters}
    \item{\code{Dark}}{a binary factor indicating whether the female has dark coloring (\code{yes} or \code{no})}
    \item{\code{GoodSpine}}{a binary factor indicating whether the female has good spine condition (\code{yes} or \code{no})}
    \item{\code{Rep1}}{a random sample with replacement from 1:173}
    \item{\code{Rep2}}{a second random sample with replacement from 1:173} 
  }
}
\details{
The variables \code{Dark} and \code{GoodSpine} are derived from the raw data. In the notation of Table 3.2 of Agresti (2007), \code{Dark = yes} corresponds to C>3 and \code{GoodSpine = yes} corresponds to S<3. The two random samples \code{Rep1} and \code{Rep2} can be used to provide random samples with replacement from the full data set. These two random samples are useful for illustrating the convergence properties of \code{glm} and \code{glm2}; see examples in the help documentation for \code{glm2}.
}
\references{
Agresti, A. (2007) An Introduction to Categorical Data Analysis (2nd ed.). Hoboken, NJ: Wiley.
}
\keyword{datasets}
