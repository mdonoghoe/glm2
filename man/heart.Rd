% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glm2-package.R
\docType{data}
\name{heart}
\alias{heart}
\title{Heart Attack Data}
\format{
A data frame with 74 observations on the following 6 variables.
\describe{ 
\item{Deaths}{number of deaths}
\item{Patients}{number of patients}
\item{AgeGroup}{categorization of the age of the patients}
\item{Severity}{severity of the heart attack}
\item{Delay}{categorization of the time from heart attack to
treatment}
\item{Region}{geographical region in which the patients
were treated} 
}
}
\description{
This data set is a cross-tabulation of data on 16949 individuals who
experienced a heart attack (ASSENT-2 Investigators, 1999).  There are 4
categorical factors each at 3 levels, together with the number of patients
and the number of deaths for each observed combination of the factors. This
data set is useful for illustrating the convergence properties of \code{glm}
and \code{glm2}.
}
\details{
The variable \code{AgeGroup} groups the age of the patients as follows: 1 =
<65 years, 2 = 65-75 years, 3 = >75 years. The variable \code{Delay} groups
the time from heart attack onset to treatment as follows: 1 = <2 hours, 2 =
2-4 hours, 3 = >4 hours. The variable \code{Severity} denotes the severity
of the heart attack using the Killip class: 1 = least severe (class I), 2 =
middle severity (class II), 3 = most severe (class III or IV). The variable
\code{Region} provides the geographical region in which the patients were
treated: 1 = Western countries, 2 = Latin America, 3 = Eastern Europe. This
data set is useful for illustrating the convergence properties of \code{glm}
and \code{glm2}; see examples in the help documentation for \code{glm2}.
}
\references{
ASSENT-2 Investigators. (1999) Single-bolus tenecteplase
compared with front-loaded alteplase in acute myocardial infraction: The
ASSENT-2 double-blind randomized trial. Lancet 354, 716-722.
}
\keyword{datasets}
