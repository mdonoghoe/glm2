#' Fitting Generalized Linear Models
#' 
#' Fits generalized linear models using the same model specification as
#' \code{glm} in the \bold{stats} package, but with a modified default fitting
#' method. The method provides greater stability for models that may fail to
#' converge using \code{glm}.
#' 
#' \tabular{ll}{ Package: \tab glm2\cr Type: \tab Package\cr Version: \tab
#' 1.2.1\cr License: \tab GPL (>=2)\cr LazyData: \tab true\cr } There are two
#' functions in the package, \code{glm2} and \code{glm.fit2}. The \code{glm2}
#' function fits generalized linear models using the same model specification
#' as \code{glm} in the \bold{stats} package. It is identical to \code{glm}
#' except for minor modifications to change the default fitting method. The
#' \code{glm.fit2} function provides the default fitting method for
#' \code{glm2}. It is identical to \code{glm.fit} in the \bold{stats} package,
#' except for modifications to the computational method that provide more
#' stable convergence. Normally only \code{glm2} would be called directly,
#' although like \code{glm.fit}, \code{glm.fit2} can be called directly. It can
#' also be passed to \code{glm} as an alternative fitting method, using the
#' \code{method} argument. See Marschner (2011) for a discussion of the need
#' for a modified fitting method.
#' 
#' @name glm2-package
#' @docType package
#' @author Ian Marschner (using code from \code{glm} and \code{glm.fit} in the
#' \bold{stats} package)
#' 
#' Maintainer: Mark W. Donoghoe \email{markdonoghoe@@gmail.com}
#' @references Marschner, I.C. (2011) glm2: Fitting generalized linear models
#' with convergence problems. The R Journal, Vol. 3/2, pp.12-15.
#' @keywords package
"_PACKAGE"

#' Horseshoe Crab Data
#' 
#' This data set is derived from Agresti (2007, Table 3.2, pp.76-77). It gives
#' 4 variables for each of 173 female horseshoe crabs. Also provided are two
#' random samples of the data with replacement, which are useful for
#' illustrating the convergence properties of \code{glm} and \code{glm2}.
#' 
#' The variables \code{Dark} and \code{GoodSpine} are derived from the raw
#' data. In the notation of Table 3.2 of Agresti (2007), \code{Dark = yes}
#' corresponds to C>3 and \code{GoodSpine = yes} corresponds to S<3. The two
#' random samples \code{Rep1} and \code{Rep2} can be used to provide random
#' samples with replacement from the full data set. These two random samples
#' are useful for illustrating the convergence properties of \code{glm} and
#' \code{glm2}; see examples in the help documentation for \code{glm2}.
#' 
#' @name crabs
#' @docType data
#' @format A data frame with 173 observations on the following 6 variables:
#' \describe{ \item{list("Satellites")}{number of male partners in addition to
#' the female's primary partner} \item{list("Width")}{width of the female in
#' centimeters} \item{list("Dark")}{a binary factor indicating whether the
#' female has dark coloring (\code{yes} or \code{no})}
#' \item{list("GoodSpine")}{a binary factor indicating whether the female has
#' good spine condition (\code{yes} or \code{no})} \item{list("Rep1")}{a random
#' sample with replacement from 1:173} \item{list("Rep2")}{a second random
#' sample with replacement from 1:173} }
#' @references Agresti, A. (2007) An Introduction to Categorical Data Analysis
#' (2nd ed.). Hoboken, NJ: Wiley.
#' @keywords datasets
NULL

#' Heart Attack Data
#' 
#' This data set is a cross-tabulation of data on 16949 individuals who
#' experienced a heart attack (ASSENT-2 Investigators, 1999).  There are 4
#' categorical factors each at 3 levels, together with the number of patients
#' and the number of deaths for each observed combination of the factors. This
#' data set is useful for illustrating the convergence properties of \code{glm}
#' and \code{glm2}.
#' 
#' The variable \code{AgeGroup} groups the age of the patients as follows: 1 =
#' <65 years, 2 = 65-75 years, 3 = >75 years. The variable \code{Delay} groups
#' the time from heart attack onset to treatment as follows: 1 = <2 hours, 2 =
#' 2-4 hours, 3 = >4 hours. The variable \code{Severity} denotes the severity
#' of the heart attack using the Killip class: 1 = least severe (class I), 2 =
#' middle severity (class II), 3 = most severe (class III or IV). The variable
#' \code{Region} provides the geographical region in which the patients were
#' treated: 1 = Western countries, 2 = Latin America, 3 = Eastern Europe. This
#' data set is useful for illustrating the convergence properties of \code{glm}
#' and \code{glm2}; see examples in the help documentation for \code{glm2}.
#' 
#' @name heart
#' @docType data
#' @format A data frame with 74 observations on the following 6 variables.
#' \describe{ \item{list("Deaths")}{number of deaths}
#' \item{list("Patients")}{number of patients}
#' \item{list("AgeGroup")}{categorization of the age of the patients}
#' \item{list("Severity")}{severity of the heart attack}
#' \item{list("Delay")}{categorization of the time from heart attack to
#' treatment} \item{list("Region")}{geographical region in which the patients
#' were treated} }
#' @references ASSENT-2 Investigators. (1999) Single-bolus tenecteplase
#' compared with front-loaded alteplase in acute myocardial infraction: The
#' ASSENT-2 double-blind randomized trial. Lancet 354, 716-722.
#' @keywords datasets
NULL
