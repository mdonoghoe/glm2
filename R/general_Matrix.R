#  These functions are based on Files
#    src/library/stats/R/aov.R
#    src/library/base/R/kappa.R
#  Part of the R package, https://www.R-project.org
#
#  Modified by Mark W. Donoghoe:
#    19/09/2025 - implementation for use with Matrix
#

#  Copyright (C) 1995-2012 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/

#' Final Aliases in a Model using Matrix-based QR decomposition
#'
#' This function identifies complete aliasing (perfect linear dependencies) among predictors
#' in a GLM model using a QR decomposition from the Matrix package. It mimics the logic of
#' base R's \code{alias()} but works with sparse or dense Matrix objects.
#'
#' @param object A fitted GLM object with Matrix-based QR decomposition, for example one fit using \code{glm2(..., method = "glm.fit2.Matrix")}.
#' @param complete Should information on complete aliasing be included?
#' @param partial Should information on partial aliasing be included?
#' @param partial.pattern Should partial aliasing be presented in a schematic way? If this is done, the results are presented in a more compact way, usually giving the deciles of the coefficients.
#' @param ... Additional arguments (currently unused).
#' 
#' @details
#' See documentation for \code{\link[stats]{alias}} for full details.
#' 
#' \code{complete = TRUE} requires the package \CRANpkg{MASS} (recommended by \CRANpkg{stats}) to be installed.
#'
#' @return A list containing aliasing information, including coefficients expressing dependent
#'   columns as linear combinations of independent ones. See documentation for \code{\link[stats]{alias}} for full details.
#'   
#' @note
#' Although some parts of this function employ Matrix-based functions, densification (\code{as.matrix}) is unavoidable and can be slow or memory intensive for large models.
#' 
#' @author Uses the code from \code{\link[stats]{alias}}, whose authors are
#' listed in its documentation. Modifications to this code were made by Mark W. Donoghoe.
#'   
#' @seealso \code{\link[stats]{alias}}
#' 
#' @examples
#' 
#' data(heart)
#' 
#' start.p <- sum(heart$Deaths)/sum(heart$Patients)
#' 
#' # Model with a linearly dependent column
#' fit1 <- glm2(cbind(Deaths,Patients-Deaths) ~ factor(AgeGroup) + factor(Severity) +
#'  I(AgeGroup > 1), data = heart, family = binomial(link="log"), 
#'  start = c(log(start.p), rep(0, 5)), method = "glm.fit2.Matrix", x = TRUE)
#' 
#' alias(fit1, partial = TRUE, partial.pattern = TRUE)
#'
#' @export
#' @method alias glm2Matrix

alias.glm2Matrix <- function(object, complete = TRUE, partial = FALSE,
                             partial.pattern = FALSE, ...) {

  CompPatt <- function(x, ...) {
    x[abs(x) < 1e-6] <- 0
    MASS::fractions(x)
  }
  PartPatt <- function(x) {
    z <- zapsmall(x) != 0
    if(any(z)) {
      xx <- abs(signif(x[z], 2))
      ll <- length(unique(xx))
      if(ll > 10L) xx <- cut(xx, 9L) else if(ll == 1L) x[] <- 1
      x[z] <- paste0(ifelse(x[z] > 0, " ", "-"), xx)
    }
    x[!z] <- ""
    collabs <- colnames(x)
    collabs <- if(length(collabs))
      abbreviate(sub(".", "", collabs, fixed=TRUE), 3L)
    else 1L:ncol(x)
    colnames(x) <- collabs
    class(x) <- "mtable"
    x
  }
  
  Model <- object$terms
  attributes(Model) <- NULL
  value <- list(Model = Model)
  R <- Matrix::qr.R(object$qr_full$qr) # Use the Matrix function
  d <- dim(R)
  rank <- object$rank
  p <- d[2L]
  if(complete) {                      # full rank, no aliasing
    value$Complete <-
      if(is.null(p) || rank == p) NULL else {
        p1 <- 1L:rank
        # Separate the linearly independent and dependent coefficients
        indep <- object$qr$pivot[p1]
        dep <- object$qr$pivot[-p1]
        X <- R[indep, indep]
        Y <-  R[indep, dep, drop = FALSE]
        # pivot as appropriate
        piv <- Matrix::invertPerm(object$qr$qr@q + 1)
        # Use the Matrix function
        beta12 <- Matrix::qr.coef(Matrix::qr(X), Y)[piv, , drop = FALSE]
        # dimnames(beta12) <- list(dn[p1], dn[ -p1])
        # Unfortunately MASS::fractions (used in CompPatt) needs a matrix
        CompPatt(as.matrix(Matrix::t(beta12)))
      }
  }
  if(partial) {
    ## We only want one aspect of the summary, which we know to be reliable
    tmp <- suppressWarnings(summary.glm2Matrix(object)$cov.unscaled)
    # Use the Matrix function
    ses <- sqrt(Matrix::diag(tmp))
    beta11 <- tmp /outer(ses, ses)
    beta11[row(beta11) >= col(beta11)] <- 0
    beta11[abs(beta11) < 1e-6] <- 0
    # Unfortunately PartPatt needs a matrix
    if(all(beta11 == 0)) beta11 <- NULL
    else if(partial.pattern) beta11 <- PartPatt(as.matrix(beta11))
    value$Partial <- beta11
  }
  class(value) <- "listof"
  value

}

#' Condition Number for glm2Matrix Fits
#'
#' Computes the condition number of the design matrix for a
#' \code{glm2Matrix} fit (a \code{glm} object fit using 
#' \code{method = "glm.fit2Matrix"}).
#'
#' @param z A fitted \code{glm2Matrix} object.
#' @param ... Additional arguments passed to internal methods.
#'   See \code{\link[base]{kappa}} for details, and the Details below for
#'   a note about the \code{exact} argument.
#'
#' @details
#' This method is based on the QR decomposition used internally to fit
#' the model.
#'
#' It calls \code{\link[base]{kappa}} on the dense R matrix, and therefore
#' does not provide any efficiency benefits for large design matrices over
#' a model fit using \code{\link[stats]{glm}} (or \code{glm2} without \code{Matrix}).
#'
#' Additionally, when \code{exact = FALSE} (the default), the returned value 
#' may not agree with that returned by \code{\link[base]{kappa.lm}} applied to
#' an equivalent dense fit. This is because the sparse QR decomposition may
#' have employed column pivoting so that the R matrix corresponds to a 
#' column-permuted design matrix. This can be checked by inspecting
#' \code{z$qr_full$qr@q}.
#'
#' @return The condition number, \eqn{\kappa}, or an approximation if \code{exact = FALSE}.
#'
#' @seealso \code{\link[base]{kappa}}, \code{\link[Matrix]{sparseQR-class}}
#'
#' @export
#' @method kappa glm2Matrix

kappa.glm2Matrix <- function(z, ...) .kappa_qr_Matrix(z$qr_full, ...)

.kappa_qr_Matrix <- function(z, ...) {
  
  stopifnot(length(d <- dim(qr <- z$qr)) == 2L)
  R <- Matrix::qr.R(qr)
  .kappa_tri(as.matrix(R), ...)
  
}
