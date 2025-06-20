#  These functions are based on Files src/library/stats/R/glm.R
#    and src/library/stats/R/vcov.R
#  Part of the R package, https://www.R-project.org
#
#  Modified by Mark W. Donoghoe:
#    31/01/2025 - version of summary.glm that employs Matrix
#    06/06/2025 - version of vcov.glm that employs Matrix

#
#  Copyright (C) 1995-2020 The R Core Team
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

#' @title Summarizing Generalized Linear Model Fits
#' @description Implementation of \code{\link[stats]{summary.glm}} for GLMs fit 
#' using \code{glm2(..., method = "\link{glm.fit2.Matrix}")}. See its documentation
#' for full details. 
#' @usage NULL
#' @method summary glm2Matrix
#' @importFrom stats coef pnorm pt residuals
#' @export

summary.glm2Matrix <- function(object, dispersion = NULL,
                               correlation = FALSE, symbolic.cor = FALSE, ...)
{
  est.disp <- FALSE
  df.r <- object$df.residual
  if(is.null(dispersion)) {	# calculate dispersion if needed
    fam <- object$family
    dispersion <-
      if (!is.null(fam$dispersion) && !is.na(fam$dispersion)) fam$dispersion
    else if(fam$family %in% c("poisson", "binomial"))  1
    else if(df.r > 0) {
      est.disp <- TRUE
      if(any(object$weights==0))
        warning("observations with zero weight not used for calculating dispersion")
      sum((object$weights*object$residuals^2)[object$weights > 0])/ df.r
    } else {
      est.disp <- TRUE
      NaN
    }
  }
  ## calculate scaled and unscaled covariance matrix
  
  aliased <- is.na(coef(object))  # used in print method
  p <- object$rank
  if (p > 0) {
    p1 <- 1L:p
    if (is.null(object$qr)) stop("lm object does not have a proper 'qr' component.\n Rank zero or should not have used lm(..., qr=FALSE).")
    Qr <- object$qr
    ## WATCHIT! doesn't this rely on pivoting not permuting 1L:p? -- that's quaranteed
    coef.p <- object$coefficients[Qr$pivot[p1]]
    # Need to back-permute the vcov matrix based on how Matrix::qr permuted R
    qrpivot <- Qr$qr@q + 1L
    qr.uns <- is.unsorted(qrpivot, strictly = TRUE)
    p2 <- Matrix::invertPerm(qrpivot)
    covmat.perm <- Matrix::chol2inv(Matrix::qr.R(Qr$qr))
    if (qr.uns) covmat.unscaled <- covmat.perm[p2, p2, drop = FALSE]
    else covmat.unscaled <- covmat.perm
    rm(covmat.perm)
    dimnames(covmat.unscaled) <- list(names(coef.p),names(coef.p))
    covmat <- dispersion*covmat.unscaled
    var.cf <- Matrix::diag(covmat)
    
    ## calculate coef table
    
    s.err <- sqrt(var.cf)
    tvalue <- coef.p/s.err
    
    dn <- c("Estimate", "Std. Error")
    if(!est.disp) { # known dispersion
      pvalue <- 2*pnorm(-abs(tvalue))
      coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
      dimnames(coef.table) <- list(names(coef.p),
                                   c(dn, "z value","Pr(>|z|)"))
    } else if(df.r > 0) {
      pvalue <- 2*pt(-abs(tvalue), df.r)
      coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
      dimnames(coef.table) <- list(names(coef.p),
                                   c(dn, "t value","Pr(>|t|)"))
    } else { # df.r == 0
      coef.table <- cbind(coef.p, NaN, NaN, NaN)
      dimnames(coef.table) <- list(names(coef.p),
                                   c(dn, "t value","Pr(>|t|)"))
    }
    df.f <- NCOL(Qr$qr)
  } else {
    coef.table <- matrix(, 0L, 4L)
    dimnames(coef.table) <-
      list(NULL, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    covmat.unscaled <- covmat <- matrix(, 0L, 0L)
    df.f <- length(aliased)
  }
  ## return answer
  
  ## these need not all exist, e.g. na.action.
  keep <- match(c("call","terms","family","deviance", "aic",
                  "contrasts", "df.residual","null.deviance","df.null",
                  "iter", "na.action"), names(object), 0L)
  ans <- c(object[keep],
           list(deviance.resid = residuals(object, type = "deviance"),
                coefficients = coef.table,
                aliased = aliased,
                dispersion = dispersion,
                df = c(object$rank, df.r, df.f),
                cov.unscaled = covmat.unscaled,
                cov.scaled = covmat))
  
  if(correlation && p > 0) {
    dd <- sqrt(Matrix::diag(covmat.unscaled))
    ans$correlation <-
      covmat.unscaled/outer(dd,dd)
    ans$symbolic.cor <- symbolic.cor
  }
  class(ans) <- c("summary.glm2Matrix", "summary.glm")
  return(ans)
}

## Allow for 'dispersion' to be passed down (see the help for vcov)

#' @title Calculate Variance-Covariance Matrix for a Fitted Model Object
#' @description Implementation of \code{\link[stats]{vcov}} for GLMs fit 
#' using \code{glm2(..., method = "\link{glm.fit2.Matrix}")}. See its documentation
#' for full details. 
#' @usage NULL
#' @method vcov glm2Matrix
#' @importFrom stats vcov
#' @export

vcov.glm2Matrix <- function(object, complete = TRUE, ...)
  vcov(summary.glm2Matrix(object, ...), complete = complete)

#' @rdname vcov.glm2Matrix
#' @usage NULL
#' @exportS3Method vcov summary.glm2Matrix

vcov.summary.glm2Matrix <- function(object, complete = TRUE, ...)
  .vcov.aliased.Matrix(object$aliased, object$cov.scaled, complete = complete)

## Augment a vcov - matrix by NA rows & cols when needed:
.vcov.aliased.Matrix <- function(aliased, vc, complete = TRUE) {
  ## Checking for "NA coef": "same" code as in print.summary.lm() in ./lm.R :
  if(complete && NROW(vc) < (P <- length(aliased)) && any(aliased)) {
	  ## add NA rows and columns in vcov
	  cn <- names(aliased)
	  VC <- Matrix::Matrix(NA_real_, P, P, dimnames = list(cn,cn), sparse = TRUE)
	  j <- which(!aliased)
	  VC[j,j] <- vc
	  VC
  } else  # default
	  vc
}