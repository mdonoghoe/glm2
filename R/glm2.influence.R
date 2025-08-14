#  This function is based on File src/library/stats/R/lm.influence.R
#  Part of the R package, https://www.R-project.org
#
#  Modified by Mark W. Donoghoe:
#    31/01/2025 - error for glm.fit2.Matrix
#    14/08/2025 - implementation for use with Matrix
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

### "lm"  *and*	 "glm"	 leave-one-out influence measures

## The following is adapted from John Fox's  "car" :


#' Regression Diagnostics
#'
#' @description An implementation of \code{\link[stats]{influence}}
#' for GLMs fit using \code{glm2(..., method = "\link{glm.fit2.Matrix}")}.
#'
#' @param model an object of class \code{glm2Matrix} as returned by \code{\link{glm2}}.
#' @param do.coef logical indicating if the changed \code{coefficients} are desired.
#' @param ... further arguments passed to or from other methods (ignored).
#'
#' @return See documentation for \code{\link[stats]{influence}}.
#'
#' @seealso \code{\link[stats]{influence}}, and its "See Also" section.
#'
#' @examples
#' ## Analysis of the life-cycle savings data
#' ## given in Belsley, Kuh and Welsch.
#' summary(glm.SR <- glm2(sr ~ pop15 + pop75 + dpi + ddpi,
#'                        data = LifeCycleSavings,
#'                        method = "glm.fit2.Matrix"))
#' utils::str(glmI <- influence(glm.SR))
#'
#' @method influence glm2Matrix
#' @importFrom stats influence
#' @export

influence.glm2Matrix <- function(model, do.coef = TRUE, ...) {
  
  res <- lm.influence.glm2Matrix(model, do.coef = do.coef, ...)
  pRes <- na.omit(residuals(model, type = "pearson"))[model$prior.weights != 0]
  pRes <- naresid(model$na.action, pRes)
  names(res)[names(res) == "wt.res"] <- "dev.res"
  c(res, list(pear.res = pRes))
  
}

#' @rdname influence.glm2Matrix
#' @export

lm.influence.glm2Matrix <- function(model, do.coef = TRUE, ...) {

  # TODO: Improve on this to use Matrix instead of matrix
  wt.res <- weighted.residuals(model)
  e <- na.omit(wt.res)
  is.mlm <- is.matrix(e) # n x q  matrix in the mulitvariate lm case
  if (model$rank == 0) {
    n <- length(wt.res) # drops 0 wt, may drop NAs
    sigma <- sqrt(deviance(model) / df.residual(model))
    res <- list(hat = rep(0, n), coefficients = matrix(0, n, 0),
                sigma = rep(sigma, n))
  } else {
    ## if we have a point with hat = 1, the corresponding e should be
    ## exactly zero.  Protect against returning Inf by forcing this
    e[abs(e) < 100 * .Machine$double.eps * median(abs(e))] <- 0
    mqr <- model$qr
    n <- as.integer(nrow(mqr$qr))
    if (is.na(n)) stop("invalid model QR matrix")
    ## in na.exclude case, omit NAs; also drop 0-weight cases
    if(NROW(e) != n)
      stop("non-NA residual length does not match cases used in fitting")
    do.coef <- as.logical(do.coef)
    tol <- 10 * .Machine$double.eps
    ## original code:
    ## res <- .Call(C_influence, mqr, e, tol)
    res <- influence.Matrix(mqr, e, tol)
    if (do.coef) {
      # NB mqr is QR after removal of linearly dependent columns
      Q <- Matrix::qr.Q(mqr$qr)
      R <- Matrix::qr.R(mqr$qr)
      hat <- res$hat
      invRQtt <- Matrix::t(Matrix::solve(R, Matrix::t(Q)))
      # Need to back-permute this
      mpivot <- mqr$qr@q + 1L
      qr.uns <- is.unsorted(mpivot, strictly = TRUE)
      if (qr.uns) invRQtt <- invRQtt[, Matrix::invertPerm(mpivot)]
      k <- NCOL(Q)
      q <- NCOL(e)
      ## NB: The following relies on recycling: diag(v) %*% A == A * v
      ## so we need a for loop for the mlm case
      res$coefficients <-
        if(is.mlm) {
          cf <- array(0, c(n,k,q))
          for (j in seq_len(q))
            cf[,,j] <- invRQtt * ifelse(hat == 1, 0, e[,j]/(1 - hat))
          cf
        } else
          invRQtt * ifelse(hat == 1, 0, e/(1 - hat))
    }

    drop1d <- function(a) { # more cautious variant of drop(.)
      d <- dim(a)
      if(length(d) == 3L && d[[3L]] == 1L)
        dim(a) <- d[-3L]
      a
    }
    if (is.null(model$na.action)) {
      if (!is.mlm) {  ## drop the 'q=1' array extent (from C)
        res$sigma <- drop(res$sigma)
        if (do.coef)
          res$coefficients <- drop1d(res$coefficients)
      }
    } else {
      hat <- naresid(model$na.action, res$hat)
      hat[is.na(hat)] <- 0      # omitted cases have 0 leverage
      res$hat <- hat
      if (do.coef) {
        coefficients <- naresid.Matrix(model$na.action, res$coefficients)
        coefficients[is.na(coefficients)] <- 0 # omitted cases have 0 change
        res$coefficients <- if(is.mlm) coefficients else drop1d(coefficients)
      }
      sigma <- naresid(model$na.action, res$sigma)
      sigma[is.na(sigma)] <- sqrt(deviance(model)/df.residual(model))
      res$sigma <- if(is.mlm) sigma else drop(sigma)
    }
  }
  res$wt.res <- naresid(model$na.action, e)
  res$hat[res$hat > 1 - 10*.Machine$double.eps] <- 1 # force 1
  names(res$hat) <- names(res$sigma) <- names(res$wt.res)
  if(do.coef) {
	  cf <- coef(model)
	  if(is.mlm) { # coef is 3d array
	    dnr <- dimnames(res$wt.res)
	    dimnames(res$coefficients) <- list(
		    dnr[[1L]],
		    rownames(cf)[!apply(cf, 1L, anyNA)],
		    dnr[[2L]])
	  } else {
	    dimnames(res$coefficients) <- list(names(res$wt.res),
					                               names(cf)[!is.na(cf)])
	  }
  }
  res[c("hat", "coefficients", "sigma", "wt.res")] # ensure order, for backward compatibility and regression tests
}

# R version of src/library/stats/src/influence.c
#   which calls src/library/stats/src/lminfl.f
# Using the Matrix package

#' @keywords internal

influence.Matrix <- function(mqr, e, tol) {

  # Diagonal of hat matrix
  Q <- Matrix::qr.Q(mqr$qr)
  h <- Matrix::rowSums(Q^2)
  h[h >= 1 - tol] <- 1

  # Estimated residual standard deviation
  n <- nrow(mqr$qr)
  k <- Matrix::qr2rankMatrix(mqr$qr)
  
  is.mlm <- is.matrix(e)
  q <- NCOL(e)
  e <- matrix(e, nrow = n, ncol = q)

  denom <- (n - k - 1)
  e2sum <- colSums(e^2)
  sigma <- matrix(sqrt(e2sum / denom), nrow = n, ncol = q, byrow = TRUE)
  
  hk <- h < 1
  # sigma[i,j] = sqrt((e2sum[j] - e[i,j]^2 / (1 - h[i])) / denom)
  s1 <- sweep(e[hk,,drop=FALSE]^2, 1, 1 - h[hk], FUN = "/")
  s2 <- sweep(-s1, 2, e2sum, FUN = "+")
  sigma[hk,] <- sqrt(s2 / denom)

  list(hat = h, sigma = sigma)

}

#  This function is based on File src/library/stats/R/nafns.R
#  Part of the R package, https://www.R-project.org
#
#  Modified by Mark W. Donoghoe:
#    14/08/2025 - applied to x of type Matrix
#

#' @keywords internal

naresid.Matrix <- function(omit, x, ...) {

  if (class(omit) != "exclude")
    return(x)

  if (length(omit) == 0 || !is.numeric(omit))
    stop("invalid argument 'omit'")
  if (!inherits(x, "Matrix"))
    stop("x is not a Matrix")

  n <- NROW(x)
  keep <- rep.int(NA, n + length(omit))
  keep[-omit] <- 1L:n
  x <- x[keep, , drop = FALSE]
  temp <- rownames(x)
  if (length(temp)) {
    temp[omit] <- names(omit)
    rownames(x) <- temp
  }

  x

}

#' @rdname dfbetas.glm2Matrix
#' @method dfbeta glm2Matrix
#' @export

dfbeta.glm2Matrix <- function(model, infl = lm.influence.glm2Matrix(model, do.coef=TRUE), ...) {
  b <- infl$coefficients
  mlm <- is.matrix(wr <- infl$wt.res)
  if (!mlm) dimnames(b) <- list(names(wr), variable.names(model))
  b
}

#' Regression Deletion Diagnostics (partially implemented)
#'
#' @description A suite of functions computing regression (leave-one-out deletion)
#' diagnostics for GLMs fit using \code{glm2(..., method = "\link{glm.fit2.Matrix}")}. 
#' Currently only \code{dfbetas} and \code{dfbeta} are implemented.
#'
#' @return See documentation for \code{\link[stats]{influence.measures}}.
#'
#' @seealso \code{\link[stats]{influence.measures}}.
#'
#' @method dfbetas glm2Matrix
#' @importFrom stats dfbetas lm.influence
#' @export

dfbetas.glm2Matrix <- function(model, infl = lm.influence.glm2Matrix(model, do.coef=TRUE), ...) {
  
  qrm <- model$qr
  xxi <- Matrix::chol2inv(Matrix::qr.R(qrm$qr))
  db <- dfbeta(model, infl)
  if (length(dim(db)) == 3L) db <- aperm(db, c(1L, 3:2))
  # Might need to fix this up if sigma is a Matrix
  db / outer(infl$sigma, sqrt(Matrix::diag(xxi)))
  
}
