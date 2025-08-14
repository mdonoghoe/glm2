#  This function is based on File src/library/stats/R/lm.influence.R
#  Part of the R package, https://www.R-project.org
#
#  Modified by Mark W. Donoghoe:
#    31/01/2025 - error for glm.fit2.Matrix
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


#' Regression Diagnostics (not yet implemented)
#'
#' @description This method is a placeholder. The \code{\link[stats]{influence}}
#' method for GLMs fit using \code{glm2(..., method = "\link{glm.fit2.Matrix}")} 
#' is not yet implemented.
#' 
#' @usage NULL
#'
#' @return An error indicating the method is not implemented.
#' @method influence glm2Matrix
#' @importFrom stats influence
#' @export

influence.glm2Matrix <- function(model, do.coef = TRUE, ...) {
  
  #stop("influence measures not implemented for glm.fit2.Matrix")
  res <- lm.influence.Matrix(model, do.coef = do.coef, ...)
  pRes <- na.omit(residuals(model, type = "pearson"))[model$prior.weights != 0]
  pRes <- naresid(model$na.action, pRes)
  names(res)[names(res) == "wt.res"] <- "dev.res"
  c(res, list(pear.res = pRes))
  
}

lm.influence.Matrix <- function(model, do.coef = TRUE) {

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

# Adapted from stats:::naresid.exclude
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

#' Regression Deletion Diagnostics (not yet implemented)
#'
#' @description This method is a placeholder. The \code{\link[stats]{dfbetas}}
#' method for GLMs fit using \code{glm2(..., method = "\link{glm.fit2.Matrix}")} 
#' is not yet implemented.
#' 
#' @usage NULL
#'
#' @return An error indicating the method is not implemented.
#' @method dfbetas glm2Matrix
#' @importFrom stats dfbetas lm.influence
#' @export

dfbetas.glm2Matrix <- function (model, infl = lm.influence(model, do.coef=TRUE), ...)
{
  
  stop("dfbetas not implemented for glm.fit2.Matrix")
  
}