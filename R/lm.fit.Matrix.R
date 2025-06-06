#  This function is based on File src/library/stats/R/lm.R
#  Part of the R package, https://www.R-project.org
#
#  Modified by Ian C. Marschner, 2011-2017
#  Modified by Mark W. Donoghoe:
#    25/05/2018 - accept singular.ok argument
#    31/01/2025 - version that employs Matrix


#
#  Copyright (C) 1995-2024 The R Core Team
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


lm.fit.Matrix <- function (x, y, offset = NULL, method = "qr", tol = 1e-07, singular.ok = FALSE, 
          ...) 
{
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package \"Matrix\" must be installed to use method = \"glm.fit2.Matrix\".",
         call. = FALSE)
  }
  if (singular.ok)
    stop("'singular.ok' must be set to FALSE")
  if (is.null(n <- nrow(x))) 
    stop("'x' must be a matrix")
  if (n == 0L) 
    stop("0 (non-NA) cases")
  p <- ncol(x)
  if (p == 0L) {
    ## oops, null model
    return(list(coefficients = numeric(), residuals = y, 
                fitted.values = 0 * y, rank = 0, df.residual = length(y)))
  }
  ny <- NCOL(y)
  ## treat one-col matrix as vector
  if (is.matrix(y) && ny == 1) 
    y <- drop(y)
  if (!is.null(offset)) 
    y <- y - offset
  if (NROW(y) != n) 
    stop("incompatible dimensions")
  if (method != "qr") 
    warning(gettextf("method = '%s' is not supported. Using 'qr'", 
                     method), domain = NA)
  chkDots(...)
  # Original code calls Fortran code here
  # NB order = 0L can be very slow
  # xfit <- Matrix::qr(x, order = 0L, tol = tol)
  xfit <- Matrix::qr(x, tol = tol)
  z <- list(qr = xfit, 
            coefficients = Matrix::qr.coef(xfit, y),
            residuals = Matrix::qr.resid(xfit, y),
            effects = Matrix::qr.qty(xfit, y),
            rank = Matrix::qr2rankMatrix(xfit, tol = tol)
           )

  if (!singular.ok && z$rank < p) stop("singular fit encountered")
  coef <- z$coefficients
  ## careful here: the rank might be 0
  r1 <- seq_len(z$rank)
  dn <- if (is.null(colnames(x))) paste0("x", 1L:p) else colnames(x)
  nmeffects <- c(dn[r1], rep.int("", n - z$rank))
  r2 <- if (z$rank < p) 
    (z$rank + 1L):p
  else integer()
  if (is.matrix(y)) {
    coef[r2, ] <- NA
    dimnames(coef) <- list(dn, colnames(y))
    dimnames(z$effects) <- list(nmeffects, colnames(y))
  }
  else {
    coef[r2] <- NA
    names(coef) <- dn
    names(z$effects) <- nmeffects
  }
  z$coefficients <- coef
  r1 <- y - z$residuals
  if (!is.null(offset)) 
    r1 <- r1 + offset
  c(z[c("coefficients", "residuals", "effects", "rank")], 
    list(fitted.values = r1, assign = attr(x, "assign"), 
         qr = z$qr, df.residual = n - z$rank))
}
