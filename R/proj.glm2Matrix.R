#  This function is based on File src/library/stats/R/proj.R
#  Part of the R package, https://www.R-project.org
#
#  Modified by Mark W. Donoghoe:
#    27/03/2026 - implementation for use with Matrix
#
#  Copyright (C) 1998-2020 The R Core Team
#  Copyright (C) 1998 B. D. Ripley
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

#' Projections of glm2Matrix models
#' 
#' @description An implementation of \code{\link[stats]{proj}}
#' for GLMs fit using \code{glm2(..., method = "\link{glm.fit2.Matrix}")}
#'
#' @param object an object of class \code{glm2Matrix} as returned by \code{\link{glm2}}.
#' @param onedf A logical flag. If \code{TRUE}, a projection is returned for
#' all the columns of the model matrix. If \code{FALSE}, the single-column
#' projections are collapse by terms of the model (as represented in the
#' analysis of variance table).
#' @param unweighted.scale If the fit producing \code{object} used weights, this
#' determines if the projections correspond to weighted or unweighted observations.
#' @param column.order \code{"fit"} uses the sparse QR ordering from
#'   the fitted model. \code{"original"} recomputes the final weighted 
#'   least-squares QR in the original reduced model-matrix column order 
#'   (via \code{Matrix::qr(..., order = 0L)}); this may be much slower for 
#'   large sparse problems.
#' @param ... Swallow and ignore any other arguments.
#' 
#' @return See documentation for \code{\link[stats]{proj}}.
#' 
#' @author Adapted by Mark Donoghoe from the code for \code{proj} and \code{proj.lm}.
#'
#' @seealso \code{\link[stats]{proj}}
#'
#' @importFrom stats model.matrix proj
#' @method proj glm2Matrix
#' @export

proj.glm2Matrix <- function(object, onedf = FALSE, unweighted.scale = FALSE,
                            column.order = c("fit", "original"), ...) {
  
  column.order <- match.arg(column.order)
  
  rank <- object$rank
  rn <- names(object$residuals)
  wt <- object$weights
  used <- if(is.null(wt)) rep.int(TRUE, length(object$residuals)) else (wt > 0)
  
  if (rank > 0) {
    
    # Only get the model matrix if necessary
    Xall <- NULL
    if (!onedf || column.order == "original")
      Xall <- if (!is.null(object$x)) object$x else model.matrix(object)
    
    if (column.order == "fit") {
      
      if (!inherits(object$qr$qr, "sparseQR"))
        stop("argument does not include a sparseQR 'qr' component")
      if (is.null(object$effects))
        stop("argument does not include an 'effects' component")
      qr <- object$qr$qr
      eff <- object$effects
      
      eff.pivot <- object$effects.pivot
      
    } else {
      # Use the base projection (wo column permutation) 
      # rather than the default from Matrix::qr
      keep <- which(!is.na(object$coefficients))
      Xred <- Xall[used, keep, drop = FALSE]
      Xw <- Xred * sqrt(wt[used])
      offset <- object$offset
      if (is.null(offset)) offset <- rep.int(0, length(object$residuals))
      z <- object$linear.predictors - offset + object$residuals
      zw <- z[used] * sqrt(wt[used])
      # qr without pivoting (CAN BE SLOW)
      tol <- min(1e-07, object$control$epsilon/1000)
      qr <- Matrix::qr(Xw, order = 0L, tol = tol)
      rank0 <- Matrix::qr2rankMatrix(qr0, tol = tol)
      if (rank0 != rank)
        stop("column.order = \"original\" produced rank ", rank0,
             " but the fitted model has rank ", rank)
      eff <- Matrix::qr.qty(qr, zw)
      
      dn <- colnames(Xred)
      if (is.null(dn)) dn <- paste0("x", seq_len(ncol(Xred)))

      nmeffects <- c(dn[seq_len(rank0)],
                     rep.int("", length(eff) - rank0))
      names(eff) <- nmeffects
      
      eff.pivot <- keep
      
    }
    
    # Matrix analogue of stats:::proj.default
    RB <- c(eff[seq_len(rank)],
            rep.int(0, nrow(qr) - rank))
    dqr <- qr@Dim
    RBn <- dqr[1L]
    RBncols <- min(dqr)
    RBD <- Matrix::sparseMatrix(i = seq_len(RBncols), j = seq_len(RBncols),
                                x = RB[seq_len(RBncols)], dims = c(RBn, RBncols))
    
    prj <- Matrix::qr.qy(qr, RBD)
    dimnames(prj) <- list(rn[used], names(eff)[seq(ncol(prj))])
    prj.full <- Matrix::Matrix(0, nrow = length(object$residuals), ncol = ncol(prj),
                               dimnames = list(rn, colnames(prj)), sparse = TRUE)
    prj.full[used, ] <- prj
    prj <- prj.full
    
      
    if (onedf) {
      df <- rep.int(1, rank)
      result <- prj
    } else {
      # NB this comes from the stats::proj.lm source code, but 
      # since a glm2Matrix object (like a glm object) does not contain
      # an $assign element, it doesn't do anything
      asgn <- object$assign[eff.pivot]
      uasgn <- unique(asgn)
      nmeffect <- c("(Intercept)",
                    attr(object$terms, "term.labels"))[1 + uasgn]
      nterms <- length(uasgn)
      df <- vector("numeric", nterms)
      result <- Matrix::Matrix(0, nrow = length(object$residuals), 
                               ncol = nterms, sparse = TRUE,
                               dimnames = list(rn, nmeffect))
      for (i in seq_along(uasgn)) {
        select <- (asgn == uasgn[i])
        df[i] <- sum(select)
        result[, i] <- prj[, select, drop = FALSE] %*% rep.int(1, df[i])
      }
      
    }
      
  } else {
    result <- NULL
    df <- NULL
  }
  
  # Only divide rows actually used in the final WLS fit
  if (!is.null(wt) && unweighted.scale && !is.null(result))
    result[used, ] <- result[used, , drop = FALSE] / sqrt(wt[used])
  
  use.wt <- !is.null(wt) && !unweighted.scale
  if (object$df.residual > 0) {
    res <- if (use.wt) object$residuals * sqrt(wt) else object$residuals
    if (!inherits(result, "Matrix")) {
      result <- Matrix::Matrix(res, nrow = length(res), ncol = 1L,
                               dimnames = list(names(res), "Residuals"))
    } else {
      res_Mat <- Matrix::Matrix(res, ncol = 1, sparse = TRUE)
      colnames(res_Mat) <- "Residuals"
      result <- cbind(result, res_Mat)
      rownames(result) <- names(res)
    }
    df <- c(df, object$df.residual)
  }
  
  names(df) <- colnames(result)
  attr(result, "df") <- df
  attr(result, "formula") <- object$call$formula
  attr(result, "onedf") <- onedf
  attr(result, "column.order") <- column.order
  if (!is.null(wt)) attr(result, "unweighted.scale") <- unweighted.scale
  result
  
}