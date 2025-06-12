#  These functions are based on Files src/library/stats/R/predict.glm.R
#      and src/library/stats/R/lm.R
#  Part of the R package, https://www.R-project.org
#
#  Modified by Mark W. Donoghoe:
#    31/01/2025 - version that employs Matrix, if needed

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

predict.glm2Matrix <-
  function(object, newdata = NULL, type = c("link", "response", "terms"),
           se.fit = FALSE, dispersion = NULL, terms = NULL,
           na.action = na.pass, ...)
  {
    ## 1998/06/23 KH:  predict.lm() now merged with the version in lm.R
    
    type <- match.arg(type)
    
    # In this case we have to do things differently for Matrix
    if (se.fit || type == "terms" || !missing(newdata)) {
      
      na.act <- object$na.action
      object$na.action <- NULL # kill this for predict.lm calls
      
      if (!se.fit) {
        ## No standard errors
        if(missing(newdata)) {
          pred <- switch(type,
                         link = object$linear.predictors,
                         response = object$fitted.values,
                         terms = predict.lm.Matrix(object,  se.fit = se.fit,
                                                   scale = 1, type = "terms", terms = terms)
          )
          if(!is.null(na.act)) pred <- napredict(na.act, pred)
        } else {
          pred <- predict.lm.Matrix(object, newdata, se.fit, scale = 1,
                                    type = if(type == "link") "response" else type,
                                    terms = terms, na.action = na.action)
          switch(type,
                 response = {pred <- family(object)$linkinv(pred)},
                 link = , terms = )
        }
      } else {
        ## summary.survreg has no ... argument.
        if(inherits(object, "survreg")) dispersion <- 1.
        if(is.null(dispersion) || dispersion == 0)
          dispersion <- summary(object, dispersion=dispersion)$dispersion
        residual.scale <- as.vector(sqrt(dispersion))
        # Trick to get around the default value of newdata = NULL
        if (missing(newdata)) newdata <- substitute()
        pred <- predict.lm.Matrix(object, newdata, se.fit, scale = residual.scale,
                                  type = if(type == "link") "response" else type,
                                  terms = terms, na.action = na.action)
        fit <- pred$fit
        se.fit <- pred$se.fit
        switch(type,
               response = {
                 se.fit <- se.fit * abs(family(object)$mu.eta(fit))
                 fit <- family(object)$linkinv(fit)
               },
               link = , terms = )
        if( missing(newdata) && !is.null(na.act) ) {
          fit <- napredict(na.act, fit)
          se.fit <- napredict(na.act, se.fit)
        }
        pred <- list(fit = fit, se.fit = se.fit, residual.scale = residual.scale)
      }
      
      
    } else {
      # Trick to get around the default value of newdata = NULL
      if (missing(newdata)) newdata <- substitute()
      pred <- predict.glm(object, newdata, type, se.fit, dispersion, terms,
                            na.action, ...)
    }
    
    pred

  }

predict.lm.Matrix <-
  function(object, newdata, se.fit = FALSE, scale = NULL, df = Inf,
           interval = c("none", "confidence", "prediction"),
           level = .95,  type = c("response", "terms"),
           terms = NULL, na.action = na.pass, pred.var = res.var/weights,
           weights = 1,
           rankdeficient = c("warnif", "simple", "non-estim", "NA", "NAwarn"),
           tol = 1e-6, verbose = FALSE,
           ...)
  {
    
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      stop("Package \"Matrix\" must be installed to use method = \"predict.lm.Matrix\".",
           call. = FALSE)
    }
    tt <- terms(object)
    if(!inherits(object, "lm"))
      warning("calling predict.lm(<fake-lm-object>) ...")
    type <- match.arg(type)
    missRdef <- missing(rankdeficient)
    rankdeficient <- match.arg(rankdeficient)
    noData <- (missing(newdata) || is.null(newdata))
    if(noData) {
      mm <- X <- Matrix::sparse.model.matrix(object)
      mmDone <- TRUE
      offset <- object$offset
    }
    else {
      Terms <- delete.response(tt)
      m <- model.frame(Terms, newdata, na.action = na.action,
                       xlev = object$xlevels)
      if(!is.null(cl <- attr(Terms, "dataClasses"))) .checkMFClasses(cl, m)
      X <- Matrix::sparse.model.matrix(Terms, m, contrasts.arg = object$contrasts)
      if(type != "terms") {
        offset <- model.offset(m)
        if(!is.null(addO <- object$call$offset)) {
          addO <- eval(addO, newdata, environment(tt))
          offset <- if(length(offset)) offset + addO else addO
        }
      }
      mmDone <- FALSE
    }
    n <- length(object$residuals) # NROW(qr(object)$qr)
    p <- object$rank
    p1 <- seq_len(p)
    piv <- if(p) (qrX <- object$qr)$pivot[p1]
    hasNonest <- (p < ncol(X) && !noData)
    nonest <- integer(0L)
    ### NB: Q[p1,] %*% X[,piv] = R[p1,p1]
    if(hasNonest) {
      msg1 <- gettext("prediction from rank-deficient fit")
      if(rankdeficient == "simple") {
        warning(gettextf('%s; consider predict(., rankdeficient="NA")', msg1), domain=NA)
      } else { # rankdeficient more than "simple"
        if(verbose) message("lower-rank qr: determining non-estimable cases")
        stopifnot(is.numeric(tol), tol > 0)
        #if(!p) qrX <- object$qr_full
        qrXf <- object$qr_full
        tR <- Matrix::t(Matrix::qr.R(qrXf$qr))
        pp <- nrow(tR)
        if(verbose) cat(sprintf("  n=%d, p=%d < ncol(X)=%d; ncol(tR)=%d <?< pp=%d (=?= n)\n",
                                n, p, full_p, ncol(tR), pp))
        if (ncol(tR) < pp) { # Add extra rows & cols if needed
          tR <- cbind(tR, Matrix::Matrix(0, nrow = pp, ncol = pp - ncol(tR), sparse = TRUE))
          if(verbose)
            cat(sprintf("    new tR: ncol(tR)=%d =!?= %d = pp = nrow(tR)\n", ncol(tR), pp))
        }
        ## Pad diagonal with ones
        d <- c(pp,pp) ; tR[.row(d) > p  &  .row(d) == .col(d)] <- 1
        ## Null basis is last pp-p cols of Q in QR decomposition of tR
        nbasis <- Matrix::qr.Q(Matrix::qr(tR))[, (p+1L):pp, drop = FALSE]
        ## Determine estimability; nbasis is orthonormal :
        ## *remember* rows are in qrX$pivot order; we use ALL columns of X, not just p of them
        Xb <- X[, qrXf$pivot] %*% nbasis
        ## simple vector norm.  norm() breaks if there are NAs
        norm2 <- function(x) sqrt(sum(x*x))
        Xb.norm <- apply(Xb, 1L, norm2) # = ||N'k||
        X.norm  <- apply(X , 1L, norm2) # = ||k'k||
        ## Find indices of non-estimable cases;  estimable <==> "Xb is basically 0"
        nonest <- which(tol * X.norm <= Xb.norm)
        if(rankdeficient == "warnif" && length(nonest))# warn only if there's a case
          warning(gettextf('%s; attr(*, "non-estim") has doubtful cases', msg1), domain=NA)
        ## newdata[i, ] is doubtful for indices in attr(*, "non-estim")
      }
    } # otherwise everything is estimable
    
    beta <- object$coefficients
    if(type != "terms") { # type == "terms" re-computes {predictor, ip}
      predictor <- drop(X[, piv, drop = FALSE] %*% beta[piv])
      if (!is.null(offset))
        predictor <- predictor + offset
      if(startsWith(rankdeficient, "NA") && length(nonest)) {
        predictor[nonest] <- NA
        if(rankdeficient == "NAwarn")
          warning(gettextf("%s: NAs produced for non-estimable cases", msg1), domain=NA)
      }
      else if(rankdeficient == "non-estim" || (hasNonest && length(nonest)))
        attr(predictor, "non-estim") <- nonest
    }
    interval <- match.arg(interval)
    if (interval == "prediction") {
      if (missing(newdata))
        warning("predictions on current data refer to _future_ responses\n")
      if (missing(newdata) && missing(weights)) {
        w <-  weights.default(object)
        if (!is.null(w)) {
          weights <- w
          warning("assuming prediction variance inversely proportional to weights used for fitting\n")
        }
      }
      if (!missing(newdata) && missing(weights) && !is.null(object$weights) && missing(pred.var))
        warning("Assuming constant prediction variance even though model fit is weighted\n")
      if (inherits(weights, "formula")){
        if (length(weights) != 2L)
          stop("'weights' as formula should be one-sided")
        d <- if(noData) model.frame(object) else newdata
        weights <- eval(weights[[2L]], d, environment(weights))
      }
    }
    
    if(se.fit || interval != "none") {
      ## w is needed for interval = "confidence"
      w <- object$weights
      res.var <-
        if (is.null(scale)) {
          r <- object$residuals
          rss <- sum(if(is.null(w)) r^2 else r^2 * w)
          df <- object$df.residual
          rss/df
        } else scale^2
      if(type != "terms") {
        if(p > 0) {
          XRinv <-
            if(missing(newdata) && is.null(w))
              Matrix::qr.Q(qrX$qr)[, p1, drop = FALSE]
          else
            X[, piv] %*% Matrix::solve(Matrix::qr.R(qrX$qr, backPermute = TRUE)[p1, p1])
          #	NB:
          #	 qr.Q(qrX)[, p1, drop = FALSE] / sqrt(w)
          #	looks faster than the above, but it's slower, and doesn't handle zero
          #	weights properly
          #
          ip <- drop(XRinv^2 %*% rep(res.var, p))
        } else ip <- rep(0, n)
      }
    }
    
    if (type == "terms") { ## type == "terms" ------ re-compute {predictor, ip} from scratch
      ## FIXME: if(hasNonest)  we are *not* yet obeying `rankdeficient`
      ## {offset is *not* considered (ok ?)}
      if(!mmDone) {
        mm <- Matrix::sparse.model.matrix(object)
        mmDone <- TRUE
      }
      aa <- attr(mm, "assign")
      ll <- attr(tt, "term.labels")
      hasintercept <- attr(tt, "intercept") > 0L
      if (hasintercept) ll <- c("(Intercept)", ll)
      aaa <- factor(aa, labels = ll)
      asgn <- split(order(aa), aaa)
      if (hasintercept) {
        asgn$"(Intercept)" <- NULL
        avx <- Matrix::colMeans(mm)
        termsconst <- sum(avx[piv] * beta[piv])
      }
      nterms <- length(asgn)
      if(nterms > 0) {
        predictor <- as(Matrix::Matrix(NA, ncol = nterms, nrow = NROW(X), sparse = TRUE),
                        "dMatrix")
        dimnames(predictor) <- list(rownames(X), names(asgn))
        
        if (se.fit || interval != "none") {
          ip <- predictor
          Rinv <- Matrix::solve(Matrix::qr.R(object$qr$qr, backPermute = TRUE)[p1, p1])
        }
        if(hasintercept)
          X <- sweep(X, 2L, avx, check.margin=FALSE)
        unpiv <- rep.int(0L, NCOL(X))
        unpiv[piv] <- p1
        ## Predicted values will be set to 0 for any term that
        ## corresponds to columns of the X-matrix that are
        ## completely aliased with earlier columns.
        for (i in seq.int(1L, nterms, length.out = nterms)) {
          iipiv <- asgn[[i]]      # Columns of X, ith term
          ii <- unpiv[iipiv]      # Corresponding rows of Rinv
          iipiv[ii == 0L] <- 0L
          predictor[, i] <-
            if(any(iipiv > 0L)) X[, iipiv, drop = FALSE] %*% beta[iipiv]
          else 0
          if (se.fit || interval != "none")
            ip[, i] <-
            if(any(iipiv > 0L))
              as.matrix(X[, iipiv, drop = FALSE] %*%
                          Rinv[ii, , drop = FALSE])^2 %*% rep.int(res.var, p)
          else 0
        }
        if (!is.null(terms)) {
          predictor <- predictor[, terms, drop = FALSE]
          if (se.fit)
            ip <- ip[, terms, drop = FALSE]
        }
      } else {                        # no terms
        predictor <- ip <- matrix(0, n, 0L)
      }
      attr(predictor, 'constant') <- if (hasintercept) termsconst else 0
    } ## type == "terms"
    
    ### Now construct elements of the list that will be returned
    
    if(interval != "none") {
      tfrac <- qt((1 - level)/2, df)
      hwid <- tfrac * switch(interval,
                             confidence = sqrt(ip),
                             prediction = sqrt(ip+pred.var)
      )
      if(type != "terms") {
        predictor <- cbind(predictor, predictor + hwid %o% c(1, -1))
        colnames(predictor) <- c("fit", "lwr", "upr")
      } else {
        if (!is.null(terms)) hwid <- hwid[, terms, drop = FALSE]
        lwr <- predictor + hwid
        upr <- predictor - hwid
      }
    }
    if(se.fit || interval != "none") {
      se <- sqrt(ip)
      if(type == "terms" && !is.null(terms) && !se.fit)
        se <- se[, terms, drop = FALSE]
    }
    if(missing(newdata) && !is.null(na.act <- object$na.action)) {
      predictor <- napredict(na.act, predictor)
      if(se.fit) se <- napredict(na.act, se)
    }
    if(type == "terms" && interval != "none") {
      if(missing(newdata) && !is.null(na.act)) {
        lwr <- napredict(na.act, lwr)
        upr <- napredict(na.act, upr)
      }
      list(fit = predictor, se.fit = se, lwr = lwr, upr = upr,
           df = df, residual.scale = sqrt(res.var))
    } else if (se.fit)
      list(fit = predictor, se.fit = se,
           df = df, residual.scale = sqrt(res.var))
    else predictor
  } ## {predict.lm}
