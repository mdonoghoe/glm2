#  This function is based on File src/library/stats/R/glm.R
#  Part of the R package, https://www.R-project.org
#
#  Modified by Ian C. Marschner, 2011-2017
#  Modified by Mark W. Donoghoe:
#    25/05/2018 - accept singular.ok argument


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

utils::globalVariables("n", add = TRUE)

### Written by Simon Davies, Dec 1995
### glm.fit modified by Thomas Lumley, Apr 1997, and then others..

## Modified by Thomas Lumley 26 Apr 97
## Added boundary checks and step halving
## Modified detection of fitted 0/1 in binomial
## Updated by KH as suggested by BDR on 1998/06/16



#' Generalized Linear Models Fitting Method
#' 
#' \code{glm.fit2} is identical to \code{glm.fit} in the \bold{stats} package,
#' except for a modification to the computational method that provides improved
#' convergence properties. It is the default fitting method for \code{glm2} and
#' can also be used as an alternative fitting method for \code{glm}, instead of
#' the default method \code{glm.fit}.
#' 
#' \code{glm.fit2} is a modified version of \code{glm.fit} in the \bold{stats}
#' package. The computational method in \code{glm.fit2} uses a stricter form of
#' step-halving to deal with numerical instability in the iteratively
#' reweighted least squares algorithm. Whereas \code{glm.fit} uses step-halving
#' to correct divergence and parameter space violations, \code{glm.fit2}
#' additionally uses step-halving to force the model deviance to decrease at
#' each iteration, which improves the convergence properties. Like
#' \code{glm.fit}, \code{glm.fit2} usually would not be called directly.
#' Instead, it is called by \code{glm2} as the default fitting method. Like
#' \code{glm.fit}, it is possible to call \code{glm.fit2} directly if the
#' response vector and design matrix have already been calculated, in which
#' case it may be more efficient than calling \code{glm2}. It can also be
#' passed to \code{glm} in the \bold{stats} package, using the \code{method}
#' argument, providing an alternative to the default fitting method
#' \code{glm.fit}. See Marschner (2011) for a discussion of the need for a
#' modified fitting method.
#' 
#' @param x as for \code{\link{glm.fit}}
#' @param y as for \code{\link{glm.fit}}
#' @param weights as for \code{\link{glm.fit}}
#' @param start as for \code{\link{glm.fit}}
#' @param etastart as for \code{\link{glm.fit}}
#' @param mustart as for \code{\link{glm.fit}}
#' @param offset as for \code{\link{glm.fit}}
#' @param family as for \code{\link{glm.fit}}
#' @param control as for \code{\link{glm.fit}}
#' @param intercept as for \code{\link{glm.fit}}
#' @param singular.ok as for \code{\link{glm.fit}}
#' @return The value returned by \code{glm.fit2} has exactly the same
#' documentation as the value returned by \code{glm.fit}.
#' @author \code{glm.fit2} uses the code from \code{glm.fit}, whose authors are
#' listed in the help documentation for the \bold{stats} package. Modifications
#' to this code were made by Ian Marschner.
#' @seealso \code{\link{glm.fit}}
#' @references Marschner, I.C. (2011) glm2: Fitting generalized linear models
#' with convergence problems. The R Journal, Vol. 3/2, pp.12-15.
#' @examples
#' 
#' library(glm2)
#' 
#' #--- logistic regression null model ---------------#
#' # (intercept estimate = log(0.75/0.25) = 1.098612)
#' 
#' y <- c(1,1,1,0)
#' x <- rep(1,4)
#' 
#' #--- divergence of glm.fit to infinite estimate ---#
#' fit1 <- glm.fit(x,y, family=binomial(link="logit"),start=-1.81)
#' fit2 <- glm.fit2(x,y, family=binomial(link="logit"),start=-1.81)
#' print.noquote(c(fit1$coef,fit2$coef))
#' 
#' 
#' @export glm.fit2
glm.fit2 <- 
function (x, y, weights = rep(1, nobs), start = NULL, etastart = NULL, 
    mustart = NULL, offset = rep(0, nobs), family = gaussian(), 
    control = list(), intercept = TRUE, singular.ok = TRUE) 
{
    control <- do.call("glm.control", control)
    x <- as.matrix(x)
    xnames <- dimnames(x)[[2L]]
    ynames <- if (is.matrix(y)) rownames(y) else names(y)
    conv <- FALSE
    nobs <- NROW(y)
    nvars <- ncol(x)
    EMPTY <- nvars == 0
    ## define weights and offset if needed
    if (is.null(weights)) 
        weights <- rep.int(1, nobs)
    if (is.null(offset)) 
        offset <- rep.int(0, nobs)
    
    ## get family functions:
    variance <- family$variance
    linkinv <- family$linkinv
    if (!is.function(variance) || !is.function(linkinv)) 
        stop("'family' argument seems not to be a valid family object", 
            call. = FALSE)
    dev.resids <- family$dev.resids
    aic <- family$aic
    mu.eta <- family$mu.eta
    unless.null <- function(x, if.null) if (is.null(x)) if.null else x
    valideta <- unless.null(family$valideta, function(eta) TRUE)
    validmu <- unless.null(family$validmu, function(mu) TRUE)
    if (is.null(mustart)) {
        ## calculates mustart and may change y and weights and set n (!)
        eval(family$initialize)
    }
    else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
    }
    if (EMPTY) {
        eta <- rep.int(0, nobs) + offset
        if (!valideta(eta)) 
            stop("invalid linear predictor values in empty model", 
                call. = FALSE)
        mu <- linkinv(eta)
        ## calculate initial deviance and coefficient
        if (!validmu(mu)) 
            stop("invalid fitted means in empty model", call. = FALSE)
        dev <- sum(dev.resids(y, mu, weights))
        w <- ((weights * mu.eta(eta)^2)/variance(mu))^0.5
        residuals <- (y - mu)/mu.eta(eta)
        good <- rep(TRUE, length(residuals))
        boundary <- conv <- TRUE
        coef <- numeric()
        iter <- 0L
    }
    else {
        coefold <- NULL
        eta <- if (!is.null(etastart)) 
            etastart
        else if (!is.null(start)) 
            if (length(start) != nvars) 
                stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s", 
                  nvars, paste(deparse(xnames), collapse = ", ")), 
                  domain = NA)
            else {
                coefold <- start
                offset + as.vector(if (NCOL(x) == 1L) x*start else x%*%start)
            }
        else family$linkfun(mustart)
        mu <- linkinv(eta)
        if (!(validmu(mu) && valideta(eta))) 
            stop("cannot find valid starting values: please specify some", 
                call. = FALSE)
        ## calculate initial deviance and coefficient
        devold <- sum(dev.resids(y, mu, weights))
        boundary <- conv <- FALSE
        
        ##------------- THE Iteratively Reweighting L.S. iteration -----------
        for (iter in 1L:control$maxit) {
            good <- weights > 0
            varmu <- variance(mu)[good]
            if (any(is.na(varmu))) 
                stop("NAs in V(mu)")
            if (any(varmu == 0)) 
                stop("0s in V(mu)")
            mu.eta.val <- mu.eta(eta)
            if (any(is.na(mu.eta.val[good]))) 
                stop("NAs in d(mu)/d(eta)")
            ## drop observations for which w will be zero
            good <- (weights > 0) & (mu.eta.val != 0)
            
            if (all(!good)) {
                conv <- FALSE
                warning("no observations informative at iteration ", 
                  iter)
                break
            }
            z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
            w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
            ngoodobs <- as.integer(nobs - sum(!good))
            # Original glm.fit calls Fortran code
            fit <- lm.fit(x=x[good, , drop = FALSE]*w, y=z*w, singular.ok=TRUE, 
                          tol=min(1e-07, control$epsilon/1000))
            fit$coefficients[is.na(fit$coefficients)] <- 0
            if (any(!is.finite(fit$coefficients))) {
                conv <- FALSE
                warning(gettextf("non-finite coefficients at iteration %d", 
                  iter), domain = NA)
                break
            }
            ## stop if not enough parameters
            if (nobs < fit$rank) 
                stop(gettextf("X matrix has rank %d, but only %d observations", 
                  fit$rank, nobs), domain = NA)
            ## calculate updated values of eta and mu with the new coef:
            start <- fit$coefficients
            eta <- drop(x %*% start)
            mu <- linkinv(eta <- eta + offset)
            dev <- sum(dev.resids(y, mu, weights))
            if (control$trace) 
                cat("Deviance =", dev, "Iterations -", iter, 
                  "\n")
            ## check for divergence
            boundary <- FALSE
            if (!is.finite(dev)) {
                if (is.null(coefold)) 
                  stop("no valid set of coefficients has been found: please supply starting values", 
                    call. = FALSE)
                warning("step size truncated due to divergence", 
                  call. = FALSE)
                ii <- 1
                while (!is.finite(dev)) {
                  if (ii > control$maxit) 
                    stop("inner loop 1; cannot correct step size", 
                      call. = FALSE)
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  eta <- drop(x %*% start)
                  mu <- linkinv(eta <- eta + offset)
                  dev <- sum(dev.resids(y, mu, weights))
                }
                boundary <- TRUE
                if (control$trace) 
                  cat("Step halved: new deviance =", dev, "\n")
            }
            ## check for fitted values outside domain.
            if (!(valideta(eta) && validmu(mu))) {
                if (is.null(coefold)) 
                  stop("no valid set of coefficients has been found: please supply starting values", 
                    call. = FALSE)
                warning("step size truncated: out of bounds", 
                  call. = FALSE)
                ii <- 1
                while (!(valideta(eta) && validmu(mu))) {
                  if (ii > control$maxit) 
                    stop("inner loop 2; cannot correct step size", 
                      call. = FALSE)
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  eta <- drop(x %*% start)
                  mu <- linkinv(eta <- eta + offset)
                }
                boundary <- TRUE
                dev <- sum(dev.resids(y, mu, weights))
                if (control$trace) 
                  cat("Step halved: new deviance =", dev, "\n")
            }
            ## check for increasing deviance
            if (((dev - devold)/(0.1 + abs(dev)) >= control$epsilon)&(iter>1)) {
                if (is.null(coefold)) 
                  stop("no valid set of coefficients has been found: please supply starting values", 
                    call. = FALSE)
                warning("step size truncated due to increasing deviance", call. = FALSE)
                ii <- 1
                while ((dev - devold)/(0.1 + abs(dev)) > -control$epsilon) {
                  if (ii > control$maxit) break
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  eta <- drop(x %*% start)
                  mu <- linkinv(eta <- eta + offset)
                  dev <- sum(dev.resids(y, mu, weights))
                }
                if (ii > control$maxit) warning("inner loop 3; cannot correct step size")
                else if (control$trace) cat("Step halved: new deviance =", dev, "\n") 
            }
            ## check for convergence
            if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) {
                conv <- TRUE
                coef <- start
                break
            }
            else {
                devold <- dev
                coef <- coefold <- start
            }
        } ##-------------- end IRLS iteration -------------------------------
        
        if (!conv) 
            warning("glm.fit2: algorithm did not converge. Try increasing the maximum iterations", call. = FALSE)
        if (boundary) 
            warning("glm.fit2: algorithm stopped at boundary value", 
                call. = FALSE)
        eps <- 10 * .Machine$double.eps
        if (family$family == "binomial") {
            if (any(mu > 1 - eps) || any(mu < eps)) 
                warning("glm.fit2: fitted probabilities numerically 0 or 1 occurred", 
                  call. = FALSE)
        }
        if (family$family == "poisson") {
            if (any(mu < eps)) 
                warning("glm.fit2: fitted rates numerically 0 occurred", 
                  call. = FALSE)
        }
        ## If X matrix was not full rank then columns were pivoted,
        ## hence we need to re-label the names ...
        ## Original code changed as suggested by BDR---give NA rather
        ## than 0 for non-estimable parameters
        if (fit$rank < nvars) {
            if (!singular.ok) stop("singular fit encountered")
            coef[fit$qr$pivot][seq.int(fit$rank + 1, nvars)] <- NA
        }
        xxnames <- xnames[fit$qr$pivot]
        ## update by accurate calculation, including 0-weight cases.
        residuals <- (y - mu)/mu.eta(eta)
        fit$qr$qr <- as.matrix(fit$qr$qr)
        nr <- min(sum(good), nvars)
        if (nr < nvars) {
            Rmat <- diag(nvars)
            Rmat[1L:nr, 1L:nvars] <- fit$qr$qr[1L:nr, 1L:nvars]
        }
        else Rmat <- fit$qr$qr[1L:nvars, 1L:nvars]
        Rmat <- as.matrix(Rmat)
        Rmat[row(Rmat) > col(Rmat)] <- 0
        names(coef) <- xnames
        colnames(fit$qr$qr) <- xxnames
        dimnames(Rmat) <- list(xxnames, xxnames)
    }
    names(residuals) <- ynames
    names(mu) <- ynames
    names(eta) <- ynames
    # for compatibility with lm, which has a full-length weights vector
    wt <- rep.int(0, nobs)
    wt[good] <- w^2
    names(wt) <- ynames
    names(weights) <- ynames
    names(y) <- ynames
    if (!EMPTY) 
        names(fit$effects) <- 
            c(xxnames[seq_len(fit$rank)], rep.int("", sum(good) - fit$rank))
    ## calculate null deviance -- corrected in glm() if offset and intercept
    wtdmu <- 
        if (intercept) sum(weights * y)/sum(weights) else linkinv(offset)
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    ## calculate df
    n.ok <- nobs - sum(weights == 0)
    nulldf <- n.ok - as.integer(intercept)
    rank <- if (EMPTY) 0 else fit$rank
    resdf <- n.ok - rank
    ## calculate AIC
    aic.model <- 
        aic(y, n, mu, weights, dev) + 2 * rank
        ##     ^^ is only initialize()d for "binomial" [yuck!]
    list(coefficients = coef, residuals = residuals, fitted.values = mu, 
        effects = if (!EMPTY) fit$effects, R = if (!EMPTY) Rmat, 
        rank = rank, qr = if (!EMPTY) structure(fit$qr[c("qr", "rank", 
            "qraux", "pivot", "tol")], class = "qr"), family = family, 
        linear.predictors = eta, deviance = dev, aic = aic.model, 
        null.deviance = nulldev, iter = iter, weights = wt, prior.weights = weights, 
        df.residual = resdf, df.null = nulldf, y = y, converged = conv, 
        boundary = boundary)
}
