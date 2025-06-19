#  This function is based on File src/library/stats/R/glm.R
#  Part of the R package, https://www.R-project.org
#
#  Modified by Ian C. Marschner, 2011-2017
#  Modified by Mark W. Donoghoe:
#    25/05/2018 - accept singular.ok argument
#    11/08/2018 - treat singular.ok differently depending on R version
#    31/01/2025 - version that employs Matrix


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

### This function fits a generalized linear model via
### iteratively reweighted least squares for any family.
### Written by Simon Davies, Dec 1995



#' Fitting Generalized Linear Models
#' 
#' Fits generalized linear models using the same model specification as
#' \code{glm} in the \bold{stats} package, but with a modified default fitting
#' method. The method provides greater stability for models that may fail to
#' converge using \code{glm}.
#' 
#' \code{glm2} is a modified version of \code{glm} in the \bold{stats} package.
#' It fits generalized linear models using the same model specification as
#' \code{glm}. It is identical to \code{glm} except for minor modifications to
#' change the default fitting method. The default method uses a stricter form
#' of step-halving to force the deviance to decrease at each iteration and is
#' implemented in \code{glm.fit2}. Like \code{glm}, user-supplied fitting
#' functions can be used with \code{glm2} by passing a function or a character
#' string naming a function to the \code{method} argument. See Marschner (2011)
#' for a discussion of the need for a modified fitting method.
#' 
#' @param formula as for \code{\link{glm}}
#' @param family as for \code{\link{glm}}
#' @param data as for \code{\link{glm}}
#' @param weights as for \code{\link{glm}}
#' @param subset as for \code{\link{glm}}
#' @param na.action as for \code{\link{glm}}
#' @param start as for \code{\link{glm}}
#' @param etastart as for \code{\link{glm}}
#' @param mustart as for \code{\link{glm}}
#' @param offset as for \code{\link{glm}}
#' @param control as for \code{\link{glm}} except by default \code{control} is
#' passed to \code{glm.fit2} instead of \code{glm.fit}
#' @param model as for \code{\link{glm}}
#' @param method the method used in fitting the model. The default method
#' \code{"glm.fit2"} uses iteratively reweighted least squares with modified
#' step-halving that forces the deviance to decrease at each iteration; see
#' help documentation for \code{glm.fit2}. As in \code{glm}, the alternative
#' method \code{"model.frame"} returns the model frame and does no fitting.
#' @param x as for \code{\link{glm}}
#' @param y as for \code{\link{glm}}
#' @param singular.ok as for \code{\link{glm}}. NB this is ignored (and
#' defaults to \code{TRUE}) for R versions < 3.5.0.
#' @param contrasts as for \code{\link{glm}}
#' @param \dots as for \code{\link{glm}}
#' @return The value returned by \code{glm2} has exactly the same documentation
#' as the value returned by \code{glm}, except for:
#' 
#' \item{method}{ the name of the fitter function used, which by default is
#' \code{"glm.fit2"}. }
#' @author \code{glm2} uses the code from \code{glm}, whose authors are listed
#' in the help documentation for the \bold{stats} package. Modifications to
#' this code were made by Ian Marschner.
#' @seealso \code{\link{glm}}
#' @references Marschner, I.C. (2011) glm2: Fitting generalized linear models
#' with convergence problems. The R Journal, Vol. 3/2, pp.12-15.
#' @examples
#' 
#' library(glm2)
#' data(crabs)
#' data(heart)
#' 
#' #==========================================================
#' # EXAMPLE 1: logistic regression null model
#' # (behaviour of glm and glm2 for different starting values)
#' #==========================================================
#' 
#' y <- c(1,1,1,0)
#' # intercept estimate = log(0.75/0.25) = 1.098612
#' 
#' #--- identical behaviour ---#
#' fit1 <- glm(y ~ 1, family=binomial(link="logit"),
#' 	control=glm.control(trace=TRUE))
#' fit2 <- glm2(y ~ 1, family=binomial(link="logit"),
#' 	control=glm.control(trace=TRUE))
#' print.noquote(c(fit1$coef,fit2$coef))
#' 
#' #--- convergence via different paths ---#
#' fit1 <- glm(y ~ 1, family=binomial(link="logit"),start=-1.75,
#' 	control=glm.control(trace=TRUE))
#' fit2 <- glm2(y ~ 1, family=binomial(link="logit"),start=-1.75,
#' 	control=glm.control(trace=TRUE))
#' print.noquote(c(fit1$coef,fit2$coef))
#' 
#' #--- divergence of glm to infinite estimate ---#
#' fit1 <- glm(y ~ 1, family=binomial(link="logit"),start=-1.81)
#' fit2 <- glm2(y ~ 1, family=binomial(link="logit"),start=-1.81)
#' print.noquote(c(fit1$coef,fit2$coef))
#' 
#' 
#' #=======================================================================
#' # EXAMPLE 2: identity link Poisson (successful boundary convergence
#' # using 4 identical approaches in glm and glm2 with the method argument) 
#' #=======================================================================
#' 
#' satellites <- crabs$Satellites
#' width.shifted <- crabs$Width - min(crabs$Width)
#' dark <- crabs$Dark
#' goodspine <- crabs$GoodSpine
#' 
#' fit1 <- glm(satellites ~ width.shifted + factor(dark) + factor(goodspine), 
#'  family = poisson(link="identity"), start = rep(1,4))
#' 
#' fit2 <- glm2(satellites ~ width.shifted + factor(dark) + factor(goodspine), 
#'  family = poisson(link="identity"), start = rep(1,4))
#' 
#' fit1.eq <- glm2(satellites ~ width.shifted + factor(dark) + factor(goodspine), 
#'  family = poisson(link="identity"), start = rep(1,4), method = "glm.fit")
#' 
#' fit2.eq <- glm(satellites ~ width.shifted + factor(dark) + factor(goodspine), 
#'  family = poisson(link="identity"), start = rep(1,4), method = "glm.fit2")
#' 
#' noquote(c("deviances: ",fit1$dev,fit2$dev,fit1.eq$dev,fit2.eq$dev))
#' noquote(c("converged: ",fit1$conv,fit2$conv,fit1.eq$conv,fit2.eq$conv))
#' noquote(c("boundary:  ",fit1$bound,fit2$bound,fit1.eq$bound,fit2.eq$bound))
#' 
#' #===================================================================
#' # EXAMPLE 3: identity link Poisson (periodic non-convergence in glm)
#' #===================================================================
#' 
#' R1 <- crabs$Rep1
#' satellites <- crabs$Satellites[R1]
#' width.shifted <- crabs$Width[R1] - min(crabs$Width)
#' dark <- crabs$Dark[R1]
#' goodspine <- crabs$GoodSpine[R1]
#' 
#' fit1 <- glm(satellites ~ width.shifted + factor(dark) + factor(goodspine), 
#'  family = poisson(link="identity"), start = rep(1,4), 
#'  control = glm.control(trace=TRUE))
#' 
#' fit2 <- glm2(satellites ~ width.shifted + factor(dark) + factor(goodspine), 
#'  family = poisson(link="identity"), start = rep(1,4), 
#'  control = glm.control(trace=TRUE))
#' 
#' noquote(c("deviances: ",fit1$dev,fit2$dev))
#' noquote(c("converged: ",fit1$conv,fit2$conv))
#' 
#' #===============================================================
#' # EXAMPLE 4: log link binomial (periodic non-convergence in glm)
#' #===============================================================
#' 
#' patients <- heart$Patients
#' deaths <- heart$Deaths
#' agegroup <- heart$AgeGroup
#' severity <-heart$Severity
#' delay <- heart$Delay
#' region <- heart$Region
#' start.p <- sum(deaths)/sum(patients)
#' 
#' fit1 <- glm(cbind(deaths,patients-deaths) ~ factor(agegroup) + factor(severity)
#'  + factor(delay) + factor(region), family = binomial(link="log"), 
#'  start = c(log(start.p), rep(0,8)), control = glm.control(trace=TRUE,maxit=100))
#' 
#' fit2 <- glm2(cbind(deaths,patients-deaths) ~ factor(agegroup) + factor(severity)
#'  + factor(delay) + factor(region), family = binomial(link="log"), 
#'  start = c(log(start.p), rep(0,8)), control = glm.control(trace=TRUE))
#' 
#' noquote(c("deviances: ",fit1$dev,fit2$dev))
#' noquote(c("converged: ",fit1$conv,fit2$conv))
#' 
#' #====================================================================
#' # EXAMPLE 5: identity link Poisson (aperiodic non-convergence in glm)
#' #====================================================================
#' 
#' R2 <- crabs$Rep2
#' satellites <- crabs$Satellites[R2]
#' width.shifted <- crabs$Width[R2] - min(crabs$Width)
#' dark <- crabs$Dark[R2]
#' goodspine <- crabs$GoodSpine[R2]
#' 
#' fit1 <- glm(satellites ~ width.shifted + factor(dark) + factor(goodspine), 
#'  family = poisson(link="identity"), start = rep(1,4), 
#'  control = glm.control(trace=TRUE))
#' 
#' fit2 <- glm2(satellites ~ width.shifted + factor(dark) + factor(goodspine), 
#'  family = poisson(link="identity"), start = rep(1,4), 
#'  control = glm.control(trace=TRUE))
#' 
#' noquote(c("deviances: ",fit1$dev,fit2$dev))
#' noquote(c("converged: ",fit1$conv,fit2$conv))
#' 
#' 
#' @importFrom stats .getXlevels is.empty.model model.extract model.matrix
#' model.offset model.response model.weights
#' @export glm2
glm2 <- 
function (formula, family = gaussian, data, weights, subset, 
    na.action, start = NULL, etastart, mustart, offset, control = list(...), 
    model = TRUE, method = "glm.fit2", x = FALSE, y = TRUE, 
    singular.ok = TRUE, contrasts = NULL,
    ...) 
{
    call <- match.call()
    ## family
    if (is.character(family)) 
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family)) 
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
    
    ## extract x, y, etc from the model formula and frame
    if (missing(data)) 
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "etastart", "mustart", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    if (identical(method, "model.frame")) return(mf)
    
    if (!is.character(method) && !is.function(method)) 
        stop("invalid 'method' argument")
    ## for back-compatibility in return result
    if (identical(method, "glm.fit2") || identical(method, "glm.fit2.Matrix")) 
        control <- do.call("glm.control", control)
    
    mt <- attr(mf, "terms") # allow model.frame to have updated it
    
    Y <- model.response(mf, "any") # e.g. factors are allowed
    ## avoid problems with 1D arrays, but keep names
    if (length(dim(Y)) == 1L) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm)) 
            names(Y) <- nm
    }
    ## null model support
    if (!is.empty.model(mt)) {
        if (identical(method, "glm.fit2.Matrix")) {
            if (!requireNamespace("Matrix", quietly = TRUE)) {
                stop("Package \"Matrix\" must be installed to use method = \"glm.fit2.Matrix\".",
                     call. = FALSE)
            }
            X <- Matrix::sparse.model.matrix(mt, mf, contrasts)
            if (control$trace) cat("Created sparse model matrix\n")
        } else X <- model.matrix(mt, mf, contrasts)
    } else X <- matrix(, NROW(Y), 0L)
    ## avoid any problems with 1D or nx1 arrays by as.vector.
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights)) 
        stop("'weights' must be a numeric vector")
    ## check weights and offset
    if (!is.null(weights) && any(weights < 0)) 
        stop("negative weights not allowed")
    
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) {
        if (length(offset) != NROW(Y)) 
            stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
                length(offset), NROW(Y)), domain = NA)
    }
    ## these allow starting values to be expressed in terms of other vars.
    mustart <- model.extract(mf, "mustart")
    etastart <- model.extract(mf, "etastart")
    
    ## Change singular.ok behaviour depending on R version
    R.maj <- as.numeric(R.version$major)
    R.min <- as.numeric(unlist(strsplit(R.version$minor, ".", TRUE))[1])
    if (R.maj > 3 | (R.maj == 3 & R.min >= 5)) {
      ## We want to set the name on this call and the one below for the
      ## sake of messages from the fitter function
      fit <- eval(call(if (is.function(method)) "method" else method, 
          x = X, y = Y, weights = weights, start = start, etastart = etastart, 
          mustart = mustart, offset = offset, family = family, 
          control = control, intercept = attr(mt, "intercept") > 
              0L, singular.ok = singular.ok))
    } else {
      if (!missing(singular.ok)) warning("singular.ok is ignored (defaults to TRUE) for R version < 3.5.0")
      fit <- eval(call(if (is.function(method)) "method" else method, 
                       x = X, y = Y, weights = weights, start = start, etastart = etastart, 
                       mustart = mustart, offset = offset, family = family, 
                       control = control, intercept = attr(mt, "intercept") > 
                         0L))
    }
    
    ## This calculated the null deviance from the intercept-only model
    ## if there is one, otherwise from the offset-only model.
    ## We need to recalculate by a proper fit if there is intercept and
    ## offset.
    ##
    ## The glm.fit calculation could be wrong if the link depends on the
    ## observations, so we allow the null deviance to be forced to be
    ## re-calculated by setting an offset (provided there is an intercept).
    ## Prior to 2.4.0 this was only done for non-zero offsets.
    if (length(offset) && attr(mt, "intercept") > 0L) {
        fit2 <- eval(call(if (is.function(method)) "method" else method, 
            x = X[, "(Intercept)", drop = FALSE], y = Y, 
            ## starting values potentially required (PR#16877):
            mustart = fit$fitted.values,
            weights = weights, offset = offset, family = family, 
            control = control, intercept = TRUE))
        ## That fit might not have converged ....
        if (!fit2$converged) 
            warning("fitting to calculate the null deviance did not converge -- increase maxit?")
        fit$null.deviance <- fit2$deviance
    }
    if (model) fit$model <- mf
    fit$na.action <- attr(mf, "na.action")
    if (x) fit$x <- X
    if (!y) fit$y <- NULL
    fit <- c(fit, list(call = call, formula = formula, terms = mt, 
        data = data, offset = offset, control = control, method = method, 
        contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, 
            mf)))
    class(fit) <- c(fit$class, c("glm", "lm"))
    fit
}
