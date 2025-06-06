
<!-- README.md is generated from README.Rmd. Please edit that file -->

# glm2

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/glm2)](https://cran.r-project.org/package=glm2)

`glm2` fits generalized linear models using the same model specification
as `glm` in the `stats` package, but with a modified default fitting
method. The method provides greater stability for models that may fail
to converge using `glm`.

## Example

An example of periodic non-convergence using `glm`:

``` r
require(glm2, quietly = TRUE)
data(crabs)

fit1 <- glm(Satellites ~ I(Width - min(Width)) + factor(Dark) + factor(GoodSpine),
            family = poisson(link = "identity"), start = rep(1, 4),
            data = crabs, subset = Rep1, trace = TRUE)
#> Deviance = 679.309 Iterations - 1
#> Deviance = 704.3898 Iterations - 2
#> Deviance = 671.6167 Iterations - 3
#> Deviance = 687.4029 Iterations - 4
#> Deviance = 673.1342 Iterations - 5
#> Deviance = 690.9305 Iterations - 6
#> Deviance = 673.2107 Iterations - 7
#> Deviance = 691.0754 Iterations - 8
#> Deviance = 673.1994 Iterations - 9
#> Deviance = 691.049 Iterations - 10
#> Deviance = 673.1997 Iterations - 11
#> Deviance = 691.05 Iterations - 12
#> Deviance = 673.1998 Iterations - 13
#> Deviance = 691.0501 Iterations - 14
#> Deviance = 673.1998 Iterations - 15
#> Deviance = 691.0501 Iterations - 16
#> Deviance = 673.1998 Iterations - 17
#> Deviance = 691.0501 Iterations - 18
#> Deviance = 673.1998 Iterations - 19
#> Deviance = 691.0501 Iterations - 20
#> Deviance = 673.1998 Iterations - 21
#> Deviance = 691.0501 Iterations - 22
#> Deviance = 673.1998 Iterations - 23
#> Deviance = 691.0501 Iterations - 24
#> Deviance = 673.1998 Iterations - 25
#> Warning: glm.fit: algorithm did not converge
```

The same model converges successfully with `glm2`:

``` r
fit2 <- glm2(Satellites ~ I(Width - min(Width)) + factor(Dark) + factor(GoodSpine),
             family = poisson(link = "identity"), start = rep(1, 4),
             data = crabs, subset = Rep1, trace = TRUE)
#> Deviance = 679.309 Iterations - 1 
#> Deviance = 704.3898 Iterations - 2
#> Warning: step size truncated due to increasing deviance
#> Step halved: new deviance = 665.461 
#> Deviance = 664.8571 Iterations - 3 
#> Deviance = 672.5481 Iterations - 4
#> Warning: step size truncated due to increasing deviance
#> Step halved: new deviance = 658.3299 
#> Deviance = 658.6043 Iterations - 5
#> Warning: step size truncated due to increasing deviance
#> Step halved: new deviance = 656.3354 
#> Deviance = 656.339 Iterations - 6
#> Warning: step size truncated due to increasing deviance
#> Step halved: new deviance = 656.3126 
#> Deviance = 656.3116 Iterations - 7 
#> Deviance = 656.3116 Iterations - 8

noquote(c("deviances: ", fit1$deviance, fit2$deviance))
#> [1] deviances:       673.199771143023 656.311586320185
noquote(c("converged: ", fit1$converged, fit2$converged))
#> [1] converged:  FALSE       TRUE
```

## Installation

Get the released version from CRAN:

``` r
install.packages("glm2")
```

Or the development version from github:

``` r
# install.packages("devtools")
devtools::install_github("mdonoghoe/glm2")
```

## References

- Marschner, I. C. (2011). glm2: Fitting generalized linear models with
  convergence problems. *The R Journal* **3**(2): 12–15.
