# glm2 1.2.1

* Changed `glm2` to call the fit function without `singular.ok` when the R version is less than 3.5.0, allowing backwards compatibility to 3.2.0. (Thanks to Paul Bailey for the suggestion)

# glm2 1.2

* `glm2` and `glm.fit2` now accept `singular.ok` argument, for compatibility with `glm` in R 3.5.0.

# glm2 1.1.3

* Removed some examples for compatibility with changes to `glm` in R 3.5.0.

# glm2 1.1.2

* Pass `singular.ok = FALSE` to `lm.fit` so that singular models will stop with an error.

# glm2 1.1.1

* Extract the `qr` component correctly.