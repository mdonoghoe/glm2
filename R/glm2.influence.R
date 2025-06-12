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
influence.glm2Matrix <- function(model, do.coef = TRUE, ...) {
  
  stop("influence measures not implemented for glm.fit2.Matrix")
  
}

dfbetas.glm2Matrix <- function (model, infl = lm.influence(model, do.coef=TRUE), ...)
{
  
  stop("dfbetas not implemented for glm.fit2.Matrix")
  
}