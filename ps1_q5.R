#############################################
## Shooting algorithm for Q5
##
## Author: Martin Sielfeld
## Created: 02/10/2026
## Last edition: 02/10/2026
##
## Source: Greenwood, J., Hercowitz, Z., & Krusell, P. (1997)
##
#############################################

## Settings:
rm(list = ls())
options(scipen = 999)

## Packages:
packages = c("deSolve", "rootSolve")
installed_packages = packages %in% installed.packages()[, 1]
if (any(installed_packages == F)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = T))


###################################
## Model primitives / parameters ##
###################################

par = list(
  theta = 0.4,
  beta = 0.95,
  alpha_e = 0.17,
  alpha_s = 0.13,
  alpha = 0.3,
  delta_s = 0.056,
  delta_e = 0.124,
  g_gross = 1.0124,
  g_z = 1.0124 - 1,
  rho = -log(0.9)
)

## Building a single depreciation rate delta for this PS (since paper has d_s and d_e)
w_e = par$alpha_e / par$alpha
w_s = par$alpha_s / par$alpha
par$delta = w_e * par$delta_e + w_s * par$delta_s

cat("Depreciation:\n") ## Check values
cat(sprintf("delta_e (paper) = %.3f\n", par$delta_e))
cat(sprintf("delta_s (paper) = %.3f\n", par$delta_s))
cat(sprintf("delta   (PS)    = %.6f\n", par$delta))


###########################################################################
## Now let's compute the steady state from Question 4 with endogenous l* ##
###########################################################################

## Fist lets build gamma for the one-capital repres.
par$Gamma = (par$alpha_e / par$alpha)^par$alpha_e *
  (par$alpha_s / par$alpha)^par$alpha_s

## ktilde as a function of l from the steady-state MPK condition:
ktilde_from_l = function(l, par) {
  A = (par$delta + par$rho + par$g_z) /
    (par$alpha * par$Gamma * l^(1 - par$alpha))
  A^(1 / (par$alpha - 1))
}

## Given l, compute (ytilde, ctilde) implied by steady-state kdot=0:
yc_from_l = function(l, par) {
  ktilde = ktilde_from_l(l, par)
  ytilde = par$Gamma * ktilde^par$alpha * l^(1 - par$alpha)
  ctilde = ytilde - (par$delta + par$g_z) * ktilde
  list(ktilde = ktilde, ytilde = ytilde, ctilde = ctilde)
}

## Intratemporal residual as a function of l alone:
intra_resid = function(l, par) {
  # keep away from boundaries
  if (!is.finite(l) || l <= 0 || l >= 1) {
    return(NA_real_)
  }

  tmp = yc_from_l(l, par)
  if (!is.finite(tmp$ctilde) || tmp$ctilde <= 0) {
    return(NA_real_)
  } # infeasible

  lhs = (1 - par$theta) / par$theta * (tmp$ctilde / (1 - l))
  rhs = (1 - par$alpha) * (tmp$ytilde / l)
  lhs - rhs
}

## Helper: find a valid bracket automatically
find_bracket = function(par, grid = seq(0.05, 0.95, by = 0.01)) {
  vals = sapply(grid, intra_resid, par = par)
  ok = is.finite(vals)
  grid = grid[ok]
  vals <- vals[ok]
  if (length(grid) < 2) {
    stop("No feasible points found for l in (0,1). Check parameters.")
  }

  sgn = sign(vals)
  idx = which(sgn[-1] * sgn[-length(sgn)] < 0)
  if (length(idx) == 0) {
    stop(
      "Could not find a sign change for intratemporal residual. Check grid or parameters."
    )
  }

  c(grid[idx[1]], grid[idx[1] + 1])
}

br = find_bracket(par)

l_star = uniroot(
  intra_resid,
  lower = br[1],
  upper = br[2],
  par = par,
  tol = 1e-12
)$root
tmp = yc_from_l(l_star, par)

ktilde_star = tmp$ktilde
ytilde_star = tmp$ytilde
ctilde_star = tmp$ctilde

## Q4 steady state (detrended, endogenous labor):
cat(sprintf("l*      = %.10f\n", l_star))
cat(sprintf("ktilde* = %.10f\n", ktilde_star))
cat(sprintf("ytilde* = %.10f\n", ytilde_star))
cat(sprintf("ctilde* = %.10f\n\n", ctilde_star))

## For Q5, fix labor at the steady state level:
par$l = l_star
