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


## Model primitives / parameters
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

cat("Depreciation:\n")
cat(sprintf("delta_e (paper) = %.3f\n", par$delta_e))
cat(sprintf("delta_s (paper) = %.3f\n", par$delta_s))
cat(sprintf("delta   (PS)    = %.6f\n", par$delta))
