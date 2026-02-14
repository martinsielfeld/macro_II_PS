#############################################
##
## Q7 Steady stat with increasing q(t)
## Robust RK4
##
## Author: Martin Sielfeld
## Created: 02/14/2026
## Last edition: 02/14/2026
##
#############################################

rm(list = ls())
options(scipen = 999)

packages <- c("data.table", "ggplot2")
installed_packages <- packages %in% installed.packages()[, 1]
if (any(!installed_packages)) {
  install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))

###################################
## Model primitives / parameters ##
###################################

theta <- 0.4
ae <- 0.17
as <- 0.13
alpha <- ae + as
de <- 0.124
ds <- 0.056
rho <- -log(0.9)
gz <- 0.0124
gq <- log(1.032)
gy <- gz + (ae / (1 - alpha)) * gq

## Constants:
B <- (ae / (de + rho + gy - gq))^ae * (as / (ds + rho + gy))^as
A <- B^(1 / (1 - alpha)) # so ytilde* = A * l*

## Functions:
ytilde_from_l <- function(l) A * l
ks_from_l <- function(l) (as / (ds + rho + gy)) * ytilde_from_l(l)
ke_from_l <- function(l) (ae / (de + rho + gy - gq)) * ytilde_from_l(l)
ctilde_from_l <- function(l) {
  (theta / (1 - theta)) * (1 - alpha) * ytilde_from_l(l) * (1 - l) / l
}

## Solver:
resid_l <- function(l) {
  y <- ytilde_from_l(l)
  ks <- ks_from_l(l)
  ke <- ke_from_l(l)
  c <- ctilde_from_l(l)
  is <- (ds + gy) * ks
  ie <- (de + gy + gq) * ke
  y - c - is - ie
}

l_star <- uniroot(resid_l, lower = 1e-6, upper = 1 - 1e-6)$root

## Get SS values:
y_star <- ytilde_from_l(l_star)
ks_star <- ks_from_l(l_star)
ke_star <- ke_from_l(l_star)
c_star <- ctilde_from_l(l_star)
is_star <- (ds + gy) * ks_star
ie_star <- (de + gy + gq) * ke_star

## Print:
cat(paste0(
  '\n',
  'Detrended Ks ss:',
  ks_star,
  '\n',
  'Detrended Ke ss:',
  ke_star,
  '\n',
  'Detrended C ss:',
  c_star,
  '\n',
  'Detrended Y ss:',
  y_star,
  '\n',
  'Detrended L ss:',
  l_star,
  '\n',
  'Detrended Is ss:',
  is_star,
  '\n',
  'Detrended Ie ss:',
  ie_star,
  '\n'
))
