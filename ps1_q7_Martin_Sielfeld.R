#############################################
##
## Q7 Steady stat with increasing q(t)
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

## Original (component) deltas:
de <- 0.124
ds <- 0.056

## Representative depreciation (impose d = de = ds in the system)
w_e <- ae / alpha
w_s <- as / alpha
d <- w_e * de + w_s * ds

rho <- -log(0.9)
gz <- 0.0124
gq <- log(1.032)
gy <- gz + (ae / (1 - alpha)) * gq

cat(paste0(
  "Representative depreciation:\n",
  "  w_e = ",
  round(w_e, 6),
  "\n",
  "  w_s = ",
  round(w_s, 6),
  "\n",
  "  d   = ",
  round(d, 6),
  "\n\n"
))

## Constants (replace de, ds by d everywhere)
B <- (ae / (d + rho + gy - gq))^ae * (as / (d + rho + gy))^as
A <- B^(1 / (1 - alpha)) # so ytilde* = A * l*

## Functions:
ytilde_from_l <- function(l) A * l

ks_from_l <- function(l) (as / (d + rho + gy)) * ytilde_from_l(l)

ke_from_l <- function(l) (ae / (d + rho + gy - gq)) * ytilde_from_l(l)

ctilde_from_l <- function(l) {
  (theta / (1 - theta)) * (1 - alpha) * ytilde_from_l(l) * (1 - l) / l
}

## Solver:
resid_l <- function(l) {
  y <- ytilde_from_l(l)
  ks <- ks_from_l(l)
  ke <- ke_from_l(l)
  c <- ctilde_from_l(l)

  ## Investment flows in detrended units (replace de, ds by d)
  is <- (d + gy) * ks
  ie <- (d + gy + gq) * ke

  y - c - is - ie
}

l_star <- uniroot(resid_l, lower = 1e-6, upper = 1 - 1e-6)$root

## Get SS values:
y_star <- ytilde_from_l(l_star)
ks_star <- ks_from_l(l_star)
ke_star <- ke_from_l(l_star)
c_star <- ctilde_from_l(l_star)
is_star <- (d + gy) * ks_star
ie_star <- (d + gy + gq) * ke_star

###################
## Print results ##
###################

cat(paste0(
  "\n",
  "Detrended Ks ss: ",
  ks_star,
  "\n",
  "Detrended Ke ss: ",
  ke_star,
  "\n",
  "Detrended C  ss: ",
  c_star,
  "\n",
  "Detrended Y  ss: ",
  y_star,
  "\n",
  "Detrended L  ss: ",
  l_star,
  "\n",
  "Detrended Is ss: ",
  is_star,
  "\n",
  "Detrended Ie ss: ",
  ie_star,
  "\n"
))
