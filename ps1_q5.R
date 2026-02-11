#############################################
##
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

## We fix for Q5 labor at the steady state level:
par$l = l_star


###############################################################
## Define the ODE system (ktilde, ctilde) with l fixed at l* ##
## and build the shooting objective that chooses ctilde(0)   ##
## so that the path converges to the steady state.           ##
###############################################################

## RHS for the 2D system with fixed l:
rhs = function(state, par) {
  k = state[1]
  c = state[2]
  l = par$l

  # Return NULL if infeasible
  if (!is.finite(k) || !is.finite(c) || k <= 0 || c <= 0 || l <= 0 || l >= 1) {
    return(NULL)
  }

  y = par$Gamma * k^par$alpha * l^(1 - par$alpha)

  dk = y - c - (par$delta + par$g_z) * k
  dc = c * (par$alpha * (y / k) - par$delta - par$rho - par$g_z)

  c(dk, dc)
}

## One RK4 step:
rk4_step = function(state, dt, par) {
  k1 = rhs(state, par)
  if (is.null(k1)) {
    return(NULL)
  }
  k2 = rhs(state + 0.5 * dt * k1, par)
  if (is.null(k2)) {
    return(NULL)
  }
  k3 = rhs(state + 0.5 * dt * k2, par)
  if (is.null(k3)) {
    return(NULL)
  }
  k4 = rhs(state + dt * k3, par)
  if (is.null(k4)) {
    return(NULL)
  }

  state_next = state + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

  # Feasibility check
  if (any(!is.finite(state_next)) || state_next[1] <= 0 || state_next[2] <= 0) {
    return(NULL)
  }

  state_next
}

## Simulate path with RK4 but STOP when it "crashes" or "explodes":
simulate_path_rk4 = function(
  k0,
  c0,
  par,
  T_end = 400,
  dt = 0.01,
  k_min = 1e-10,
  c_min = 1e-10,
  k_max = 1e6,
  c_max = 1e6
) {
  n = as.integer(T_end / dt)
  k = numeric(n + 1)
  c = numeric(n + 1)
  t = numeric(n + 1)

  k[1] = k0
  c[1] = c0
  t[1] = 0
  state = c(k0, c0)

  for (j in 1:n) {
    # If already out of bounds, stop and return what we have
    if (
      !is.finite(state[1]) ||
        !is.finite(state[2]) ||
        state[1] <= k_min ||
        state[2] <= c_min ||
        state[1] >= k_max ||
        state[2] >= c_max
    ) {
      k = k[1:j]
      c = c[1:j]
      t = t[1:j]
      return(data.frame(time = t, ktilde = k, ctilde = c))
    }

    state_next = rk4_step(state, dt, par)

    # If the step fails (infeasible), treat as crash: stop and return last valid
    if (is.null(state_next)) {
      k = k[1:j]
      c = c[1:j]
      t = t[1:j]
      return(data.frame(time = t, ktilde = k, ctilde = c))
    }

    state = state_next
    t[j + 1] = j * dt
    k[j + 1] = state[1]
    c[j + 1] = state[2]
  }

  data.frame(time = t, ktilde = k, ctilde = c)
}

## Terminal error: ALWAYS returns a number (no NA) using last simulated point
terminal_error = function(c0_guess, k0, par, k_star, T_end = 400, dt = 0.01) {
  if (!is.finite(c0_guess) || c0_guess <= 0) {
    return(NA_real_)
  }

  out = simulate_path_rk4(k0, c0_guess, par, T_end = T_end, dt = dt)

  kT = out$ktilde[nrow(out)]
  # if it crashed to tiny capital => very negative error
  if (!is.finite(kT) || kT <= 1e-10) {
    return(-1e6)
  }

  (kT - k_star) / k_star
}

## Bracket finder: now we KEEP crash cases (they give big negative errors)
find_c_bracket = function(k0, par, k_star, c_star, T_end = 400, dt = 0.01) {
  factors = c(
    1e-6,
    3e-6,
    1e-5,
    3e-5,
    1e-4,
    3e-4,
    1e-3,
    3e-3,
    1e-2,
    3e-2,
    1e-1,
    3e-1,
    1,
    3,
    10,
    30,
    100
  )

  grid = c_star * factors

  vals = sapply(
    grid,
    terminal_error,
    k0 = k0,
    par = par,
    k_star = k_star,
    T_end = T_end,
    dt = dt
  )

  ok = is.finite(vals)
  grid = grid[ok]
  vals = vals[ok]
  if (length(grid) < 2) {
    stop("Could not evaluate terminal errors. Try different dt/T_end.")
  }

  sgn = sign(vals)
  idx = which(sgn[-1] * sgn[-length(sgn)] < 0)

  if (length(idx) == 0) {
    cat("No sign change found.\n")
    cat(sprintf("Min err = %.6f, Max err = %.6f\n", min(vals), max(vals)))
    cat(sprintf("Grid range tried: [%.6g, %.6g]\n", min(grid), max(grid)))
    stop("Still no sign change. Increase factors range or change k0.")
  }

  c(grid[idx[1]], grid[idx[1] + 1])
}

ktilde0 = 0.99 * ktilde_star

br_c = find_c_bracket(
  k0 = ktilde0,
  par = par,
  k_star = ktilde_star,
  c_star = ctilde_star,
  T_end = 400,
  dt = 0.01
)

cat(sprintf("Bracket for ctilde(0): [%.10f, %.10f]\n", br_c[1], br_c[2]))

ctilde0_star = uniroot(
  f = terminal_error,
  lower = br_c[1],
  upper = br_c[2],
  k0 = ktilde0,
  par = par,
  k_star = ktilde_star,
  T_end = 400,
  dt = 0.01,
  tol = 1e-10,
  maxiter = 200
)$root

cat(sprintf("Shooting solution: ctilde(0) = %.10f\n", ctilde0_star))
