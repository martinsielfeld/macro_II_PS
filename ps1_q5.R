#############################################
##
## Q5 Shooting algorithm
## Robust RK4
##
## Author: Martin Sielfeld
## Created: 02/10/2026
## Last edition: 02/12/2026
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

par <- list(
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

# single depreciation
w_e <- par$alpha_e / par$alpha
w_s <- par$alpha_s / par$alpha
par$delta <- w_e * par$delta_e + w_s * par$delta_s

cat("Depreciation:\n")
cat(sprintf("delta_e (paper) = %.3f\n", par$delta_e))
cat(sprintf("delta_s (paper) = %.3f\n", par$delta_s))
cat(sprintf("delta   (PS)    = %.6f\n\n", par$delta))


################################################
## Steady state (detrended) with l endogenous ##
################################################

par$Gamma <- (par$alpha_e / par$alpha)^par$alpha_e *
  (par$alpha_s / par$alpha)^par$alpha_s

ktilde_from_l <- function(l, par) {
  A <- (par$delta + par$rho + par$g_z) /
    (par$alpha * par$Gamma * l^(1 - par$alpha))
  A^(1 / (par$alpha - 1))
}

yc_from_l <- function(l, par) {
  ktilde <- ktilde_from_l(l, par)
  ytilde <- par$Gamma * ktilde^par$alpha * l^(1 - par$alpha)
  ctilde <- ytilde - (par$delta + par$g_z) * ktilde
  list(ktilde = ktilde, ytilde = ytilde, ctilde = ctilde)
}

intra_resid <- function(l, par) {
  if (!is.finite(l) || l <= 0 || l >= 1) {
    return(NA_real_)
  }
  tmp <- yc_from_l(l, par)
  if (!is.finite(tmp$ctilde) || tmp$ctilde <= 0) {
    return(NA_real_)
  }
  lhs <- (1 - par$theta) / par$theta * (tmp$ctilde / (1 - l))
  rhs <- (1 - par$alpha) * (tmp$ytilde / l)
  lhs - rhs
}

find_bracket_l <- function(par, grid = seq(0.05, 0.95, by = 0.01)) {
  vals <- sapply(grid, intra_resid, par = par)
  ok <- is.finite(vals)
  grid <- grid[ok]
  vals <- vals[ok]
  sgn <- sign(vals)
  idx <- which(sgn[-1] * sgn[-length(sgn)] < 0)
  if (length(idx) == 0) {
    stop("No sign change for intratemporal residual.")
  }
  c(grid[idx[1]], grid[idx[1] + 1])
}

br_l <- find_bracket_l(par)
l_star <- uniroot(
  intra_resid,
  lower = br_l[1],
  upper = br_l[2],
  par = par,
  tol = 1e-12
)$root
tmp <- yc_from_l(l_star, par)

ktilde_star <- tmp$ktilde
ytilde_star <- tmp$ytilde
ctilde_star <- tmp$ctilde

cat("Steady state (detrended):\n")
cat(sprintf("l*      = %.10f\n", l_star))
cat(sprintf("ktilde* = %.10f\n", ktilde_star))
cat(sprintf("ytilde* = %.10f\n", ytilde_star))
cat(sprintf("ctilde* = %.10f\n\n", ctilde_star))

par$l <- l_star


#####################################################
## Dynamics (ktilde,ctilde), disinvestment allowed ##
#####################################################

rhs <- function(state, par) {
  k <- state[1]
  c <- state[2]
  l <- par$l
  if (!is.finite(k) || !is.finite(c) || k <= 0 || c <= 0) {
    return(c(NA_real_, NA_real_))
  }
  y <- par$Gamma * k^par$alpha * l^(1 - par$alpha)
  dk <- y - c - (par$delta + par$g_z) * k
  dc <- c * (par$alpha * (y / k) - par$delta - par$rho - par$g_z)
  c(dk, dc)
}

# Check SS consistency with ODEs:
cat("SS residuals (dk, dc):\n")
print(rhs(c(ktilde_star, ctilde_star), par))
cat("\n")

rk4_step <- function(state, dt, par) {
  k1 <- rhs(state, par)
  if (anyNA(k1)) {
    return(NULL)
  }
  k2 <- rhs(state + 0.5 * dt * k1, par)
  if (anyNA(k2)) {
    return(NULL)
  }
  k3 <- rhs(state + 0.5 * dt * k2, par)
  if (anyNA(k3)) {
    return(NULL)
  }
  k4 <- rhs(state + dt * k3, par)
  if (anyNA(k4)) {
    return(NULL)
  }
  s_next <- state + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
  if (any(!is.finite(s_next)) || s_next[1] <= 0 || s_next[2] <= 0) {
    return(NULL)
  }
  s_next
}

simulate_path_rk4 <- function(
  k0,
  c0,
  par,
  T_end = 200,
  dt = 0.001,
  k_min = 1e-12,
  c_min = 1e-12,
  k_max = 1e8,
  c_max = 1e8
) {
  n <- as.integer(T_end / dt)
  out <- data.frame(
    time = numeric(n + 1),
    ktilde = numeric(n + 1),
    ctilde = numeric(n + 1)
  )

  state <- c(k0, c0)
  out[1, ] <- c(0, state[1], state[2])

  status <- "success"
  last_j <- 1

  for (j in 1:n) {
    # bounds stop
    if (
      !is.finite(state[1]) ||
        !is.finite(state[2]) ||
        state[1] <= k_min ||
        state[2] <= c_min ||
        state[1] >= k_max ||
        state[2] >= c_max
    ) {
      status <- if (state[1] >= k_max || state[2] >= c_max) {
        "exploded"
      } else {
        "hit_lower"
      }
      last_j <- j
      break
    }

    nxt <- rk4_step(state, dt, par)
    if (is.null(nxt)) {
      status <- "failed"
      last_j <- j
      break
    }

    state <- nxt
    out[j + 1, ] <- c(j * dt, state[1], state[2])
    last_j <- j + 1
  }

  out <- out[1:last_j, ]
  attr(out, "status") <- status
  out
}


########################
## Shooting objective ##
########################

# Key fix: if simulation stops early, use sign(k_last - k*) to classify
terminal_error <- function(c0, k0, par, k_star, T_end, dt) {
  if (!is.finite(c0) || c0 <= 0) {
    return(-1e6)
  }

  out <- simulate_path_rk4(k0, c0, par, T_end = T_end, dt = dt)
  status <- attr(out, "status")
  t_last <- tail(out$time, 1)
  k_last <- tail(out$ktilde, 1)

  # base sign from where k ended up
  base <- (k_last - k_star) / k_star

  # if it ended early, amplify but DO NOT force negative
  if (status != "success") {
    ampl <- 1 + 50 * (T_end - t_last) / T_end
    if (!is.finite(base)) {
      base <- -1
    }
    return(ampl * base)
  }

  base
}

# Bracket search: expand multiplicatively around c*
find_c_bracket <- function(
  k0,
  par,
  k_star,
  c_star,
  T_end,
  dt,
  mult_grid = exp(seq(log(1e-6), log(1e6), length.out = 121))
) {
  grid <- c_star * mult_grid
  vals <- sapply(
    grid,
    terminal_error,
    k0 = k0,
    par = par,
    k_star = k_star,
    T_end = T_end,
    dt = dt
  )

  ok <- is.finite(vals)
  grid <- grid[ok]
  vals <- vals[ok]

  sgn <- sign(vals)
  idx <- which(sgn[-1] * sgn[-length(sgn)] < 0)

  cat(sprintf(
    "\nBracket scan: finite=%d/%d, f range=[%.3g, %.3g]\n",
    length(vals),
    length(mult_grid),
    min(vals),
    max(vals)
  ))

  if (length(idx) == 0) {
    imn <- which.min(vals)
    imx <- which.max(vals)
    cat(sprintf("  min f=%.4g at c0=%.6g\n", vals[imn], grid[imn]))
    cat(sprintf("  max f=%.4g at c0=%.6g\n", vals[imx], grid[imx]))
    stop(
      "No sign change found in bracket scan. This usually indicates Euler/detrending mismatch or too-short T_end."
    )
  }

  c(grid[idx[1]], grid[idx[1] + 1])
}


##########################################
## Solve for c0 and simulate transition ##
##########################################

T_end <- 40
dt <- 0.001
k0 <- 0.99 * ktilde_star

br_c <- find_c_bracket(k0, par, ktilde_star, ctilde_star, T_end, dt)
cat(sprintf("Bracket for c0: [%.10f, %.10f]\n\n", br_c[1], br_c[2]))

c0_star <- uniroot(
  f = function(c0) terminal_error(c0, k0, par, ktilde_star, T_end, dt),
  lower = br_c[1],
  upper = br_c[2],
  tol = 1e-10,
  maxiter = 200
)$root

cat(sprintf("Shooting solution c0 = %.10f\n\n", c0_star))

path <- simulate_path_rk4(k0, c0_star, par, T_end = T_end, dt = dt)
cat("RK4 status:", attr(path, "status"), " t_last=", max(path$time), "\n")


######################################################
## Levels + dots (analytic) + ke/ks levels and dots ##
######################################################

# levels
path$ytilde <- par$Gamma * path$ktilde^par$alpha * par$l^(1 - par$alpha)
path$itilde <- path$ytilde - path$ctilde

# numeric diagnostic (already had)
path$kdot_num <- c(NA_real_, diff(path$ktilde) / dt)

# analytic dots from model equations
path$kdot <- path$ytilde - path$ctilde - (par$delta + par$g_z) * path$ktilde
path$cdot <- path$ctilde *
  (par$alpha * (path$ytilde / path$ktilde) - par$delta - par$rho - par$g_z)

# ydot via chain rule (l fixed): y = Gamma k^alpha l^(1-alpha)
path$ydot <- (par$alpha * path$ytilde / path$ktilde) * path$kdot

# idot = ydot - cdot
path$idot <- path$ydot - path$cdot

# recovered capital types (PS one-capital simplification)
w_e <- par$alpha_e / par$alpha
w_s <- par$alpha_s / par$alpha

path$ktilde_e <- w_e * path$ktilde
path$ktilde_s <- w_s * path$ktilde

path$kdot_e <- w_e * path$kdot
path$kdot_s <- w_s * path$kdot


###################################################
## Steady-state itilde*, ke*, ks* (and dots = 0) ##
###################################################

ktilde_e_star <- w_e * ktilde_star
ktilde_s_star <- w_s * ktilde_star
itilde_star <- ytilde_star - ctilde_star

cat(sprintf("\nSteady-state implied objects:\n"))
cat(sprintf("itilde*    = %.10f\n", itilde_star))
cat(sprintf("ktilde_e*  = %.10f\n", ktilde_e_star))
cat(sprintf("ktilde_s*  = %.10f\n", ktilde_s_star))
cat(sprintf(
  "check sum  = %.10f (should equal ktilde*)\n",
  ktilde_e_star + ktilde_s_star
))


####################################
## Terminal diagnostics (levels)  ##
####################################

kT <- tail(path$ktilde, 1)
cT <- tail(path$ctilde, 1)
yT <- tail(path$ytilde, 1)
iT <- tail(path$itilde, 1)

cat(sprintf(
  "Terminal: k(T)=%.10f (rel err %.3e)\n",
  kT,
  (kT - ktilde_star) / ktilde_star
))
cat(sprintf(
  "          c(T)=%.10f (rel err %.3e)\n",
  cT,
  (cT - ctilde_star) / ctilde_star
))
cat(sprintf(
  "          y(T)=%.10f (rel err %.3e)\n",
  yT,
  (yT - ytilde_star) / ytilde_star
))
cat(sprintf(
  "          i(T)=%.10f (rel err %.3e)\n\n",
  iT,
  (iT - itilde_star) / itilde_star
))


#############
## Figures ##
#############

## Collapse data:
path_long = melt.data.table(
  data.table(path),
  id.vars = "time",
  variable.name = "variable"
)

## Include dots in the steady-state reference lines as 0:
ss = data.table(
  variable = c(
    "ktilde",
    "ktilde_e",
    "ktilde_s",
    "ctilde",
    "ytilde",
    "itilde",
    "kdot_num",
    "kdot",
    "cdot",
    "ydot",
    "idot",
    "kdot_e",
    "kdot_s"
  ),
  value = c(
    ktilde_star,
    ktilde_e_star,
    ktilde_s_star,
    ctilde_star,
    ytilde_star,
    itilde_star,
    0,
    0,
    0,
    0,
    0,
    0,
    0
  )
)

## Factors:
path_long[,
  lab := factor(
    variable,
    levels = c(
      'ktilde',
      'ktilde_e',
      'ktilde_s',
      'ctilde',
      'ytilde',
      'itilde',
      'kdot_num',
      'kdot',
      'kdot_e',
      'kdot_s',
      'cdot',
      'ydot',
      'idot'
    ),
    labels = c(
      'Detrended Capital',
      'Detrended Equipment Capital',
      'Detrended Structure Capital',
      'Detrended Consumption',
      'Detrended Production',
      'Detrended Investment',
      'Diagnostic K dot',
      'K dot (model)',
      'Equipment K dot (model)',
      'Structure K dot (model)',
      'C dot (model)',
      'Y dot (model)',
      'I dot (model)'
    )
  )
]
ss[,
  lab := factor(
    variable,
    levels = c(
      'ktilde',
      'ktilde_e',
      'ktilde_s',
      'ctilde',
      'ytilde',
      'itilde',
      'kdot_num',
      'kdot',
      'kdot_e',
      'kdot_s',
      'cdot',
      'ydot',
      'idot'
    ),
    labels = c(
      'Detrended Capital',
      'Detrended Equipment Capital',
      'Detrended Structure Capital',
      'Detrended Consumption',
      'Detrended Production',
      'Detrended Investment',
      'Diagnostic K dot',
      'K dot (model)',
      'Equipment K dot (model)',
      'Structure K dot (model)',
      'C dot (model)',
      'Y dot (model)',
      'I dot (model)'
    )
  )
]


## Figures:
plot01 =
  ggplot() +
  geom_hline(
    data = ss[!variable %like% 'dot'],
    mapping = aes(yintercept = value, color = "Detrended Steady State")
  ) +
  geom_line(
    data = path_long[!variable %like% 'dot'],
    mapping = aes(x = time, y = value)
  ) +
  facet_wrap(~lab, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "bottom", plot.title.position = "plot") +
  labs(
    title = NULL,
    x = "Time",
    y = NULL,
    color = NULL
  )

plot01

ggsave(
  plot = plot01,
  filename = 'Transition Path in Detrended Levels.png',
  height = 4,
  width = 8,
  dpi = 300
)

plot02 =
  ggplot() +
  geom_hline(
    data = ss[variable %like% 'dot' & !variable %like% 'num'],
    mapping = aes(yintercept = value, color = "Detrended Steady State")
  ) +
  geom_line(
    data = path_long[variable %like% 'dot' & !variable %like% 'num'],
    mapping = aes(x = time, y = value)
  ) +
  facet_wrap(~lab, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "bottom", plot.title.position = "plot") +
  labs(
    title = NULL,
    x = "Time",
    y = NULL,
    color = NULL
  )

plot02

ggsave(
  plot = plot02,
  filename = 'Transition Dynamics in Detrended Time Derivatives.png',
  height = 4,
  width = 8,
  dpi = 300
)

## Checks / objects you might want to print:
ss
l_star
