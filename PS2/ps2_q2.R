#############################################
##
## Problem Set II
## Question 2
##
## Author: Martin Sielfeld
## Created: 02/24/2026
## Last edition: 02/24/2026
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
  r = 0.05,
  lambda = 0.1,
  alpha = 1.5,
  c = 20,
  mu = 3.5,
  sigma = 0.5,
  g1 = 0.03,
  g2 = 0
)

############################################################################
## Q3: Compute the reservation wage and expected duration of unemployment ##
############################################################################

## We work with ln(w) ~ N(mu, sigma^2).
## We truncate the support at extreme quantiles only for numerical
## integration of expectations.

## First we create the PDF grid with small deltas:
p_low = 1e-8
p_high = 1 - 1e-8

w_min = qlnorm(p_low, meanlog = par$mu, sdlog = par$sigma)
w_max = qlnorm(p_high, meanlog = par$mu, sdlog = par$sigma)

n_grid = 8000
w_grid = seq(w_min, w_max, length.out = n_grid)
dw = w_grid[2] - w_grid[1]

pdf_w = dlnorm(w_grid, meanlog = par$mu, sdlog = par$sigma)
F_min = plnorm(w_min, meanlog = par$mu, sdlog = par$sigma)
F_max = plnorm(w_max, meanlog = par$mu, sdlog = par$sigma)
Z = F_max - F_min
pdf_w = pdf_w / Z

plot(pdf_w)

trapz = function(x, y) {
  ## Numerical integration rule:
  sum((y[-1] + y[-length(y)]) / 2 * diff(x))
}

## Given U, solve the employed HJB ODE for W(w):
## (r+lambda) W(w) = w + g w W_w(w) + lambda U
## <=> g w W_w(w) = (r+lambda)W(w) - w - lambda U
## This is a first-order linear ODE. We solve it backward from w_max using a boundary
## condition based on the linear particular solution for large w:
## W(w) ~ A w + B,  A = 1/(r+lambda-g),  B = lambda U/(r+lambda).

solve_W_given_U <- function(U, g, r, lambda, w_grid) {
  n <- length(w_grid)
  W <- numeric(n)

  if (abs(g) < 1e-12) {
    ## If g = 0, HJB is algebraic: rW = w + lambda(U - W)
    ## => (r+lambda)W = w + lambda U
    W <- (w_grid + lambda * U) / (r + lambda)
    return(W)
  }

  ## boundary at w_max using asymptotic linear solution
  A <- 1 / (r + lambda - g)
  B <- (lambda * U) / (r + lambda)
  W[n] <- A * w_grid[n] + B

  ## backward Euler / upwind solve from high w to low w
  ## ODE: W'(w) = ((r+lambda)W(w) - w - lambda U) / (g w)
  for (i in (n - 1):1) {
    w_next <- w_grid[i + 1]
    W_next <- W[i + 1]
    deriv_next <- ((r + lambda) * W_next - w_next - lambda * U) / (g * w_next)
    W[i] <- W_next - deriv_next * (w_grid[i + 1] - w_grid[i])
  }

  W
}

## Given arrays (w_grid, W), find reservation wage w_R from W(w_R) = U
find_wR <- function(U, w_grid, W) {
  diffWU <- W - U

  ## If even the smallest offer is worth accepting, cutoff is at the bottom
  if (diffWU[1] >= 0) {
    return(w_grid[1])
  }

  ## If even the largest offer is not worth accepting (should not happen typically),
  ## return top and let the outer solver deal with it.
  if (diffWU[length(diffWU)] <= 0) {
    return(w_grid[length(w_grid)])
  }

  ## Find first index where W crosses U from below
  j <- which(diffWU >= 0)[1]
  j0 <- j - 1

  ## Linear interpolation between (w_j0, W_j0) and (w_j, W_j)
  w0 <- w_grid[j0]
  w1 <- w_grid[j]
  W0 <- W[j0]
  W1 <- W[j]

  wR <- w0 + (U - W0) * (w1 - w0) / (W1 - W0)
  wR
}

## Compute W_w(w_R) by finite differences on the grid
Ww_at_wR <- function(wR, w_grid, W) {
  n <- length(w_grid)

  ## locate nearest index
  k <- findInterval(wR, w_grid, all.inside = TRUE)

  if (k <= 1) {
    ## forward diff
    return((W[2] - W[1]) / (w_grid[2] - w_grid[1]))
  }
  if (k >= n) {
    ## backward diff
    return((W[n] - W[n - 1]) / (w_grid[n] - w_grid[n - 1]))
  }

  ## central diff around k
  (W[k + 1] - W[k - 1]) / (w_grid[k + 1] - w_grid[k - 1])
}

## Unemployment HJB residual for a given U:
## rU - c - alpha * ∫ max{W(w)-U,0} dF(w) = 0
U_residual <- function(U, g, r, lambda, alpha, c, w_grid, pdf_w) {
  W <- solve_W_given_U(U, g, r, lambda, w_grid)
  wR <- find_wR(U, w_grid, W)

  gain <- pmax(W - U, 0) * pdf_w
  integ <- trapz(w_grid, gain)

  res <- r * U - c - alpha * integ

  list(res = res, W = W, wR = wR)
}

## Solve for (U, w_R) given g using a 1D root-find on U
solve_model <- function(g, par, w_grid, pdf_w) {
  r <- par$r
  lambda <- par$lambda
  alpha <- par$alpha
  c <- par$c

  ## Bracket for U. We choose a wide bracket and expand if needed.
  U_lo <- 0
  U_hi <- 1e6

  f_lo <- U_residual(U_lo, g, r, lambda, alpha, c, w_grid, pdf_w)$res
  f_hi <- U_residual(U_hi, g, r, lambda, alpha, c, w_grid, pdf_w)$res

  ## Expand upper bound until sign change (or hit a cap)
  cap <- 1e10
  while (f_lo * f_hi > 0 && U_hi < cap) {
    U_hi <- U_hi * 2
    f_hi <- U_residual(U_hi, g, r, lambda, alpha, c, w_grid, pdf_w)$res
  }

  if (f_lo * f_hi > 0) {
    stop("Could not bracket the root for U. Increase cap or check parameters.")
  }

  root <- uniroot(
    f = function(U) U_residual(U, g, r, lambda, alpha, c, w_grid, pdf_w)$res,
    lower = U_lo,
    upper = U_hi,
    tol = 1e-10
  )

  U_star <- root$root
  out <- U_residual(U_star, g, r, lambda, alpha, c, w_grid, pdf_w)
  W_star <- out$W
  wR_star <- out$wR

  ## Check reservation wage condition from the PS:
  ## rU = w_R + g w_R W_w(w_R)
  Ww_star <- Ww_at_wR(wR_star, w_grid, W_star)

  ## Acceptance probability (true lognormal CDF, not truncated)
  accept_prob <- 1 - plnorm(wR_star, meanlog = par$mu, sdlog = par$sigma)

  ## Expected unemployment duration:
  ## E[T_U] = 1 / (alpha * (1 - F(w_R)))
  ETU <- 1 / (par$alpha * accept_prob)

  list(
    g = g,
    U = U_star,
    wR = wR_star,
    Ww_wR = Ww_star,
    check_reservation = r * U_star - (wR_star + g * wR_star * Ww_star),
    accept_prob = accept_prob,
    ETU = ETU
  )
}

## Run for g = 0.03 and g = 0
sol_g1 <- solve_model(g = par$g1, par = par, w_grid = w_grid, pdf_w = pdf_w)
sol_g2 <- solve_model(g = par$g2, par = par, w_grid = w_grid, pdf_w = pdf_w)

cat(sprintf(
  "g = %.4f: U = %.6f, w_R = %.6f, E[T_U] = %.6f\n",
  sol_g1$g,
  sol_g1$U,
  sol_g1$wR,
  sol_g1$ETU
))
cat(sprintf(
  "reservation check (should be ~0): %.3e\n",
  sol_g1$check_reservation
))

cat(sprintf(
  "g = %.4f: U = %.6f, w_R = %.6f, E[T_U] = %.6f\n",
  sol_g2$g,
  sol_g2$U,
  sol_g2$wR,
  sol_g2$ETU
))
cat(sprintf(
  "reservation check (should be ~0): %.3e\n",
  sol_g2$check_reservation
))


############################################################################
## Q4: Simulate 1000 workers for T = 50 years, g = 0.03                   ##
############################################################################

## New parameters
set.seed(123)
N <- 1000
T <- 50

## reservation wage from part (c)
wR <- sol_g1$wR

## convert continuous-time rates to 1-year transition probabilities
p_sep <- 1 - exp(-par$lambda) # separation prob per year
p_offer <- 1 - exp(-par$alpha) # offer-arrival prob per year (while unemployed)

## state vectors
employed <- rep(FALSE, N)
wage <- rep(NA_real_, N)

for (t in 1:T) {
  ## unemployed: offer arrival and accept/reject
  idx_u <- which(!employed)
  if (length(idx_u) > 0) {
    offer_arrives <- runif(length(idx_u)) < p_offer
    idx_offer <- idx_u[offer_arrives]

    if (length(idx_offer) > 0) {
      w_draw <- rlnorm(length(idx_offer), meanlog = par$mu, sdlog = par$sigma)
      accept <- w_draw >= wR

      idx_accept <- idx_offer[accept]
      if (length(idx_accept) > 0) {
        employed[idx_accept] <- TRUE
        wage[idx_accept] <- w_draw[accept]
      }
    }
  }

  ## employed: wage growth then separation
  idx_e <- which(employed)
  if (length(idx_e) > 0) {
    wage[idx_e] <- wage[idx_e] * exp(par$g1)

    sep <- runif(length(idx_e)) < p_sep
    idx_sep <- idx_e[sep]
    if (length(idx_sep) > 0) {
      employed[idx_sep] <- FALSE
      wage[idx_sep] <- NA_real_
    }
  }
}

## counts after T years
cat(sprintf(
  "After %d years: employed = %d, unemployed = %d\n",
  T,
  sum(employed),
  N - sum(employed)
))

## histogram of employed wages
w_emp <- wage[employed]
df <- data.table(w = w_emp)

plot01 =
  ggplot(df, aes(x = w)) +
  geom_histogram(bins = 50, color = 'white', fill = 'royalblue') +
  labs(x = "Wage at t = 50", y = "Count") +
  theme_bw()

plot01

ggsave(
  plot01,
  filename = 'PS2/employed_wages_freq.png',
  height = 4,
  width = 6,
  dpi = 300
)

## Log-log survival plot for the right tail
if (length(w_emp) > 50) {
  w_sorted <- sort(w_emp)
  surv <- rev(seq_along(w_sorted)) / length(w_sorted) # empirical survival P(W >= w)
  tail_df <- data.table(w = w_sorted, surv = surv)

  plot02 =
    ggplot(tail_df, aes(x = log(w), y = log(surv))) +
    geom_point(alpha = 0.4, color = '#3737D3') +
    labs(x = "log(w)", y = "log P(W >= w)") +
    theme_bw()

  plot02

  ggsave(
    plot02,
    filename = 'PS2/log-log_survival_for_right_tail.png',
    height = 4,
    width = 6,
    dpi = 300
  )
}
