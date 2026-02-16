"""
Q5 Shooting algorithm (Python translation of ps1_q5.R)

- Implements the SAME detrended one-capital model as your R code:
    y~ = Gamma * k~^alpha * l^(1-alpha)
    k~dot = y~ - c~ - (delta + g_z) k~
    c~dot = c~ * ( alpha*(y~/k~) - delta - rho - g_z )
  with endogenous l* from intratemporal condition at SS:
    ((1-theta)/theta) * (c~/(1-l)) = (1-alpha) * (y~/l)

- Uses:
  * robust RK4 integrator (manual, like your R)
  * bracket search for l* and for c0
  * univariate root-finding (brentq) for l* and for c0

- Produces and saves two figures:
    1) "Transition Path in Detrended Levels.png"
    2) "Transition Dynamics in Detrended Time Derivatives.png"

Requirements:
    pip install numpy pandas scipy matplotlib
"""

from __future__ import annotations
import matplotlib
matplotlib.use("Agg")  # headless backend; prevents UI render calls
import matplotlib.pyplot as plt
plt.ioff()  # turn interactive mode off

from dataclasses import dataclass
import math
import numpy as np
import pandas as pd
from scipy.optimize import brentq
import matplotlib.pyplot as plt


# =========================
# Parameters (mirror the R)
# =========================

@dataclass
class Params:
    theta: float = 0.4
    beta: float = 0.95  # not used directly (R notes this too)
    alpha_e: float = 0.17
    alpha_s: float = 0.13
    alpha: float = 0.30
    delta_s: float = 0.056
    delta_e: float = 0.124
    g_gross: float = 1.0124
    g_z: float = 1.0124 - 1.0
    rho: float = -math.log(0.9)  # NOTE: matches your R script (even if beta=0.95)


def build_par() -> dict:
    par = Params().__dict__.copy()

    # single depreciation (weighted sum, as in your R)
    w_e = par["alpha_e"] / par["alpha"]
    w_s = par["alpha_s"] / par["alpha"]
    par["delta"] = w_e * par["delta_e"] + w_s * par["delta_s"]

    print("Depreciation:")
    print(f"delta_e (paper) = {par['delta_e']:.3f}")
    print(f"delta_s (paper) = {par['delta_s']:.3f}")
    print(f"delta   (PS)    = {par['delta']:.6f}\n")

    # Gamma (one-capital reduced form)
    par["Gamma"] = (par["alpha_e"] / par["alpha"]) ** par["alpha_e"] * \
                   (par["alpha_s"] / par["alpha"]) ** par["alpha_s"]

    return par


# ======================================
# Steady state with endogenous l (detr.)
# ======================================

def ktilde_from_l(l: float, par: dict) -> float:
    A = (par["delta"] + par["rho"] + par["g_z"]) / \
        (par["alpha"] * par["Gamma"] * (l ** (1.0 - par["alpha"])))
    return A ** (1.0 / (par["alpha"] - 1.0))


def yc_from_l(l: float, par: dict) -> dict:
    ktilde = ktilde_from_l(l, par)
    ytilde = par["Gamma"] * (ktilde ** par["alpha"]) * (l ** (1.0 - par["alpha"]))
    ctilde = ytilde - (par["delta"] + par["g_z"]) * ktilde
    return {"ktilde": ktilde, "ytilde": ytilde, "ctilde": ctilde}


def intra_resid(l: float, par: dict) -> float:
    if (not np.isfinite(l)) or (l <= 0.0) or (l >= 1.0):
        return np.nan
    tmp = yc_from_l(l, par)
    if (not np.isfinite(tmp["ctilde"])) or (tmp["ctilde"] <= 0.0):
        return np.nan

    lhs = (1.0 - par["theta"]) / par["theta"] * (tmp["ctilde"] / (1.0 - l))
    rhs = (1.0 - par["alpha"]) * (tmp["ytilde"] / l)
    return lhs - rhs


def find_bracket_l(par: dict, grid: np.ndarray | None = None) -> tuple[float, float]:
    if grid is None:
        grid = np.arange(0.05, 0.95 + 1e-12, 0.01)
    vals = np.array([intra_resid(x, par) for x in grid], dtype=float)
    ok = np.isfinite(vals)
    grid = grid[ok]
    vals = vals[ok]
    sgn = np.sign(vals)
    idx = np.where(sgn[1:] * sgn[:-1] < 0)[0]
    if len(idx) == 0:
        raise RuntimeError("No sign change for intratemporal residual.")
    i = idx[0]
    return float(grid[i]), float(grid[i + 1])


# ====================
# Dynamics and RK4
# ====================

def rhs(state: np.ndarray, par: dict) -> np.ndarray:
    k, c = float(state[0]), float(state[1])
    l = par["l"]
    if (not np.isfinite(k)) or (not np.isfinite(c)) or (k <= 0.0) or (c <= 0.0):
        return np.array([np.nan, np.nan], dtype=float)

    y = par["Gamma"] * (k ** par["alpha"]) * (l ** (1.0 - par["alpha"]))
    dk = y - c - (par["delta"] + par["g_z"]) * k
    dc = c * (par["alpha"] * (y / k) - par["delta"] - par["rho"] - par["g_z"])
    return np.array([dk, dc], dtype=float)


def rk4_step(state: np.ndarray, dt: float, par: dict) -> np.ndarray | None:
    k1 = rhs(state, par)
    if np.any(~np.isfinite(k1)):
        return None

    k2 = rhs(state + 0.5 * dt * k1, par)
    if np.any(~np.isfinite(k2)):
        return None

    k3 = rhs(state + 0.5 * dt * k2, par)
    if np.any(~np.isfinite(k3)):
        return None

    k4 = rhs(state + dt * k3, par)
    if np.any(~np.isfinite(k4)):
        return None

    s_next = state + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    if np.any(~np.isfinite(s_next)) or (s_next[0] <= 0.0) or (s_next[1] <= 0.0):
        return None
    return s_next


def simulate_path_rk4(
    k0: float,
    c0: float,
    par: dict,
    T_end: float = 200.0,
    dt: float = 0.001,
    k_min: float = 1e-12,
    c_min: float = 1e-12,
    k_max: float = 1e8,
    c_max: float = 1e8,
) -> tuple[pd.DataFrame, str]:
    n = int(T_end / dt)
    time = np.empty(n + 1)
    ktilde = np.empty(n + 1)
    ctilde = np.empty(n + 1)

    state = np.array([k0, c0], dtype=float)
    time[0], ktilde[0], ctilde[0] = 0.0, state[0], state[1]

    status = "success"
    last = 0

    for j in range(1, n + 1):
        # bounds stop
        if (not np.isfinite(state[0])) or (not np.isfinite(state[1])) or \
           (state[0] <= k_min) or (state[1] <= c_min) or \
           (state[0] >= k_max) or (state[1] >= c_max):
            status = "exploded" if (state[0] >= k_max or state[1] >= k_max) else "hit_lower"
            last = j - 1
            break

        nxt = rk4_step(state, dt, par)
        if nxt is None:
            status = "failed"
            last = j - 1
            break

        state = nxt
        time[j] = j * dt
        ktilde[j] = state[0]
        ctilde[j] = state[1]
        last = j

    df = pd.DataFrame({"time": time[: last + 1], "ktilde": ktilde[: last + 1], "ctilde": ctilde[: last + 1]})
    return df, status


# =========================
# Shooting objective + bracket
# =========================

def terminal_error(c0: float, k0: float, par: dict, k_star: float, T_end: float, dt: float) -> float:
    if (not np.isfinite(c0)) or (c0 <= 0.0):
        return -1e6

    out, status = simulate_path_rk4(k0, c0, par, T_end=T_end, dt=dt)
    t_last = float(out["time"].iloc[-1])
    k_last = float(out["ktilde"].iloc[-1])

    base = (k_last - k_star) / k_star

    if status != "success":
        ampl = 1.0 + 50.0 * (T_end - t_last) / T_end
        if not np.isfinite(base):
            base = -1.0
        return ampl * base

    return base


def find_c_bracket(
    k0: float,
    par: dict,
    k_star: float,
    c_star: float,
    T_end: float,
    dt: float,
    mult_grid: np.ndarray | None = None,
) -> tuple[float, float]:
    if mult_grid is None:
        mult_grid = np.exp(np.linspace(np.log(1e-6), np.log(1e6), 121))

    grid = c_star * mult_grid
    vals = np.array([terminal_error(c, k0, par, k_star, T_end, dt) for c in grid], dtype=float)

    ok = np.isfinite(vals)
    grid = grid[ok]
    vals = vals[ok]

    sgn = np.sign(vals)
    idx = np.where(sgn[1:] * sgn[:-1] < 0)[0]

    print(f"\nBracket scan: finite={len(vals)}/{len(mult_grid)}, f range=[{np.min(vals):.3g}, {np.max(vals):.3g}]")

    if len(idx) == 0:
        imn = int(np.argmin(vals))
        imx = int(np.argmax(vals))
        print(f"  min f={vals[imn]:.4g} at c0={grid[imn]:.6g}")
        print(f"  max f={vals[imx]:.4g} at c0={grid[imx]:.6g}")
        raise RuntimeError(
            "No sign change found in bracket scan. "
            "This usually indicates Euler/detrending mismatch or too-short T_end."
        )

    i = idx[0]
    return float(grid[i]), float(grid[i + 1])


# =========================
# Plot helpers (matplotlib)
# =========================

LABELS = {
    "ktilde": "Detrended Capital",
    "ktilde_e": "Detrended Equipment Capital",
    "ktilde_s": "Detrended Structure Capital",
    "ctilde": "Detrended Consumption",
    "ytilde": "Detrended Production",
    "itilde": "Detrended Investment",
    "kdot_num": "Diagnostic K dot",
    "kdot": "K dot (model)",
    "kdot_e": "Equipment K dot (model)",
    "kdot_s": "Structure K dot (model)",
    "cdot": "C dot (model)",
    "ydot": "Y dot (model)",
    "idot": "I dot (model)",
}

ORDER = [
    "ktilde", "ktilde_e", "ktilde_s", "ctilde", "ytilde", "itilde",
    "kdot_num", "kdot", "kdot_e", "kdot_s", "cdot", "ydot", "idot"
]


def facet_plot(df: pd.DataFrame, ss: dict, variables: list[str], filename: str, ncols: int = 4):
    n = len(variables)
    nrows = int(math.ceil(n / ncols))
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(2.2 * ncols, 2.0 * nrows), sharex=False)
    axes = np.array(axes).reshape(-1)

    for ax_i, var in enumerate(variables):
        ax = axes[ax_i]
        ax.plot(df["time"], df[var])
        if var in ss:
            ax.axhline(ss[var], linestyle="--")
        ax.set_title(LABELS.get(var, var), fontsize=9)
        ax.set_xlabel("Time", fontsize=8)
        ax.tick_params(axis="both", labelsize=8)

    # Hide unused axes
    for j in range(n, len(axes)):
        axes[j].axis("off")

    plt.tight_layout()
    fig.savefig(filename, dpi=300)
    plt.close(fig)


# =========================
# Main
# =========================

def main():
    par = build_par()

    # --- solve SS for l* ---
    br_l = find_bracket_l(par)
    l_star = brentq(lambda x: intra_resid(x, par), br_l[0], br_l[1], xtol=1e-12, maxiter=5000)
    tmp = yc_from_l(l_star, par)
    ktilde_star, ytilde_star, ctilde_star = tmp["ktilde"], tmp["ytilde"], tmp["ctilde"]

    print("Steady state (detrended):")
    print(f"l*      = {l_star:.10f}")
    print(f"ktilde* = {ktilde_star:.10f}")
    print(f"ytilde* = {ytilde_star:.10f}")
    print(f"ctilde* = {ctilde_star:.10f}\n")

    par["l"] = float(l_star)

    # Check SS residuals (dk, dc)
    print("SS residuals (dk, dc):")
    print(rhs(np.array([ktilde_star, ctilde_star], dtype=float), par))
    print()

    # --- shooting ---
    T_end = 40.0
    dt = 0.001
    k0 = 0.99 * ktilde_star

    br_c = find_c_bracket(k0, par, ktilde_star, ctilde_star, T_end, dt)
    print(f"Bracket for c0: [{br_c[0]:.10f}, {br_c[1]:.10f}]\n")

    c0_star = brentq(
        lambda c0: terminal_error(c0, k0, par, ktilde_star, T_end, dt),
        br_c[0], br_c[1],
        xtol=1e-10, maxiter=2000
    )

    print(f"Shooting solution c0 = {c0_star:.10f}\n")

    path, status = simulate_path_rk4(k0, c0_star, par, T_end=T_end, dt=dt)
    print(f"RK4 status: {status}  t_last={path['time'].iloc[-1]:.6f}\n")

    # --- levels + dots (analytic) + recovered ke/ks ---
    l = par["l"]
    path["ytilde"] = par["Gamma"] * (path["ktilde"] ** par["alpha"]) * (l ** (1.0 - par["alpha"]))
    path["itilde"] = path["ytilde"] - path["ctilde"]

    path["kdot_num"] = np.nan
    if len(path) >= 2:
        path.loc[path.index[1:], "kdot_num"] = np.diff(path["ktilde"].values) / dt

    path["kdot"] = path["ytilde"] - path["ctilde"] - (par["delta"] + par["g_z"]) * path["ktilde"]
    path["cdot"] = path["ctilde"] * (
        par["alpha"] * (path["ytilde"] / path["ktilde"]) - par["delta"] - par["rho"] - par["g_z"]
    )
    path["ydot"] = (par["alpha"] * path["ytilde"] / path["ktilde"]) * path["kdot"]
    path["idot"] = path["ydot"] - path["cdot"]

    w_e = par["alpha_e"] / par["alpha"]
    w_s = par["alpha_s"] / par["alpha"]
    path["ktilde_e"] = w_e * path["ktilde"]
    path["ktilde_s"] = w_s * path["ktilde"]
    path["kdot_e"] = w_e * path["kdot"]
    path["kdot_s"] = w_s * path["kdot"]

    # --- SS implied objects ---
    ktilde_e_star = w_e * ktilde_star
    ktilde_s_star = w_s * ktilde_star
    itilde_star = ytilde_star - ctilde_star

    print("Steady-state implied objects:")
    print(f"itilde*    = {itilde_star:.10f}")
    print(f"ktilde_e*  = {ktilde_e_star:.10f}")
    print(f"ktilde_s*  = {ktilde_s_star:.10f}")
    print(f"check sum  = {(ktilde_e_star + ktilde_s_star):.10f} (should equal ktilde*)\n")

    # --- terminal diagnostics ---
    kT, cT, yT, iT = (path[col].iloc[-1] for col in ["ktilde", "ctilde", "ytilde", "itilde"])
    print(f"Terminal: k(T)={kT:.10f} (rel err {(kT - ktilde_star)/ktilde_star:.3e})")
    print(f"          c(T)={cT:.10f} (rel err {(cT - ctilde_star)/ctilde_star:.3e})")
    print(f"          y(T)={yT:.10f} (rel err {(yT - ytilde_star)/ytilde_star:.3e})")
    print(f"          i(T)={iT:.10f} (rel err {(iT - itilde_star)/itilde_star:.3e})\n")

    # --- Figures (facet-style grids) ---
    ss_levels = {
        "ktilde": ktilde_star,
        "ktilde_e": ktilde_e_star,
        "ktilde_s": ktilde_s_star,
        "ctilde": ctilde_star,
        "ytilde": ytilde_star,
        "itilde": itilde_star,
    }
    ss_dots = {
        "kdot": 0.0, "cdot": 0.0, "ydot": 0.0, "idot": 0.0, "kdot_e": 0.0, "kdot_s": 0.0,
        # "kdot_num" omitted from SS lines (R also excludes it in plot02)
    }

    # Plot 1: detrended levels (exclude any *dot* series)
    vars_levels = [v for v in ORDER if ("dot" not in v)]
    facet_plot(path, ss_levels, vars_levels, "Transition Path in Detrended Levels Py.png", ncols=3)

    # Plot 2: detrended derivatives (dot series, exclude kdot_num)
    vars_dots = [v for v in ORDER if (("dot" in v) and (v != "kdot_num"))]
    facet_plot(path, ss_dots, vars_dots, "Transition Dynamics in Detrended Time Derivatives Py.png", ncols=3)

    print("Saved figures:")
    print(" - Transition Path in Detrended Levels.png")
    print(" - Transition Dynamics in Detrended Time Derivatives.png")


if __name__ == "__main__":
    main()
