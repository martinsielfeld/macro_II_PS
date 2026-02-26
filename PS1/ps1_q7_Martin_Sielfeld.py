"""
Q7 steady state with increasing q(t) — Python version matching your UPDATED R parametrization.

Implements exactly the same objects as the R code you pasted:

Outputs:
  - prints the same steady state objects as your R print block
  - optional: saves a residual-vs-l plot and a bar chart into ./figures_q7/
"""

from __future__ import annotations

import os
from dataclasses import dataclass
import numpy as np
from scipy.optimize import brentq

# Headless-safe plotting (prevents UI/render timeouts)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


@dataclass
class Params:
    theta: float = 0.4
    ae: float = 0.17
    aS: float = 0.13  # "as" is a Python keyword-ish in some contexts, so use aS
    de: float = 0.124
    ds: float = 0.056
    rho: float = -np.log(0.9)
    gz: float = 0.0124
    gq: float = np.log(1.032)

    make_plots: bool = True
    outdir: str = "figures_q7"


def solve_q7(par: Params) -> dict:
    theta, ae, aS = par.theta, par.ae, par.aS
    alpha = ae + aS

    # Representative depreciation
    w_e = ae / alpha
    w_s = aS / alpha
    d = w_e * par.de + w_s * par.ds

    # Growth rates
    gy = par.gz + (ae / (1.0 - alpha)) * par.gq

    # Constants (with d used everywhere)
    B = (ae / (d + par.rho + gy - par.gq)) ** ae * (aS / (d + par.rho + gy)) ** aS
    A = B ** (1.0 / (1.0 - alpha))

    # Functions
    def ytilde(l: float) -> float:
        return A * l

    def ks_from_l(l: float) -> float:
        return (aS / (d + par.rho + gy)) * ytilde(l)

    def ke_from_l(l: float) -> float:
        return (ae / (d + par.rho + gy - par.gq)) * ytilde(l)

    def ctilde_from_l(l: float) -> float:
        return (theta / (1.0 - theta)) * (1.0 - alpha) * ytilde(l) * (1.0 - l) / l

    def resid_l(l: float) -> float:
        y = ytilde(l)
        ks = ks_from_l(l)
        ke = ke_from_l(l)
        c = ctilde_from_l(l)
        i_s = (d + gy) * ks
        i_e = (d + gy + par.gq) * ke
        return y - c - i_s - i_e

    # Bracket root safely by scanning (more robust than assuming uniroot bracket works)
    grid = np.linspace(1e-6, 1 - 1e-6, 4000)
    vals = np.array([resid_l(x) for x in grid])
    sgn = np.sign(vals)
    idx = np.where(sgn[:-1] * sgn[1:] < 0)[0]
    if len(idx) == 0:
        raise RuntimeError(
            "Could not bracket root for l in (0,1). "
            f"Residual range: [{vals.min():.6g}, {vals.max():.6g}]"
        )
    a, b = float(grid[idx[0]]), float(grid[idx[0] + 1])

    l_star = brentq(resid_l, a, b, xtol=1e-12, maxiter=5000)

    # SS values
    y_star = ytilde(l_star)
    ks_star = ks_from_l(l_star)
    ke_star = ke_from_l(l_star)
    c_star = ctilde_from_l(l_star)
    is_star = (d + gy) * ks_star
    ie_star = (d + gy + par.gq) * ke_star

    out = {
        # Params/derived
        "theta": float(theta),
        "ae": float(ae),
        "as": float(aS),
        "alpha": float(alpha),
        "de": float(par.de),
        "ds": float(par.ds),
        "w_e": float(w_e),
        "w_s": float(w_s),
        "d": float(d),
        "rho": float(par.rho),
        "gz": float(par.gz),
        "gq": float(par.gq),
        "gy": float(gy),
        "B": float(B),
        "A": float(A),
        # Steady state
        "l_star": float(l_star),
        "y_star": float(y_star),
        "ks_star": float(ks_star),
        "ke_star": float(ke_star),
        "c_star": float(c_star),
        "is_star": float(is_star),
        "ie_star": float(ie_star),
        "resource_residual": float(y_star - c_star - is_star - ie_star),
        # For plotting
        "_grid": grid,
        "_vals": vals,
    }
    return out


def maybe_save_plots(res: dict, outdir: str) -> None:
    os.makedirs(outdir, exist_ok=True)

    # Residual vs l
    fig1 = plt.figure(figsize=(7, 4))
    plt.plot(res["_grid"], res["_vals"])
    plt.axhline(0.0, linestyle="--")
    plt.axvline(res["l_star"], linestyle="--")
    plt.xlabel("l")
    plt.ylabel("resid(l) = y~ - c~ - i_s~ - i_e~")
    plt.title("Q7 steady-state residual")
    fig1.savefig(os.path.join(outdir, "q7_residual_vs_l.png"), dpi=200, bbox_inches="tight")
    plt.close(fig1)

    # Bars of SS quantities
    labels = ["Ks~", "Ke~", "C~", "Y~", "Is~", "Ie~"]
    values = [res["ks_star"], res["ke_star"], res["c_star"], res["y_star"], res["is_star"], res["ie_star"]]
    fig2 = plt.figure(figsize=(7, 4))
    plt.bar(labels, values)
    plt.title("Q7 detrended steady state")
    fig2.savefig(os.path.join(outdir, "q7_steady_state_bars.png"), dpi=200, bbox_inches="tight")
    plt.close(fig2)


def main():
    par = Params()
    res = solve_q7(par)

    print("Representative depreciation:")
    print(f"  w_e = {res['w_e']:.6f}")
    print(f"  w_s = {res['w_s']:.6f}")
    print(f"  d   = {res['d']:.6f}\n")

    print("Detrended Ks ss:", res["ks_star"])
    print("Detrended Ke ss:", res["ke_star"])
    print("Detrended C  ss:", res["c_star"])
    print("Detrended Y  ss:", res["y_star"])
    print("Detrended L  ss:", res["l_star"])
    print("Detrended Is ss:", res["is_star"])
    print("Detrended Ie ss:", res["ie_star"])
    print("Resource residual:", res["resource_residual"])

    if par.make_plots:
        maybe_save_plots(res, par.outdir)
        print(f"\nSaved plots to ./{par.outdir}/")


if __name__ == "__main__":
    main()
