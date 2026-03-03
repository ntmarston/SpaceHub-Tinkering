#!/usr/bin/env python3
"""
generate_pagn_disk.py — Generate a disk CSV from pagn (Sirko & Goodman 2003)
in the format consumed by SpaceHub's disk-model.hpp.

Key differences from disktab.py:
  - pagn uses tabulated opacities (Semenov 2003 + Badnell 2005 combined),
    disktab uses analytic Kramers + electron scattering
  - pagn uses mu ~ 1 (ct.massU = 1.66e-27 kg), disktab uses mu = 0.6
  - pagn uses Sigma = 2*rho*H (slab model) everywhere;
    disktab uses sqrt(2*pi)*rho*H (Gaussian) in the standard zone

Usage:
    python generate_pagn_disk.py --M_msun 1e8 --alpha 0.01 --le 1.0 \\
        --output ../../SpaceHub/src/interaction/disk_tab/disk_stae321_pagn.csv
"""
#Claude was used here to:
# write boilerplate / CLI
# structure the conversion pipeline

import argparse
import io
import sys
from contextlib import redirect_stdout

import numpy as np
import pandas as pd
from pagn import Sirko
import pagn.constants as ct


# ── SpaceHub unit conversions (must match disktab.py spacehub_pretab()) ──────
AU_m = 1.495978707e11        # meters per AU
Msun_kg = 1.98847e30         # kg per solar mass
year_s = 365.25636042 * 24 * 3600  # seconds per year
T_unit = year_s / (2 * np.pi)      # SpaceHub time unit in seconds

CONV_LENGTH          = 1 / AU_m                      # m → AU
CONV_DENSITY         = AU_m**3 / Msun_kg             # kg/m³ → Msun/AU³
CONV_PRESSURE        = AU_m * T_unit**2 / Msun_kg    # Pa → Msun/(AU·T²)
CONV_VELOCITY        = T_unit / AU_m                  # m/s → AU/T
CONV_VISCOSITY       = T_unit / AU_m**2              # m²/s → AU²/T
CONV_SURFACE_DENSITY = AU_m**2 / Msun_kg             # kg/m² → Msun/AU²


def log_gradient(y, x):
    """Compute -d(ln y)/d(ln x) using central differences."""
    return -np.gradient(np.log(y), np.log(x))


def _smooth_zone_gradients(df, n_gap=100, w=100):
    """
    Smooth grad_P and grad_Sigma at the Standard/Self-Reg zone boundary
    using a C¹-continuous cubic Hermite interpolant.

    Replaces 2*n_gap rows centred on the boundary with a cubic bridge
    that matches both value and slope at the anchor points.
    Original values are saved in grad_P_raw / grad_Sigma_raw columns.

    Parameters
    ----------
    df    : pd.DataFrame  — must have: R/Rg, zone, grad_P, grad_Sigma.
    n_gap : int           — rows removed on each side of boundary (default 100).
    w     : int           — half-window for slope estimation (default 100).
    """
    import warnings
    from scipy.interpolate import CubicHermiteSpline

    sr_idx = np.where(df['zone'].values == 'Self-Reg')[0]
    if len(sr_idx) == 0:
        return

    ibnd  = int(sr_idx[0])
    log_R = np.log10(df['R/Rg'].values)
    N     = len(df)
    i_lo  = ibnd - n_gap
    i_hi  = ibnd + n_gap

    if i_lo <= w or i_hi + w >= N:
        warnings.warn(
            f"Zone boundary smoothing skipped: gap [{i_lo}, {i_hi}) with slope "
            f"window w={w} exceeds data range [0, {N}). Reduce n_gap or w.",
            stacklevel=2,
        )
        return

    for col in ('grad_P', 'grad_Sigma'):
        y = df[col].values.copy()
        df[f'{col}_raw'] = y.copy()

        x0, y0 = log_R[i_lo - 1], y[i_lo - 1]
        x1, y1 = log_R[i_hi],     y[i_hi]
        dy0 = np.mean(np.gradient(y[i_lo - w - 1 : i_lo],         log_R[i_lo - w - 1 : i_lo]))
        dy1 = np.mean(np.gradient(y[i_hi         : i_hi + w + 1], log_R[i_hi : i_hi + w + 1]))

        chs = CubicHermiteSpline([x0, x1], [y0, y1], [dy0, dy1])
        y[i_lo:i_hi] = chs(log_R[i_lo:i_hi])
        df[col] = y

    warnings.warn(
        f"Zone boundary smoothing applied (CubicHermiteSpline, n_gap={n_gap}, w={w}): "
        f"{2 * n_gap} rows [{i_lo}, {i_hi}) around index {ibnd} "
        f"(R/Rg ≈ {df['R/Rg'].iloc[ibnd]:.1f}) replaced for grad_P and grad_Sigma. "
        f"Originals preserved in grad_P_raw / grad_Sigma_raw.",
        stacklevel=2,
    )


def generate_pagn_csv(M_msun, alpha, le, eps, N, output_path):
    """Generate a pagn disk model and write it as a SpaceHub-format CSV."""

    Mbh = M_msun * ct.MSun
    Rg = ct.G * Mbh / (ct.c ** 2)

    # Suppress pagn's verbose init output
    buf = io.StringIO()
    with redirect_stdout(buf):
        sk = Sirko.SirkoAGN(Mbh=Mbh, alpha=alpha, le=le, eps=eps)
        sk.solve_disk(N=int(N))

    # Print pagn's summary to stderr so it's visible but doesn't pollute stdout
    print(buf.getvalue(), file=sys.stderr)

    # ── Derived quantities in SI ─────────────────────────────────────────────
    R       = sk.R          # meters
    R_Rg    = R / Rg
    Tc      = sk.T          # Kelvin
    rho     = sk.rho        # kg/m³
    cs      = sk.cs         # m/s
    H       = sk.h          # meters
    Q       = sk.Q          # dimensionless

    P     = rho * cs**2                  # Pa (cs² = P_total/rho by pagn convention)
    visc  = alpha * cs * H               # m²/s
    Sigma = 2 * rho * H                  # kg/m² (slab model, pagn convention)

    # ── Zone labels ──────────────────────────────────────────────────────────
    isf = sk.isf  # index where star formation (Q=1) begins; -1 if never
    if isf > 0:
        zone = np.array(["Standard"] * len(R))
        zone[isf:] = "Self-Reg"
    else:
        zone = np.array(["Standard"] * len(R))

    # ── Log gradients (zone-by-zone to avoid boundary artifacts) ─────────────
    grad_T     = np.zeros_like(R)
    grad_Sigma = np.zeros_like(R)
    grad_P     = np.zeros_like(R)

    if isf > 0 and isf < len(R) - 1:
        # Standard zone (indices 0 to isf-1, inclusive)
        std = slice(0, isf)
        grad_T[std]     = log_gradient(Tc[std],    R[std])
        grad_Sigma[std] = log_gradient(Sigma[std], R[std])
        grad_P[std]     = log_gradient(P[std],     R[std])

        # Self-Reg zone (indices isf to end)
        sr = slice(isf, None)
        grad_T[sr]     = log_gradient(Tc[sr],    R[sr])
        grad_Sigma[sr] = log_gradient(Sigma[sr], R[sr])
        grad_P[sr]     = log_gradient(P[sr],     R[sr])
    else:
        # Single zone
        grad_T     = log_gradient(Tc,    R)
        grad_Sigma = log_gradient(Sigma, R)
        grad_P     = log_gradient(P,     R)

    # ── Convert to SpaceHub units ────────────────────────────────────────────
    df = pd.DataFrame({
        "R":          R * CONV_LENGTH,
        "R/Rg":       R_Rg,
        "Tc":         Tc,                          # K (unchanged)
        "rho":        rho * CONV_DENSITY,
        "P":          P * CONV_PRESSURE,
        "cs":         cs * CONV_VELOCITY,
        "H":          H * CONV_LENGTH,
        "visc":       visc * CONV_VISCOSITY,
        "Sigma":      Sigma * CONV_SURFACE_DENSITY,
        "Q":          Q,
        "grad_T":     grad_T,
        "grad_Sigma": grad_Sigma,
        "grad_P":     grad_P,
        "zone":       zone,
    })

    # Smooth grad_P / grad_Sigma across the zone boundary (n_gap=100, w=100)
    if isf > 0:
        _smooth_zone_gradients(df, n_gap=100, w=100)

    # Reorder: 13 numeric base columns → zone → _raw columns.
    # C++ reads by position and stops at grad_P (col 13), so _raw cols
    # placed after zone are invisible to disk-model.hpp.
    raw_cols  = [c for c in df.columns if c.endswith('_raw')]
    base_cols = [c for c in df.columns if c != 'zone' and not c.endswith('_raw')]
    df = df[base_cols + ['zone'] + raw_cols]

    df.to_csv(output_path, index=False)
    print(f"Wrote {len(df)} rows to {output_path}")
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate pagn (Sirko-Goodman) disk CSV for SpaceHub."
    )
    parser.add_argument("--M_msun", type=float, default=1e8,
                        help="Central BH mass in solar masses (default: 1e8)")
    parser.add_argument("--alpha", type=float, default=0.01,
                        help="Viscosity parameter (default: 0.01)")
    parser.add_argument("--le", type=float, default=1.0,
                        help="Eddington ratio L/L_Edd (default: 1.0)")
    parser.add_argument("--eps", type=float, default=0.1,
                        help="Radiative efficiency (default: 0.1)")
    parser.add_argument("--N", type=int, default=5000,
                        help="Number of radial grid points (default: 5000)")
    parser.add_argument("--output", type=str, required=True,
                        help="Output CSV file path")

    args = parser.parse_args()
    generate_pagn_csv(args.M_msun, args.alpha, args.le, args.eps, args.N, args.output)
