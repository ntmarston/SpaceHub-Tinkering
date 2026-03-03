#!/usr/bin/env python3
"""
generate_pure_ss_disk.py — Generate a pure Shakura-Sunyaev disk CSV (no Q=1 zone)
for SpaceHub's disk-model.hpp.

This is test B7: testing whether the absence of the S&G zone boundary discontinuity
in grad_P resolves the Bulirsch-Stoer integrator hang at ~172,000 yr.

Parameters match the existing pagn disk (disk_stae321_pagn.csv):
  M = 1e8 Msun, alpha = 0.01, le = 1.0, eps = 0.1

Usage:
    python3 generate_pure_ss_disk.py
"""

import sys
import os

# Add disktab module to path
sys.path.insert(0, os.path.abspath('../../disktab'))

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.constants import G, c, m_p, sigma_T
import disktab as dt

OUTPUT = '../../SpaceHub/src/interaction/disk_tab/disk_pure_ss_1e8.csv'

# ── Parameters (matching disk_stae321_pagn.csv) ──────────────────────────────
M_msun = 1e8
alpha  = 0.01
le     = 1.0   # Eddington ratio L/L_Edd
eps    = 0.1   # radiative efficiency

# ── Build disk object ─────────────────────────────────────────────────────────
disk = dt.AGNDisk()
disk.M     = M_msun * u.Msun
disk.alpha = alpha

# Eddington accretion rate
kappa    = sigma_T / m_p
L_Edd    = 4 * np.pi * G * disk.M * c / kappa
Mdot_Edd = (L_Edd / (eps * c**2)).to(u.kg / u.s)
disk.Mdot = le * Mdot_Edd

print(f"M     = {M_msun:.2e} Msun")
print(f"Mdot  = {disk.Mdot:.4e}")
print(f"alpha = {alpha}")
print(f"le    = {le}, eps = {eps}")

# ── Generate pure Shakura-Sunyaev solution (no Q=1 zone) ─────────────────────
rvals = np.logspace(1, 7, 1000)   # 10 to 10^7 Rg, 1000 log-spaced points
print(f"\nGenerating pure S-S disk ({len(rvals)} points, 10–10^7 Rg) ...")
disk.generate(rvals=rvals, sirko_goodman=False)

# CRITICAL: drop 'tau' column before spacehub_pretab().
# When sirko_goodman=False, _full_disk_solution() keeps 'tau' in the DataFrame.
# spacehub_pretab() does NOT drop it, so without this step the CSV would have
# a 'tau' column between 'Sigma' and 'Q', breaking disk-model.hpp's CSV parser
# (it would misread tau as grad_T, grad_T as grad_Sigma, grad_Sigma as grad_P).
if 'tau' in disk.model.columns:
    disk.model = disk.model.drop(columns=['tau'])
    print("Dropped 'tau' column from model.")

# ── Convert to SpaceHub units ─────────────────────────────────────────────────
# No rvals argument → uses existing model (prevents re-calling generate() with
# default sirko_goodman=True which would undo what we just did).
df = disk.spacehub_pretab()

# ── Sanity checks ─────────────────────────────────────────────────────────────
expected_cols = ['R', 'R/Rg', 'Tc', 'rho', 'P', 'cs', 'H', 'visc',
                 'Sigma', 'Q', 'grad_T', 'grad_Sigma', 'grad_P']
numeric_cols  = [c for c in df.columns if c != 'zone']

print(f"\nColumns:        {list(df.columns)}")
print(f"Numeric cols:   {len(numeric_cols)} (expected 13)")
assert numeric_cols == expected_cols, \
    f"Column mismatch!\n  Got:      {numeric_cols}\n  Expected: {expected_cols}"

nan_count = df[numeric_cols].isna().sum().sum()
inf_count = np.isinf(df[numeric_cols].values).sum()
print(f"NaN count:      {nan_count}  (expected 0)")
print(f"Inf count:      {inf_count}  (expected 0)")
assert nan_count == 0, "NaN values in output!"
assert inf_count == 0, "Inf values in output!"

print(f"\nR range:        {df['R'].iloc[0]:.2f} – {df['R'].iloc[-1]:.2e} AU")
print(f"Q range:        {df['Q'].min():.4f} – {df['Q'].max():.4e}")
print(f"grad_P range:   {df['grad_P'].min():.4f} – {df['grad_P'].max():.4f}")
max_jump = df['grad_P'].diff().abs().max()
print(f"grad_P max consecutive jump: {max_jump:.6f}  (pagn had ~0.43)")

# ── Write CSV ─────────────────────────────────────────────────────────────────
df.to_csv(OUTPUT, index=False)
print(f"\nWrote {len(df)} rows to {OUTPUT}")
print("Done.")
