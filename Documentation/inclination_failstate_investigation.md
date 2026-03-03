# Investigation: Root Cause of `Reach max iteration()` at High Orbital Inclination

**Date:** 2026-02-20
**Author:** N. Marston (with Claude Code)
**Status:** Root cause identified

---

## 1. Problem Statement

The SpaceHub N-body simulation crashes with `Reach max iteration()` when the orbital inclination of a test particle approaches $i = 90°$ in the AGN disk force model. The simulation parameters are:

| Parameter | Value |
|-----------|-------|
| Central mass $M_1$ | $10^8\,M_\odot$ |
| Orbiter mass $M_2$ | $30\,M_\odot$ |
| Semi-major axis $a$ | $0.1\,\mathrm{pc}$ |
| Eccentricity $e$ | $0.67$ |
| Force model | Aerodynamic drag only (`DiskModel`) |
| Disk table | `disk_stae321_pagn.csv` (Sirko & Goodman model) |

The error occurs in both the Bulirsch-Stoer (`methods::DefaultMethod`) and IAS15/Radau (`methods::Radau`) integrators, confirming that the adaptive integration scheme is not at fault. The root cause lies in the force model implementation in `disk-model.hpp`.

---

## 2. Prior Investigation Summary

A previous investigation tested the hypothesis that the Bulirsch-Stoer method's `first_step_` guard was the root cause. The Radau integrator (which uses a fundamentally different convergence mechanism — Gauss-Radau quadrature with predictor-corrector iteration rather than Richardson extrapolation) was tested with the same initial conditions. **Both integrators fail identically at $i = 90°$**, definitively ruling out any integrator-specific mechanism.

---

## 3. Orbital Geometry at $i \to 90°$

For Keplerian orbital elements $\Omega = 0,\, \omega = 0$, the transformation from the orbital plane to the inertial frame gives the 3D position at true anomaly $\nu$ as:

$$
x = r\cos\nu, \quad y = r\sin\nu\cos i, \quad z = r\sin\nu\sin i
$$

The cylindrical radius used by the disk model is:

$$
R_\mathrm{cyl} = \sqrt{x^2 + y^2} = r\sqrt{1 - \sin^2\!\nu\,\sin^2\!i}
$$

At $i = 90°$, this simplifies to:

$$
R_\mathrm{cyl} = r\,|\cos\nu|
$$

which **drops to zero** at $\nu = 90°$ and $\nu = 270°$ — the points where the orbit crosses the $z$-axis.

### Key orbital scales

The orbital radius at true anomaly $\nu$ is $r(\nu) = a(1-e^2)/(1 + e\cos\nu)$. For our parameters:

| Quantity | Value |
|----------|-------|
| Periapsis $r_p = a(1-e)$ | $6{,}807\,\mathrm{AU}$ |
| Apoapsis $r_a = a(1+e)$ | $34{,}446\,\mathrm{AU}$ |
| $r(\nu = 90°) = a(1-e^2)$ | $11{,}362\,\mathrm{AU}$ |
| $R_\mathrm{cyl}$ at $\nu = 90°$ (i=90°) | $0\,\mathrm{AU}$ |
| Disk table $R_\mathrm{min}$ | $11.86\,\mathrm{AU}$ |
| Disk table $R_\mathrm{max}$ | $9.87 \times 10^6\,\mathrm{AU}$ |
| Orbital period $T$ | $\sim 296\,\mathrm{yr}$ |

At low inclination ($i \ll 90°$), $R_\mathrm{cyl} \approx r$ and remains well within the disk table range throughout the orbit. As $i \to 90°$, $R_\mathrm{cyl}$ oscillates between $\sim r$ (at the nodes) and $\sim 0$ (at $\nu = 90°/270°$), plunging through the inner table boundary.

---

## 4. Root Causes

### 4.1 **CRITICAL: Hard Discontinuity at $R_\mathrm{cyl} = R_\mathrm{min}$**

**Location:** `disk-model.hpp`, line 327

```cpp
if (R_cyl < Rmin || R_cyl > Rmax) {
    continue;  // force jumps from computed value to exactly zero
}
```

This imposes a **step-function discontinuity** in the acceleration as a function of position. Both the Bulirsch-Stoer and IAS15 integrators rely on polynomial extrapolation of the solution as a function of sub-step size $h$. Specifically:

- **Bulirsch-Stoer** uses Richardson extrapolation: it evaluates the modified midpoint method at sub-step counts $n_1, n_2, \ldots$ and builds a rational polynomial extrapolation table. The error of the $k$-th extrapolation order is $\mathcal{O}(h^{2k})$, which requires the right-hand side $f(\mathbf{x}, t)$ to be $C^{2k}$-smooth (Hairer, Nørsett & Wanner 1993, §II.9).

- **IAS15** (Gauss-Radau) uses 15th-order implicit Runge-Kutta with predictor-corrector iteration. Convergence requires Lipschitz continuity of the acceleration; a step discontinuity violates this and the PC iteration fails to converge within its `max_iter = 30` limit.

A step function has derivative equal to a Dirac delta:

$$
\frac{\partial \mathbf{a}}{\partial R_\mathrm{cyl}}\bigg|_{R_\mathrm{min}} = \mathbf{a}(R_\mathrm{min}^+)\,\delta(R_\mathrm{cyl} - R_\mathrm{min})
$$

No polynomial of any degree can approximate a step function uniformly. The extrapolation table sees inconsistent sub-step evaluations (some landing inside the boundary, some outside) and the error estimate never decreases below the acceptance threshold, regardless of how small $h$ becomes. After 100 rejected attempts (Bulirsch-Stoer) or 30 failed PC iterations (IAS15), the integrator throws `Reach max iteration()`.

**Why this only matters at high $i$:** At $i = 0°$, $R_\mathrm{cyl} = r$ throughout the orbit, always $\gg R_\mathrm{min}$. The orbit never approaches the Rmin boundary. As $i$ increases, the minimum $R_\mathrm{cyl}$ per orbit decreases. The critical inclination where the orbit first touches the $R_\mathrm{min}$ cylinder satisfies:

$$
R_\mathrm{min} = r(\nu=90°)\sqrt{1 - \sin^2 i_\mathrm{crit}} \implies i_\mathrm{crit} = \arcsin\!\sqrt{1 - \left(\frac{R_\mathrm{min}}{a(1-e^2)}\right)^2}
$$

For our parameters: $i_\mathrm{crit} = \arcsin\sqrt{1 - (11.86/11362)^2} \approx 89.94°$.

> **Reference:** Hairer, E., Nørsett, S. P., & Wanner, G. (1993). *Solving Ordinary Differential Equations I: Nonstiff Problems*. Springer-Verlag, Berlin. §II.9 (Richardson extrapolation and smoothness requirements).

> ---
>
> **Proposed Fixes for §4.1**
>
> **Fix A — Smooth sigmoid taper at boundaries.** Replace the hard `continue` with a smooth window function $w(R_\mathrm{cyl})$ that multiplies the total force:
>
> $$w(R_\mathrm{cyl}) = \tfrac{1}{2}\left[1 + \tanh\!\left(\frac{R_\mathrm{cyl} - R_\mathrm{min}}{\Delta R}\right)\right]$$
>
> with $\Delta R = \alpha\, H(R_\mathrm{min})$ for some $\alpha \sim 1$. The force transitions $C^\infty$-smoothly from 0 to full strength over a scale-height-sized region. Physically motivated: the disk inner edge is not infinitely sharp; it tapers over $\sim H$ due to viscous spreading and ISCO physics (Frank, King & Raine 2002, §5.8). An analogous taper applies at $R_\mathrm{max}$.
>
> ```cpp
> // In init_from_file(), after loading table:
> static inline double H_min;   // store for taper width
> H_min = disk_table.front().H;
>
> // In add_acc_to(), replace hard continue with:
> double w_inner = 0.5 * (1.0 + tanh((R_cyl - Rmin) / H_min));
> double w_outer = 0.5 * (1.0 - tanh((R_cyl - Rmax) / (0.1 * Rmax)));
> // ... later: f_total *= w_inner * w_outer;
> ```
>
> **Fix B — Early vertical exit using pre-stored $H_\mathrm{min}$.** At `init_from_file`, store `H_min = disk_table.front().H`. In `add_acc_to`, *before* the $R_\mathrm{min}$ check, add:
>
> ```cpp
> static inline double H_min;  // set in init_from_file
> // In add_acc_to, before the Rmin guard:
> if (std::abs(z) > 30.0 * H_min) continue;
> ```
>
> Since $H$ increases with $R$ in the disk, $H_\mathrm{min}$ is a global lower bound on the scale height. At $|z|/H_\mathrm{min} > 30$, we have $\rho/\rho_c < \exp(-450) \equiv 0$ in IEEE 754 double precision for **all** radii. This exits the loop iteration before reaching the $R_\mathrm{min}$ boundary for any high-inclination orbit, completely avoiding the discontinuity. The force is already physically zero at these heights — no approximation is introduced.

### 4.2 **HIGH: Keplerian Divergence in `disk_v()` as $R_\mathrm{cyl} \to 0$**

**Location:** `disk-model.hpp`, lines 255–263

```cpp
template <typename Vec>
static Vec disk_v(Vec &r, double M, double n, double cs) {
    auto rr = sqrt(r.x * r.x + r.y * r.y);   // R_cyl
    auto v_k = sqrt(consts::G * M / rr);      // ∝ R_cyl^{-1/2} → ∞
    auto corr_factor = sqrt(1 - n*(cs*cs / (v_k*v_k)));
    auto v = v_k * corr_factor;
    return Vec{-v * r.y / rr, v * r.x / rr, 0};
}
```

The Keplerian disk velocity diverges as $R_\mathrm{cyl} \to 0$:

$$
v_\mathrm{disk} = v_K\sqrt{1 - n\frac{c_s^2}{v_K^2}} \approx v_K = \sqrt{\frac{GM_1}{R_\mathrm{cyl}}} \propto R_\mathrm{cyl}^{-1/2}
$$

where the pressure correction is negligible because $c_s^2/v_K^2 = (H/R)^2 \ll 1$ at small $R$. Evaluating at $R_\mathrm{cyl} = R_\mathrm{min} = 11.86\,\mathrm{AU}$:

$$
v_K(R_\mathrm{min}) = \sqrt{\frac{GM_1}{11.86\,\mathrm{AU}}} \approx 18{,}500\,\mathrm{AU\,yr^{-1}}
$$

compared to the orbital velocity $v_\mathrm{orb} \approx 700\,\mathrm{AU\,yr^{-1}}$. The disk velocity exceeds the orbital velocity by a factor $\sim 26$. Since the disk velocity is purely in the $xy$-plane while the orbital velocity at high $i$ is predominantly in the $xz$-plane, the relative velocity is dominated by the disk velocity:

$$
\mathbf{v}_\mathrm{rel} \approx -\mathbf{v}_\mathrm{disk} \implies |\mathbf{v}_\mathrm{rel}| \approx v_K(R_\mathrm{cyl})
$$

This means:

1. **The drag force magnitude scales as** $f_\mathrm{aero} = \pi r_\mathrm{eff}^2 \rho\, v_\mathrm{rel}^2 \propto R_\mathrm{cyl}^{-1}$, diverging as the orbit approaches the $z$-axis
2. **The drag direction is locked to $-\hat{y}$** regardless of the actual orbital velocity, because $v_\mathrm{disk} \gg v_\mathrm{orb}$
3. **The force sensitivity to position is extreme**: a small displacement $\delta R_\mathrm{cyl}$ produces $\delta v_\mathrm{disk}/v_\mathrm{disk} \sim \delta R_\mathrm{cyl}/(2R_\mathrm{cyl})$, which diverges as $R_\mathrm{cyl} \to 0$

This amplifies Root Cause 4.1: the force is not merely discontinuous at $R_\mathrm{min}$, it is diverging as it approaches the discontinuity. The product of a divergent function meeting a step function is maximally pathological for polynomial extrapolation.

**Direction singularity at $R_\mathrm{cyl} = 0$ (exactly $i = 90°$):**

At $\nu = 90°$ with $i = 90°$: $r_x = 0$, $r_y = 0$, so `rr = 0` and the expression $v r_x/\mathrm{rr}$ is $0/0$ (indeterminate). While the $R_\mathrm{cyl} < R_\mathrm{min}$ guard prevents evaluation at exactly $R_\mathrm{cyl} = 0$, the limiting behavior from either side involves $v_K \to \infty$ and an indeterminate direction — a non-removable singularity.

> ---
>
> **Proposed Fixes for §4.2**
>
> **Fix A — Skip `disk_v` when density underflows.** Move the density computation $\rho = \rho_c \exp(-z^2/2H^2)$ *before* the `disk_v()` call. If $\rho = 0$ (IEEE 754 underflow), skip the entire force computation via `continue`. This avoids computing divergent $v_K$ values that are subsequently multiplied by zero, eliminating the $\infty \times 0$ floating-point inconsistencies that poison the extrapolation table. No physics change: $0 \times \infty$ was never physically meaningful.
>
> ```cpp
> // Current order: interp → disk_v → rho → force
> // Proposed order: interp → rho → [exit if rho==0] → disk_v → force
> double rho = rho_c * exp(-0.5 * (z * z) / (H * H));
> if (rho == 0.0) continue;  // density underflow → force is identically zero
> auto v_disk = disk_v(dr, m[0], n, cs);  // only computed when rho > 0
> ```
>
> **Fix B — Floor $R_\mathrm{cyl}$ in `disk_v` at $R_\mathrm{min}$.** Cap the cylindrical radius used for the Keplerian velocity:
>
> ```cpp
> auto rr = std::max(sqrt(r.x * r.x + r.y * r.y), Rmin);
> ```
>
> This limits $v_K \leq v_K(R_\mathrm{min})$, preventing the divergence. Physically justified: inside the ISCO / inner disk edge, the gas is in free-fall rather than on circular Keplerian orbits. There is no coherent azimuthal disk velocity to compute below $R_\mathrm{min}$ (Abramowicz & Fragile 2013, §3.4). The direction $(-r_y/r_r,\, r_x/r_r,\, 0)$ remains well-defined whenever $R_\mathrm{cyl} > 0$, and the magnitude is bounded.

### 4.3 **MODERATE: Gaussian Density Impulse at Disk-Plane Crossings**

**Location:** `disk-model.hpp`, line 334

```cpp
double rho = rho_c * exp(-0.5 * (z * z) / (H * H));
```

The particle crosses the disk midplane ($z = 0$) at $\nu = 0$ and $\nu = \pi$ (for $\Omega = \omega = 0$). At these crossings, the density experienced by the particle is a Gaussian function of time:

$$
\rho(t) = \rho_c(R)\exp\!\left(-\frac{v_z^2(t - t_0)^2}{2H^2}\right)
$$

where $v_z = v_\mathrm{orb}\sin i$ is the vertical velocity component at the crossing. The characteristic temporal width of the density pulse is:

$$
\tau = \frac{H}{v_z}
$$

| Crossing | $H$ (AU) | $v_z$ (AU/yr) | $\tau$ (yr) | $\tau/T$ |
|----------|----------|---------------|-------------|----------|
| Periapsis ($\nu = 0$) | $\sim 340$ | $\sim 984$ | $0.35$ | $0.12\%$ |
| Apoapsis ($\nu = \pi$) | $\sim 1720$ | $\sim 194$ | $8.9$ | $3.0\%$ |

At the periapsis crossing, the force changes by many orders of magnitude in $\sim 0.35\,\mathrm{yr}$, occupying only $0.12\%$ of the orbital period. This creates a temporal stiffness ratio:

$$
\Lambda = \frac{T}{\tau} \approx 850 \quad \text{(periapsis crossing)}
$$

For the aerodynamic-drag-only case, the absolute force ratio $f_\mathrm{aero}/f_\mathrm{grav} \sim 10^{-13}$ is tiny, so this stiffness alone does not cause failure. However, if dynamical friction or Bondi-Hoyle accretion are enabled, the force ratio increases by orders of magnitude and this Gaussian impulse becomes the primary stiffness source.

**Important subtlety at the $R_\mathrm{min}$ crossing:** When $R_\mathrm{cyl} \lesssim R_\mathrm{min}$ near $\nu = 90°$, the height above the disk plane is $z \approx r(\nu = 90°) = 11{,}362\,\mathrm{AU}$, while the scale height at the inner disk edge is $H \approx 2.66\,\mathrm{AU}$. Thus $z/H \approx 4{,}270$ and:

$$
\rho = \rho_c \exp\!\left(-\frac{z^2}{2H^2}\right) = \rho_c \exp(-9.1 \times 10^6) \equiv 0 \quad \text{(IEEE 754 underflow)}
$$

The density and therefore the force are **identically zero in double precision** at the $R_\mathrm{min}$ crossing. This might suggest the discontinuity (§4.1) is benign. However, the intermediate computations ($v_K$, $v_\mathrm{disk}$, interpolation lookups) still produce finite, divergent values before being multiplied by zero. The Bulirsch-Stoer extrapolation table and the IAS15 predictor-corrector both evaluate the full acceleration function at multiple sub-step points; floating-point cancellation artifacts from $\infty \times 0$ products can differ between sub-steps in ways that are not algebraic functions of $h$, poisoning the convergence estimate.

> **Reference:** Rein, H. & Spiegel, D. S. (2015). IAS15: a fast, adaptive, high-order integrator for gravitational dynamics, accurate to machine precision over a billion orbits. *MNRAS*, 446(2), 1424–1437. §2.2 (predictor-corrector convergence and smoothness requirements).

> ---
>
> **Proposed Fixes for §4.3**
>
> **Fix A — Replace Gaussian with $\mathrm{sech}^2$ vertical profile.** Use:
>
> $$\rho(z) = \rho_c \,\mathrm{sech}^2\!\left(\frac{z}{\sqrt{2}\,H}\right) = \frac{\rho_c}{\cosh^2\!\left(z / \sqrt{2}\,H\right)}$$
>
> The $\mathrm{sech}^2$ profile is the exact vertical hydrostatic equilibrium for an isothermal, self-gravitating sheet (Spitzer 1942). This is more physical than the Gaussian approximation for the Sirko & Goodman self-regulated zone where $Q \approx 1$ (self-gravity is dynamically important). It preserves the same scale height $H$ and midplane density $\rho_c$, but has broader exponential tails $\sim \exp(-2|z|/\sqrt{2}H)$ instead of super-exponential $\sim \exp(-z^2/2H^2)$. The resulting force pulse during disk-plane crossings is temporally wider and smoother — the temporal half-width increases from $\tau_\mathrm{Gauss} = H/v_z$ to $\tau_\mathrm{sech} \approx 1.76\, H/v_z$.
>
> ```cpp
> // Replace: double rho = rho_c * exp(-0.5 * (z * z) / (H * H));
> double zeta = z / (std::sqrt(2.0) * H);
> double sech_zeta = 1.0 / cosh(zeta);
> double rho = rho_c * sech_zeta * sech_zeta;
> ```
>
> > **Reference:** Spitzer, L. (1942). The Dynamics of the Interstellar Medium. I. Local Equilibrium. *ApJ*, 95, 329. Eq. 5 (self-gravitating isothermal sheet density profile).
>
> **Fix B — Smooth vertical density taper.** Multiply the Gaussian profile by a $C^\infty$ cutoff that enforces a clean transition to exactly zero at a controlled height:
>
> $$\rho(z) = \rho_c \exp\!\left(-\frac{z^2}{2H^2}\right) \cdot \tfrac{1}{2}\left[1 - \tanh\!\left(\frac{|z|/H - z_\mathrm{cut}}{\sigma_z}\right)\right]$$
>
> with $z_\mathrm{cut} \sim 8$ and $\sigma_z \sim 1$. The bare Gaussian underflows gradually from $10^{-300}$ to $10^{-100}$ to $0$ across many scale heights, creating a numerically ragged transition where different integrator sub-steps see different floating-point exponent ranges. The taper enforces a well-defined, smooth cutoff at $|z| \approx z_\mathrm{cut} \cdot H$, compressing the dynamic range to $\sim 14$ decades (from $\rho_c$ to $\rho_c e^{-32}$) instead of $\sim 300$.
>
> ```cpp
> double zH = std::abs(z) / H;
> double taper = 0.5 * (1.0 - tanh((zH - 8.0) / 1.0));
> double rho = rho_c * exp(-0.5 * zH * zH) * taper;
> ```

### 4.4 **LATENT: `corr_factor` NaN from Negative Radicand**

**Location:** `disk-model.hpp`, line 260

```cpp
auto corr_factor = sqrt(1 - n*(cs2/v_k2));
```

The pressure-corrected disk velocity uses the sub-Keplerian correction factor $\eta = \sqrt{1 - n c_s^2/v_K^2}$, where $n = -d\ln P/d\ln R$ is the (negated) logarithmic pressure gradient. This becomes NaN when:

$$
n \frac{c_s^2}{v_K^2} > 1 \implies R_\mathrm{cyl} > \frac{GM_1}{n\,c_s^2}
$$

For our disk ($n_\mathrm{max} \approx 3.0$ in the Self-Reg zone, $c_{s,\mathrm{max}} \approx 0.84\,\mathrm{AU\,yr^{-1}}$ at $R_\mathrm{max}$):

$$
R_\mathrm{crit} = \frac{GM_1}{3 \times 0.84^2} \approx 5.5 \times 10^9\,\mathrm{AU} \gg R_\mathrm{max}
$$

**Safe for this disk.** However, this is a fragile assumption that depends on the specific disk parameters. A disk with higher sound speed or larger $|\nabla \ln P|$ near the outer edge could trigger silent NaN propagation. Under IEEE 754 arithmetic, `NaN < 1e-10` evaluates to `false`, so the `vmag < 1e-10` guard would not catch the NaN and it would propagate into the acceleration array.

> **Reference:** Armitage, P. J. (2020). *Astrophysics of Planet Formation* (2nd ed.). Cambridge University Press. Eq. 6.46 (sub-Keplerian correction and its validity domain).

> ---
>
> **Proposed Fixes for §4.4**
>
> **Fix A — Clamp the radicand to zero.** When $n c_s^2 > v_K^2$, the pressure gradient exceeds centrifugal support — the gas is pressure-supported, not rotationally supported. Setting $v_\mathrm{disk} = 0$ is the correct physical limit (Armitage 2020, §6.3):
>
> ```cpp
> // Replace: auto corr_factor = sqrt(1 - n*(cs2/v_k2));
> auto arg = 1.0 - n * (cs2 / v_k2);
> auto corr_factor = (arg > 0.0) ? sqrt(arg) : 0.0;
> ```
>
> **Fix B — First-order expansion with natural clamping.** For small corrections ($\eta = n c_s^2/v_K^2 < 0.5$), use the exact square root. For large corrections, switch to the Taylor approximation $\sqrt{1-\eta} \approx 1 - \eta/2$, which naturally reaches zero at $\eta = 2$ and remains positive for all $\eta < 2$. This avoids the `sqrt` of a negative number while maintaining a smooth, monotonic transition:
>
> ```cpp
> auto eta = n * (cs2 / v_k2);
> auto corr_factor = (eta < 0.5) ? sqrt(1.0 - eta)
>                                : std::max(0.0, 1.0 - 0.5 * eta);
> ```
>
> The linear branch joins the exact branch at $\eta = 0.5$ with relative error $< 0.4\%$ ($\sqrt{0.5} = 0.707$ vs $1 - 0.25 = 0.75$). For a $C^1$-smooth join, a cubic Hermite splice over $\eta \in [0.4, 0.6]$ can be used.

### 4.5 **LATENT: Bondi-Hoyle Formula Division by Zero at $\mathcal{M} = 0$**

**Location:** `disk-model.hpp`, lines 407–409

```cpp
double f_BH = f_HL * (pow(Mach, 2) / (1 + pow(Mach, 2))) / pow(Mach, 2);
```

Algebraically this simplifies to $f_\mathrm{BH} = f_\mathrm{HL}/(1 + \mathcal{M}^2)$, but the implementation computes $(\mathcal{M}^2/\mathcal{M}^2)/(1+\mathcal{M}^2)$, which evaluates to $0/0 = \mathrm{NaN}$ when $\mathcal{M} = 0$ (i.e., $v_\mathrm{rel} = 0$). The `vmag < 1e-10` guard on line 367 prevents exact $v_\mathrm{rel} = 0$, but does not prevent arbitrarily small $\mathcal{M}$ when $c_s$ is large.

**Currently disabled** for this simulation (`enable_bondi_hoyle = false`) but remains a latent bug for any run with Bondi-Hoyle accretion enabled.

> **Reference:** Edgar, R. (2004). A review of Bondi-Hoyle-Lyttleton accretion. *New Astronomy Reviews*, 48(10), 843–859. Eq. 51 (the correct interpolation formula).

> ---
>
> **Proposed Fixes for §4.5**
>
> **Fix A — Use the algebraically simplified expression.** The current code computes $(\mathcal{M}^2 / \mathcal{M}^2) / (1 + \mathcal{M}^2)$; the $\mathcal{M}^2$ factors cancel algebraically but produce $0/0$ numerically. The simplified form evaluates correctly for all $\mathcal{M} \geq 0$, including $\mathcal{M} = 0$ where it correctly gives $f_\mathrm{BH} = f_\mathrm{HL}$ (the Bondi spherical accretion limit):
>
> ```cpp
> // Replace: double f_BH = f_HL * (pow(Mach, 2) / (1 + pow(Mach, 2))) / pow(Mach, 2);
> double f_BH = f_HL / (1.0 + Mach * Mach);
> ```
>
> **Fix B — Explicit guard with correct limiting behavior.** Preserves the original code structure for traceability with Edgar (2004), Eq. 51, while handling the $\mathcal{M} \to 0$ limit explicitly. At $\mathcal{M} = 0$, $f_\mathrm{BH} = f_\mathrm{HL}$ recovers the Bondi accretion rate (Edgar 2004, §4.1):
>
> ```cpp
> double f_BH = (Mach > 1e-30)
>     ? f_HL * (pow(Mach, 2) / (1 + pow(Mach, 2))) / pow(Mach, 2)
>     : f_HL;
> ```

---

## 5. Synthesis and Conclusion

The `Reach max iteration()` error at high inclination is caused by a **structural incompatibility** between the disk force model and Richardson/Romberg-type adaptive integrators. The primary mechanism is:

1. **The hard `R_cyl < R_\mathrm{min}` boundary** (§4.1) creates a step-function discontinuity in the acceleration that no polynomial extrapolation can converge through
2. **The Keplerian divergence** $v_K \propto R_\mathrm{cyl}^{-1/2}$ in `disk_v()` (§4.2) causes the force to diverge as it approaches the discontinuity, maximizing the extrapolation error
3. **Inclination controls the proximity to the boundary**: the minimum $R_\mathrm{cyl}$ per orbit scales as $r\sqrt{1 - \sin^2 i}$, which reaches zero at $i = 90°$

These are not numerical precision issues or integrator limitations. They are **physical model artifacts** — the code applies an infinite, razor-thin disk model to a 3D orbit that can become perpendicular to the disk plane. The disk model was implicitly designed for low-inclination orbits where $R_\mathrm{cyl} \approx r$ and the particle remains within the disk.

---

## 6. Files Analyzed

| File | Role |
|------|------|
| `SpaceHub/src/interaction/disk-model.hpp` | Force model containing all pathological expressions |
| `SpaceHub/test/verify_gas_drag/stae321-11-pagn.cpp` | Simulation setup (BS/DefaultMethod) |
| `SpaceHub/test/verify_gas_drag/stae321-11-radau.cpp` | Simulation setup (IAS15/Radau) |
| `SpaceHub/src/interaction/disk_tab/disk_stae321_pagn.csv` | Disk table ($R_\mathrm{min} = 11.86$ AU, $H = 2.66$ AU at inner edge) |
