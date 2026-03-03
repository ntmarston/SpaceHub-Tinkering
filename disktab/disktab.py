#Claude was used here to: 
# Write docstrings 
# make my comments intelligible
# optimize performance
# Properly add decorators
from astropy import units as u
from astropy.constants import G, c, k_B, sigma_sb, m_p
from collections import namedtuple
import simple_disk as sd
from matplotlib import pyplot as plt
import numpy as np
from numpy import pi
import pandas as pd
import seaborn as sns

# Pre-computed coefficients for Newton solver (depend on R only, not Tc)
RadiusCoeffs = namedtuple('RadiusCoeffs', ['f', 'Omega', 'A', 'B', 'C', 'D', 'E', 'F'])

def _log_gradient(y, x):
    """Compute -d(ln y)/d(ln x) using central differences."""
    return -np.gradient(np.log(y), np.log(x))


def _smooth_zone_gradients(model, n_gap=100, w=100):
    """
    Smooth grad_P and grad_Sigma at the Standard/Self-Reg zone boundary
    using a C¹-continuous cubic Hermite interpolant.

    Replaces 2*n_gap rows centred on the boundary (n_gap from each zone)
    with a cubic bridge that matches both value and slope at the anchors.
    Original values are saved in grad_P_raw / grad_Sigma_raw columns.

    Parameters
    ----------
    model : pd.DataFrame
        Must contain columns: R/Rg, zone, grad_P, grad_Sigma.
    n_gap : int
        Rows removed on each side of the boundary (default 100).
    w : int
        Half-window for slope estimation outside the gap (default 100).
    """
    import warnings
    from scipy.interpolate import CubicHermiteSpline

    if 'zone' not in model.columns:
        return

    sr_idx = np.where(model['zone'].values == 'Self-Reg')[0]
    if len(sr_idx) == 0:
        return  # single-zone disk — nothing to do

    ibnd  = int(sr_idx[0])
    log_R = np.log10(model['R/Rg'].values)
    N     = len(model)
    i_lo  = ibnd - n_gap
    i_hi  = ibnd + n_gap

    if i_lo <= w or i_hi + w >= N:
        warnings.warn(
            f"Zone boundary smoothing skipped: gap [{i_lo}, {i_hi}) with slope "
            f"window w={w} exceeds data range [0, {N}). Reduce n_gap or w.",
            stacklevel=3,
        )
        return

    for col in ('grad_P', 'grad_Sigma'):
        y = model[col].values.copy()
        model[f'{col}_raw'] = y.copy()

        x0, y0 = log_R[i_lo - 1], y[i_lo - 1]  # last Standard anchor
        x1, y1 = log_R[i_hi],     y[i_hi]        # first Self-Reg anchor
        dy0 = np.mean(np.gradient(y[i_lo - w - 1 : i_lo],     log_R[i_lo - w - 1 : i_lo]))
        dy1 = np.mean(np.gradient(y[i_hi         : i_hi + w + 1], log_R[i_hi : i_hi + w + 1]))

        chs = CubicHermiteSpline([x0, x1], [y0, y1], [dy0, dy1])
        y[i_lo:i_hi] = chs(log_R[i_lo:i_hi])
        model[col] = y

    warnings.warn(
        f"Zone boundary smoothing applied (CubicHermiteSpline, n_gap={n_gap}, w={w}): "
        f"{2 * n_gap} rows [{i_lo}, {i_hi}) around index {ibnd} "
        f"(R/Rg ≈ {model['R/Rg'].iloc[ibnd]:.1f}) replaced for grad_P and grad_Sigma. "
        f"Originals preserved in grad_P_raw / grad_Sigma_raw.",
        stacklevel=3,
    )


class AGNDisk:
    """
    Shakura-Sunyaev AGN disk profile generator using Newton's method.

    Attributes (to be set before calling generate()):
        M: Central black hole mass (astropy Quantity with mass units)
        Mdot: Mass accretion rate (astropy Quantity with mass/time units)
        alpha: Shakura-Sunyaev viscosity parameter (dimensionless)
        mu: Mean molecular weight (dimensionless, default 0.6)
        k0: Electron scattering opacity coefficient (default 0.4 cm^2/g)
        k1: Kramers opacity coefficient (default 6.4e22 cm^5/g^2 K^(7/2))
        Rstar: Inner disk boundary (defaults to 6 * gRad if not set)

    Generated Attributes:
        gRad: Gravitational radius GM/c^2 (computed from M)
        rvals: List of radii used for computation (set by generate())
        model: pandas DataFrame containing the disk profile (SI units)
        model_spacehub: pandas DataFrame in SpaceHub units (set by spacehub_pretab())
        selfreg_threshold: Radius where Q first drops to 1 (astropy Quantity in meters,
                           numerically equal to R/Rg_threshold × gRad). Set by generate()
                           when sirko_goodman=True; None otherwise.
    """

    def __init__(self):
        # Modifiable attributes - set these before calling generate()
        self.M = None
        self.Mdot = None
        self.alpha = None
        self.mu = 0.6  # default mean molecular weight
        self.k0 = 0.4 * u.cm**2 / u.g  # electron scattering opacity
        self.k1 = 6.4e22 * u.cm**5 / u.g**2 * u.K**(7/2)  # Kramers opacity
        self.rvals = None
        self._Rstar_override = None  # internal: user-specified Rstar

        # Computed attributes
        self.gRad = None
        self.Rstar = None
        self.model = None  # DataFrame with standard disk profile
        self.model_ss = None  # DataFrame with standard shakura-sunyaev
        self.model_spacehub = None  # DataFrame in SpaceHub units
        self.selfreg_threshold = None  # Radius where Q=1 (astropy Quantity, meters = X * gRad)

        # Solver settings
        self.max_iter = 10
        self.abs_err_tol = 10
        self.verbose = False

    @property
    def Rstar(self):
        return self._Rstar

    @Rstar.setter
    def Rstar(self, value):
        if value is not None:
            self._Rstar_override = value
        self._Rstar = value

    def generate(self, rvals=None, sirko_goodman=True):
        """
        Generate the disk profile using the NewtonLastPoint method.

        This computes central temperature, density, pressure, sound speed,
        scale height, viscosity, surface density, optical depth, and Toomre Q
        for each radius in rvals.

        The result is stored in self.model as a pandas DataFrame.

        Parameters:
            rvals: List of radii to compute profile at. Can be specified as:
                   - Astropy Quantities with length units (e.g., [10*u.AU, 100*u.AU])
                   - Multiples of gravitational radius using the string "Rg"
                     (computed after M is set), e.g., np.logspace(1, 7, 200) for 10-10^7 Rg
                   If None (default), uses 200 log-spaced points from 10 to 10^7 Rg.
            sirko_goodman: If True (default), also compute the Sirko-Goodman (Q=1)
                          solution for the outer disk where Q < 1, stored in self.model_sg.
        """
        # Validate required attributes
        if self.M is None:
            raise ValueError("M (black hole mass) must be set before calling generate()")
        if self.Mdot is None:
            raise ValueError("Mdot (accretion rate) must be set before calling generate()")
        if self.alpha is None:
            raise ValueError("alpha (viscosity parameter) must be set before calling generate()")

        # Compute gravitational radius
        self.gRad = (G * self.M / c / c).to(u.m)

        # Set Rstar: use override if provided, otherwise default to 6 * gRad
        if self._Rstar_override is not None:
            self._Rstar = self._Rstar_override
        else:
            self._Rstar = 5.99 * self.gRad #prevents divide by zero errors at R=6Rg

        # Set rvals: use parameter if provided, otherwise default
        if rvals is not None:
            # If rvals are dimensionless numbers, interpret as multiples of Rg
            if not hasattr(rvals[0], 'unit'):
                self.rvals = [r * self.gRad for r in rvals]
            else:
                self.rvals = rvals
        elif self.rvals is None:
            self.rvals = [r * self.gRad for r in np.logspace(1, 7, 200)]

        # Pre-compute numerical constants in SI units for performance
        self._G_val = G.si.value
        self._c_val = c.si.value
        self._k_B_val = k_B.si.value
        self._sigma_sb_val = sigma_sb.si.value
        self._m_p_val = m_p.si.value
        self._M_val = self.M.si.value
        self._Mdot_val = self.Mdot.si.value
        self._Rstar_val = self.Rstar.si.value
        self._k0_val = self.k0.si.value
        self._k1_val = self.k1.si.value
        # Pre-computed combined constants
        self._A_val = self._k_B_val / (self.mu * self._m_p_val)
        self._B_val = 4 * self._sigma_sb_val / (3 * self._c_val)

        # 1) Generate full disk solution
        self.model_ss = self._full_disk_solution()
        
        
        # Generate Sirko-Goodman solution unless disabled
        if sirko_goodman:
            idxSR = (self.model_ss["Q"] <= 1).idxmax()
            self.selfreg_threshold = self.model_ss["R/Rg"][idxSR] * self.gRad
            SG_modelStartRow = self.model_ss.loc[idxSR-1]
            sg_rvals = [r * self.gRad for r in self.model_ss["R/Rg"][idxSR:]]
            tc_init = SG_modelStartRow["Tc"] * u.K
            model_sg = self._self_regulating_solution(sg_rvals, init_Tc=tc_init, Qmin=1)
            model_ss = self.model_ss.copy()
            model_ss = model_ss.loc[:idxSR].drop(columns=['tau'])
            model_ss["zone"] = "Standard"
            model_sg["zone"] = "Self-Reg"
            model_combined = pd.concat([model_ss, model_sg.iloc[1:]], ignore_index=True)

            self.model = model_combined
        else:
            self.model = self.model_ss

        # Add logarithmic gradients (zone-separated to avoid boundary artifacts)
        if sirko_goodman and 'zone' in self.model.columns:
            ss_mask = self.model['zone'] == 'Standard'
            sg_mask = self.model['zone'] == 'Self-Reg'

            # Standard zone: numerical gradients (smooth within-zone; one-sided at edge)
            ss_R = self.model.loc[ss_mask, 'R'].values
            for col, gcol in [('Tc','grad_T'), ('Sigma','grad_Sigma'), ('P','grad_P')]:
                self.model.loc[ss_mask, gcol] = _log_gradient(
                    self.model.loc[ss_mask, col].values, ss_R)

            # Self-Reg zone: analytical gradients via implicit differentiation of _fq
            # (see derivations/gradient_derivations.tex for full derivation)
            # Define xi = sqrt(Rstar/R), f = 1 - xi
            # grad_Sigma = 3/2 - xi/(6f)
            # grad_P     = 3   - xi/(3f)
            # grad_T     = -[C3(-3 + xi/(3f)) + 3*C2*Tc] / (4*C1*Tc^4 + C2*Tc)
            #   where C1 = 4*sigma_sb/(3c), C2 = rho*k_B/(mu*m_p), C3 = C1*Tc^4 + C2*Tc
            R_sg  = self.model.loc[sg_mask, 'R'].values
            Tc_sg = self.model.loc[sg_mask, 'Tc'].values
            xi = np.sqrt(self._Rstar_val / R_sg)
            f  = 1.0 - xi

            self.model.loc[sg_mask, 'grad_Sigma'] = 1.5 - xi / (6.0 * f)
            self.model.loc[sg_mask, 'grad_P']     = 3.0 - xi / (3.0 * f)

            Omega_sg = np.sqrt(self._G_val * self._M_val / R_sg**3)
            rho_sg   = Omega_sg**2 / (2.0 * pi * self._G_val)
            C1 = self._B_val                      # 4*sigma_sb/(3c), constant
            C2 = rho_sg * self._A_val             # rho * k_B/(mu*m_p)
            C3 = C1 * Tc_sg**4 + C2 * Tc_sg      # from fq = 0
            xi_3f = xi / (3.0 * f)
            self.model.loc[sg_mask, 'grad_T'] = -(
                C3 * (-3.0 + xi_3f) + 3.0 * C2 * Tc_sg
            ) / (4.0 * C1 * Tc_sg**4 + C2 * Tc_sg)
        else:
            self.model['grad_T'] = _log_gradient(self.model['Tc'].values, self.model['R'].values)
            self.model['grad_Sigma'] = _log_gradient(self.model['Sigma'].values, self.model['R'].values)
            self.model['grad_P'] = _log_gradient(self.model['P'].values, self.model['R'].values)

        # Smooth grad_P and grad_Sigma across the zone boundary.
        # Window is the smaller of 100 or 10% of the total number of radial points.
        _n = int(min(100, 0.1 * len(self.model)))
        _smooth_zone_gradients(self.model, n_gap=_n, w=_n)

        return self.model

    def _compute_coefficients(self, R):
        """
        Compute radius-dependent coefficients for Newton solver.
        These values depend only on R, not on Tc, so they can be computed once per radius.
        """
        f = (1 - np.sqrt(self.Rstar / R))
        Omega = np.sqrt(G * self.M / R**3)

        A = k_B / (self.mu * m_p)  # coefficient of Tc
        B = 4 * sigma_sb / (3 * c)  # coefficient of Tc^4
        C = (self.Mdot * f * Omega**2 / (3 * pi * self.alpha))**(2/3)  # coefficient of Tc^0

        D = self.k0  # coefficient of Tc^0 in opacity
        E = self.k1  # coefficient of Tc^(-7/2) in opacity
        F = ((32 * pi * sigma_sb * Omega * R**3) * (3 * pi * self.alpha) /
             (9 * G * self.M * self.Mdot * f) / (self.Mdot * f * Omega**2)**(1/3))  # coefficient of Tc^4

        return RadiusCoeffs(f=f, Omega=Omega, A=A, B=B, C=C, D=D, E=E, F=F)

    def _f1(self, Tc, R, coeffs=None):
        """
        The implicit equation for central temperature Tc at radius R.
        This function equals zero when Tc is the correct solution.

        Parameters:
            Tc: Central temperature guess
            R: Radius
            coeffs: Optional pre-computed RadiusCoeffs (for performance)
        """
        if coeffs is None:
            coeffs = self._compute_coefficients(R)

        A, B, C, D, E, F = coeffs.A, coeffs.B, coeffs.C, coeffs.D, coeffs.E, coeffs.F

        def rho(Tc):
            a = (C * E / F) * Tc**(-15/2)
            b = (C * D / F) * Tc**(-4) - (A * Tc)
            cee = -B * Tc**(4)

            rho_plus = (-b + np.sqrt(b**2 - 4 * a * cee)) / (2 * a)
            return rho_plus

        rho_val = rho(Tc)
        return rho_val**(2/3) * (D + (E * Tc**(-7/2) * rho_val)) - F * Tc**4

    def _f1prime(self, Tc, R, coeffs=None):
        """
        The derivative of f1 with respect to Tc, for Newton's method.

        Parameters:
            Tc: Central temperature guess
            R: Radius
            coeffs: Optional pre-computed RadiusCoeffs (for performance)
        """
        if coeffs is None:
            coeffs = self._compute_coefficients(R)

        A, B, C, D, E, F = coeffs.A, coeffs.B, coeffs.C, coeffs.D, coeffs.E, coeffs.F

        def rho(Tc):
            a = (C * E / F) * Tc**(-15/2)
            b = (C * D / F) * Tc**(-4) - (A * Tc)
            cee = -B * Tc**(4)

            rho_plus = (-b + np.sqrt(b**2 - 4 * a * cee)) / (2 * a)
            return rho_plus

        def dRho(Tc):
            a1 = (C * E / F)
            b1 = (C * D / F)
            b2 = A
            c1 = -B
            dscrm = 4 * a1 * c1 * Tc**(9/2)

            t1 = 17 * b2 * Tc**5 - 7 * b1

            t2radical = (b2 * Tc**5 - b1)**2 - dscrm
            t2 = (np.sqrt(t2radical)) + (b2 * Tc**5) - b1
            numerator = Tc**(5/2) * ((t1) * (t2) - (46 * a1 * c1 * Tc**(9/2)))
            denominator = 4 * a1 * np.sqrt((b2 * Tc**5 - b1)**2 - dscrm)
            return (numerator / denominator).cgs

        rhotc = rho(Tc)
        drhotc = dRho(Tc)

        top = (E * rhotc * (10 * Tc * drhotc - 21 * rhotc)) + 4 * Tc**(9/2) * (D * drhotc - 6 * F * Tc**3 * (rhotc)**(1/3))
        bottom = 6 * Tc**(9/2) * (rhotc)**(1/3)

        return top / bottom

    def _compute_coefficients_numeric(self, R_val):
        """
        Compute radius-dependent coefficients using pre-computed numerical constants.
        All inputs and outputs are raw floats (SI units).
        """
        f = 1 - np.sqrt(self._Rstar_val / R_val)
        Omega = np.sqrt(self._G_val * self._M_val / R_val**3)

        C = (self._Mdot_val * f * Omega**2 / (3 * pi * self.alpha))**(2/3)
        F = ((32 * pi * self._sigma_sb_val * Omega * R_val**3) * (3 * pi * self.alpha) /
             (9 * self._G_val * self._M_val * self._Mdot_val * f) / (self._Mdot_val * f * Omega**2)**(1/3))

        return (f, Omega, self._A_val, self._B_val, C, self._k0_val, self._k1_val, F)

    def _f1_numeric(self, Tc_val, coeffs_num):
        """
        Numeric version of _f1 using raw floats. All values in SI units.
        """
        A, B, C, D, E, F = coeffs_num[2], coeffs_num[3], coeffs_num[4], coeffs_num[5], coeffs_num[6], coeffs_num[7]

        # Compute rho
        a = (C * E / F) * Tc_val**(-15/2)
        b = (C * D / F) * Tc_val**(-4) - (A * Tc_val)
        cee = -B * Tc_val**(4)
        rho_val = (-b + np.sqrt(b**2 - 4 * a * cee)) / (2 * a)

        return rho_val**(2/3) * (D + (E * Tc_val**(-7/2) * rho_val)) - F * Tc_val**4

    def _f1prime_numeric(self, Tc_val, coeffs_num):
        """
        Numeric version of _f1prime using raw floats. All values in SI units.
        """
        A, B, C, D, E, F = coeffs_num[2], coeffs_num[3], coeffs_num[4], coeffs_num[5], coeffs_num[6], coeffs_num[7]

        # Compute rho
        a = (C * E / F) * Tc_val**(-15/2)
        b = (C * D / F) * Tc_val**(-4) - (A * Tc_val)
        cee = -B * Tc_val**(4)
        rho_val = (-b + np.sqrt(b**2 - 4 * a * cee)) / (2 * a)

        # Compute dRho/dTc
        a1 = (C * E / F)
        b1 = (C * D / F)
        b2 = A
        c1 = -B
        dscrm = 4 * a1 * c1 * Tc_val**(9/2)

        t1 = 17 * b2 * Tc_val**5 - 7 * b1
        t2radical = (b2 * Tc_val**5 - b1)**2 - dscrm
        t2 = np.sqrt(t2radical) + (b2 * Tc_val**5) - b1
        numerator = Tc_val**(5/2) * ((t1) * (t2) - (46 * a1 * c1 * Tc_val**(9/2)))
        denominator = 4 * a1 * np.sqrt((b2 * Tc_val**5 - b1)**2 - dscrm)
        drho_val = numerator / denominator

        top = (E * rho_val * (10 * Tc_val * drho_val - 21 * rho_val)) + 4 * Tc_val**(9/2) * (D * drho_val - 6 * F * Tc_val**3 * rho_val**(1/3))
        bottom = 6 * Tc_val**(9/2) * rho_val**(1/3)

        return top / bottom

    def _newton(self, tc_initial_guess, R_in, return_convergence=False):
        """
        Newton-Raphson solver for finding central temperature Tc at radius R.

        Parameters:
            tc_initial_guess: Initial guess for Tc (astropy Quantity with temperature units)
            R_in: Radius (astropy Quantity with length units)
            return_convergence: If True, also return whether the solver converged

        Returns:
            Tc: Converged central temperature
            converged: (optional) Whether the solver converged
        """
        # Extract numerical values (SI units) for fast inner loop
        tcn_val = tc_initial_guess.si.value
        tc_initial_val = tcn_val  # Save initial guess for comparison
        R_val = R_in.si.value
        abs_err_prev = 1e10
        converged = False
        monotonic_decrease_count = 0  # Track if Tc is decreasing monotonically

        # Pre-compute numeric coefficients once for this radius
        coeffs_num = self._compute_coefficients_numeric(R_val)

        for i in range(self.max_iter):
            # Use numeric methods (no astropy overhead)
            fn = self._f1_numeric(tcn_val, coeffs_num)
            dfn = self._f1prime_numeric(tcn_val, coeffs_num)

            if self.verbose:
                print(f"(n:{i}) Updating Tc_n = {tcn_val} K - {fn/dfn} K")

            # Tc_{n+1}
            tcnext_val = tcn_val - (fn / dfn)
            abs_err = np.abs(tcnext_val - tcn_val)

            # Check for convergence to unphysical values (Tc too low or negative)
            if tcnext_val <= 10:  # Tc below 10 K is unphysical for disk
                raise RuntimeError("Newton solver converging to unphysical Tc <= 10 K")

            # Track monotonic decrease (indicates convergence to wrong root)
            if tcnext_val < tcn_val:
                monotonic_decrease_count += 1
            else:
                monotonic_decrease_count = 0

            # If Tc has decreased monotonically for too many iterations and is far from initial guess
            if monotonic_decrease_count >= 5 and tcnext_val < 0.1 * tc_initial_val:
                raise RuntimeError("Newton solver converging to wrong root (monotonic decrease)")

            if abs_err <= self.abs_err_tol:
                converged = True
                if self.verbose:
                    print(f"Converged after {i} iterations! Abs error: {abs_err:.5e}")
                break

            if abs_err > abs_err_prev:
                raise RuntimeError("Newton solver diverged!")

            abs_err_prev = abs_err
            tcn_val = tcnext_val

        # Convert back to astropy Quantity for return
        result = tcn_val * u.K
        if return_convergence:
            return result, converged
        return result

    def _get_TcPowerLawSol(self, R):
        """
        Get power-law solution for central temperature at radius R.
        Used as initial guess for Newton's method.
        """
        f = (1 - np.sqrt(self.Rstar / R))
        fMdot = f * self.Mdot
        # Use zone 1 (innermost) power law solution
        return sd.T_central(alpha=self.alpha, M=self.M, fMdot=fMdot, r=R, zone=1, scale="AGN")

    def _get_zone_crossings(self):
        """
        Compute the radii where power law zones cross.
        Returns (cross21, cross32) - radii where zone 2 crosses zone 1,
        and zone 3 crosses zone 2.
        """
        # Use a fine grid to find crossings
        test_rvals = [r * self.gRad for r in np.logspace(1, 7, 500)]

        TcPowerLaw = {1: [], 2: [], 3: []}
        for r in test_rvals:
            f = (1 - np.sqrt(self.Rstar / r))
            fMdot = f * self.Mdot
            for z in [1, 2, 3]:
                TcPowerLaw[z].append(
                    sd.T_central(alpha=self.alpha, M=self.M, fMdot=fMdot, r=r, zone=z, scale="AGN")
                )

        tc = lambda z: np.array([t.si.value for t in TcPowerLaw[z]])
        mask = np.array([(r / self.gRad).value > 6 for r in test_rvals])

        def find_crossing(zone_a, zone_b):
            crossings = [test_rvals[i] for i in
                         np.where(mask[:-1] & (np.diff(np.sign(tc(zone_a) - tc(zone_b))) != 0))[0]]
            return crossings[0] if crossings else None

        cross21 = find_crossing(2, 1)
        cross32 = find_crossing(3, 2)

        return cross21, cross32

    def _NewtonLastPoint(self):
        """
        Solve for central temperature at each radius using Newton's method.
        For radii within 100 Rg, uses zone-appropriate power law as initial guess.
        For radii beyond 100 Rg, uses the previous point's solution as the initial guess.
        """
        TcVals = []
        conv_cnt = 0
        divrg_cnt = 0

        # Pre-compute zone crossings for initial guess selection
        inner_boundary_limit = 100 * self.gRad
        cross21, cross32 = self._get_zone_crossings()

        # Set default crossing values if not found
        if cross21 is None:
            cross21 = 6 * self.gRad
        if cross32 is None:
            cross32 = 1e10 * self.gRad

        for i, r in enumerate(self.rvals):
            # Determine if we should use power law (inner region or large gap)
            use_power_law = False

            # Use power law for radii within 100 Rg
            if r < inner_boundary_limit:
                use_power_law = True
            # Also use power law if there's a large gap from previous radius (ratio > 10)
            elif i > 0 and (r / self.rvals[i-1]).decompose().value > 10:
                use_power_law = True
            # Or if this is the first point
            elif i == 0:
                use_power_law = True

            if use_power_law:
                f = (1 - np.sqrt(self.Rstar / r))
                fMdot = f * self.Mdot

                if r < cross21:
                    zone = 1
                elif r <= cross32:
                    zone = 2
                else:
                    zone = 3

                tc_ig = sd.T_central(alpha=self.alpha, M=self.M, fMdot=fMdot, r=r, zone=zone, scale="AGN")

            # For nearby radii beyond 100 Rg: use previous point
            else:
                tc_ig = TcVals[-1]

            try:
                TcOut, conv = self._newton(tc_initial_guess=tc_ig, R_in=r, return_convergence=True)
                if conv:
                    conv_cnt += 1
            except Exception as e:
                # If Newton fails:
                # - If we used power law as initial guess, use that as the result
                # - Otherwise use previous value (for nearby points)
                if use_power_law:
                    TcVals.append(tc_ig)  # Power law is better than nothing for large jumps
                elif len(TcVals) > 0:
                    TcVals.append(TcVals[-1])
                else:
                    TcVals.append(tc_ig)
                divrg_cnt += 1
                continue

            TcVals.append(TcOut)

        if self.verbose:
            print(f"converged: {conv_cnt}")
            print(f"diverged: {divrg_cnt}")
            print(f"other: {len(self.rvals) - (conv_cnt + divrg_cnt)}")

        return TcVals

    def _Rho(self, Tc, R):
        """Compute density from central temperature and radius."""
        f = (1 - np.sqrt(self.Rstar / R))
        Omega = np.sqrt(G * self.M / R**3)

        A = k_B / (self.mu * m_p)
        B = 4 * sigma_sb / (3 * c)
        C = (self.Mdot * f * Omega**2 / (3 * pi * self.alpha))**(2/3)

        D = self.k0
        E = self.k1
        F = ((32 * pi * sigma_sb * Omega * R**3) * (3 * pi * self.alpha) /
             (9 * G * self.M * self.Mdot * f) / (self.Mdot * f * Omega**2)**(1/3))

        a1 = (C * E / F) * Tc**(-15/2)
        b1 = (C * D / F) * Tc**(-4) - (A * Tc)
        c1 = -B * Tc**(4)

        rho_plus = (-b1 + np.sqrt(b1**2 - 4 * a1 * c1)) / (2 * a1)
        return rho_plus

    def _P(self, Tc, rho):
        """Compute pressure from temperature and density."""
        t1 = rho * k_B * Tc / (self.mu * m_p)
        t2 = 4 * sigma_sb * Tc**4 / (3 * c)
        return (t1 + t2).si

    def _cs(self, P, rho):
        """Compute sound speed from pressure and density."""
        return np.sqrt(P / rho).si

    def _H(self, cs, R):
        """Compute scale height from sound speed and radius."""
        Omega = np.sqrt(G * self.M / R**3)
        return (cs / Omega).si

    def _visc(self, cs, H):
        """Compute kinematic viscosity from sound speed and scale height."""
        return (self.alpha * cs * H).si

    def _Sigma(self, visc, R):
        """Compute surface density from viscosity and radius."""
        f = (1 - np.sqrt(self.Rstar / R))
        RHS = self.Mdot * f / (3 * pi)
        return (RHS / visc).cgs

    def _tau(self, Sigma, rho, Tc):
        """Compute optical depth from surface density, density, and temperature."""
        prscr = self.k0 + self.k1 * rho * (Tc**(-7/2))
        return (Sigma * prscr).si

    def _toomreQ(self, Sigma, cs, R):
        """Compute Toomre Q parameter for gravitational stability."""
        Omega = np.sqrt(G * self.M / R**3)
        return (cs * Omega) / (pi * G * Sigma)

    def _full_disk_solution(self):
        """
        Compute full disk solution for all radii using Shakura-Sunyaev model.
        Returns a pandas DataFrame with all disk properties.

        Note on scale height convention:
        - H = cs/Ω (scale height of Gaussian vertical density profile)
        - ρ is the midplane (central) density
        - Σ ≈ √(2π)·ρ·H (from integrating Gaussian profile)
        - The ρ output here is NOT simply Σ/H; it's the central density from the EOS.
        """
        TcVals = self._NewtonLastPoint()

        # Pre-compute constants that don't depend on radius
        A = k_B / (self.mu * m_p)
        B = 4 * sigma_sb / (3 * c)
        D = self.k0
        E = self.k1

        rows = []
        for R, Tc in zip(self.rvals, TcVals):
            # Compute f and Omega once per radius
            f = (1 - np.sqrt(self.Rstar / R))
            Omega = np.sqrt(G * self.M / R**3)

            # Compute remaining radius-dependent coefficients
            C = (self.Mdot * f * Omega**2 / (3 * pi * self.alpha))**(2/3)
            F_coeff = ((32 * pi * sigma_sb * Omega * R**3) * (3 * pi * self.alpha) /
                 (9 * G * self.M * self.Mdot * f) / (self.Mdot * f * Omega**2)**(1/3))

            # Compute rho (inlined from _Rho)
            a1 = (C * E / F_coeff) * Tc**(-15/2)
            b1 = (C * D / F_coeff) * Tc**(-4) - (A * Tc)
            c1 = -B * Tc**(4)
            rhoi = ((-b1 + np.sqrt(b1**2 - 4 * a1 * c1)) / (2 * a1)).si

            # Compute P (inlined from _P)
            t1 = rhoi * k_B * Tc / (self.mu * m_p)
            t2 = 4 * sigma_sb * Tc**4 / (3 * c)
            Pi = (t1 + t2).si

            # Compute cs (inlined from _cs)
            csi = np.sqrt(Pi / rhoi).si

            # Compute H using pre-computed Omega (inlined from _H)
            Hi = (csi / Omega).si

            # Compute visc (inlined from _visc)
            visci = (self.alpha * csi * Hi).si

            # Compute Sigma using pre-computed f (inlined from _Sigma)
            RHS = self.Mdot * f / (3 * pi)
            Sigmai = (RHS / visci).si

            # Compute tau (inlined from _tau)
            prscr = self.k0 + self.k1 * rhoi * (Tc**(-7/2))
            taui = (Sigmai * prscr).si

            # Compute Toomre Q using pre-computed Omega (inlined from _toomreQ)
            toomreQi = ((csi * Omega) / (pi * G * Sigmai)).si

            row = {
                "R": R.si.value,
                "R/Rg": (R / self.gRad).value,
                "Tc": Tc.si.value,
                "rho": rhoi.value,
                "P": Pi.value,
                "cs": csi.value,
                "H": Hi.value,
                "visc": visci.value,
                "Sigma": Sigmai.value,
                "tau": taui.value,
                "Q": toomreQi.value
            }
            rows.append(row)

        return pd.DataFrame(rows)

    # ======================== Sirko-Goodman (Q=1) Methods ========================

    def _fq(self, Tc, R, Qmin=1):
        """
        The implicit equation for central temperature Tc at radius R
        for the Sirko-Goodman (Q=1) self-regulating disk.
        This function equals zero when Tc is the correct solution.
        """
        Omega = np.sqrt(G * self.M / R**3)
        f = (1 - np.sqrt(self.Rstar / R))
        rho = Omega**2 / (2 * pi * G * Qmin)
        C1 = 4 * sigma_sb / 3 / c
        C2 = rho * k_B / (self.mu * m_p)
        C3 = rho * ((self.Mdot * f * Omega**2) / (3 * 2 * pi * self.alpha * rho))**(2/3)
        return (C1 * Tc**4) + C2 * Tc - C3

    def _fqprime(self, Tc, R, Qmin=1):
        """
        The derivative of fq with respect to Tc, for Newton's method.
        """
        Omega = np.sqrt(G * self.M / R**3)
        rho = Omega**2 / (2 * pi * G * Qmin)
        C1 = 4 * sigma_sb / 3 / c
        C2 = rho * k_B / (self.mu * m_p)
        return 4 * C1 * Tc**3 + C2

    def _newton_q(self, tc_initial_guess, R_in, Qmin=1, return_convergence=False):
        """
        Newton-Raphson solver for finding central temperature Tc at radius R
        using the Sirko-Goodman (Q=1) prescription.

        Parameters:
            tc_initial_guess: Initial guess for Tc (astropy Quantity with temperature units)
            R_in: Radius (astropy Quantity with length units)
            Qmin: Target Toomre Q value (default 1)
            return_convergence: If True, also return whether the solver converged

        Returns:
            Tc: Converged central temperature
            converged: (optional) Whether the solver converged
        """
        tcn = tc_initial_guess
        tc_iters = [tc_initial_guess]
        abs_err_prev = 1e10
        converged = False

        for i in range(self.max_iter):
            fn = self._fq(Tc=tcn, R=R_in, Qmin=Qmin)
            dfn = self._fqprime(Tc=tcn, R=R_in, Qmin=Qmin)

            if self.verbose:
                print(f"(n:{i}) [SG] Updating Tc_n = {tcn.si} - {(fn/dfn).si}")

            # Tc_{n+1}
            tcnext = tcn - (fn / dfn)
            tc_iters.append(tcnext.si)
            abs_err = np.abs(tcnext - tcn).si.value

            if abs_err <= self.abs_err_tol:
                converged = True
                if self.verbose:
                    print(f"[SG] Converged after {i} iterations! Abs error: {abs_err:.5e}")
                break

            if abs_err > abs_err_prev:
                raise RuntimeError("Newton solver (Sirko-Goodman) diverged!")

            abs_err_prev = abs_err
            tcn = tcnext

        if return_convergence:
            return tc_iters[-1], converged
        return tc_iters[-1]

    def _NewtonLastPoint_SirkoGoodman(self, r_list, init_Tc=None):
        """
        Solve for Tc using the Sirko-Goodman (Q=1) prescription.
        Uses the previous point's solution as the initial guess for the next.

        Parameters:
            r_list: List of radii to solve at
            init_Tc: Initial temperature guess for first point (default 10^4.5 K)

        Returns:
            rvals: List of radii
            TcVals: List of central temperatures
        """
        if init_Tc is None:
            init_Tc = 10**(4.5) * u.K

        rvals = r_list
        TcVals = []
        conv_cnt = 0
        divrg_cnt = 0

        for r in rvals:
            if len(TcVals) < 1:
                tc_ig = init_Tc
            else:
                tc_ig = TcVals[-1]

            try:
                TcOut, conv = self._newton_q(
                    tc_initial_guess=tc_ig, R_in=r,
                    return_convergence=True
                )
                if conv:
                    conv_cnt += 1
            except Exception:
                TcVals.append(TcVals[-1] if TcVals else tc_ig)
                divrg_cnt += 1
                continue

            TcVals.append(TcOut)

        if self.verbose:
            print(f"[SG] converged: {conv_cnt}")
            print(f"[SG] diverged: {divrg_cnt}")
            print(f"[SG] other: {len(rvals) - (conv_cnt + divrg_cnt)}")

        return rvals, TcVals

    def _self_regulating_solution(self, r_vals, init_Tc=None, Qmin=1):
        """
        Compute disk solution for the self-regulating (Q=Qmin) outer disk.
        Analogous to _full_disk_solution() but uses the Sirko-Goodman prescription.

        Parameters:
            r_vals: List of radii to compute solution at
            init_Tc: Initial temperature guess for Newton solver
            Qmin: Target Toomre Q value (default 1)

        Returns:
            DataFrame with columns: R, R/Rg, Tc, rho, P, cs, H, visc, Sigma, Q

        Note on scale height convention (Sirko & Goodman):
        - This regime assumes a UNIFORM vertical density (slab model): ρ(z) = ρ_c
        - Σ = 2·ρ·H (factor of 2 from uniform profile, not Gaussian)
        - H = cs/Ω is still the scale height parameter
        - The factor of 2 in the cs calculation (line 765: 3 * 2 * pi * self.alpha * rho)
          comes from this alternative Σ definition. The ρ here is uniform density,
          NOT the midplane density of a Gaussian profile.
        - Therefore, H is consistent between Standard and Self-Reg regimes (both H = cs/Ω),
          but ρ has different physical meanings and creates ~1.25× discontinuity at Q=1 boundary.
        """
        rvals, TcVals = self._NewtonLastPoint_SirkoGoodman(r_list=r_vals, init_Tc=init_Tc)

        rows = []
        for R, Tc in zip(rvals, TcVals):
            nr = {}
            nr["R"] = R.si.value
            nr["R/Rg"] = (R / self.gRad).decompose().value
            nr["Tc"] = Tc.si.value

            f = (1 - np.sqrt(self.Rstar / R))
            Omega = np.sqrt(G * self.M / R**3)
            #Note extra factor of two here because S&G assume Surface Density = 2 * rho * h (Uniform slab) instead of a Gaussian density distribution 
            rho = Omega**2 / (2 * pi * G * Qmin)
            # Note: The factor of 2 in the denominator comes from Σ = 2ρH (Sirko & Goodman's uniform slab model)
            # vs Σ ≈ √(2π)ρH (Gaussian), so this ρ is uniform density, not Gaussian midplane density
            cs = ((self.Mdot * f * Omega**2) / (3 * 2 * pi * self.alpha * rho))**(1/3) 
            P = cs**2 * rho
            H = cs / Omega
            visc = self.alpha * cs * H
            Sigma = self.Mdot * f / (3 * pi * visc)
            Q = cs * Omega / (pi * G * Sigma)

            nr["rho"] = rho.si.value
            nr["P"] = P.si.value
            nr["cs"] = cs.si.value
            nr["H"] = H.si.value
            nr["visc"] = visc.si.value
            nr["Sigma"] = Sigma.si.value
            nr["Q"] = Q.decompose().value

            rows.append(nr)

        return pd.DataFrame(rows)

    # ==============================================================================
    #                            OUTPUT TOOLS
    # ==============================================================================
    def layout(self, include_powerlaw=True, include_standard=False, figsize=(10, 13)):
        """
        Create a layout of the disk profile with parameters displayed at the top.

        Parameters:
            include_powerlaw: If True, overlay power-law solutions for comparison
            include_standard: If True, also plot the standard Shakura-Sunyaev model
                             (model_ss) alongside the Sirko-Goodman model for comparison.
            figsize: Figure size tuple (width, height)

        Returns:
            fig, axs: matplotlib figure and axes array
        """
        from matplotlib.gridspec import GridSpec

        if self.model is None:
            raise ValueError("Must call generate() before layout()")

        plot_data = self.model
        plot_label = "Sirko-Goodman"
        cols = ["Tc", "rho", "P", "cs", "H", "visc", "Sigma", "Q"]
        ylabels = [
            r"Central Temperature $T_c$ (K)",
            r"Density $\rho$ (kg/m$^3$)",
            r"Pressure $P$ (Pa)",
            r"Sound Speed $c_s$ (m/s)",
            r"Scale Height $H$ (m)",
            r"Viscosity $\nu$ (m$^2$/s)",
            r"Surface Density $\Sigma$ (kg/m$^2$)",
            r"Toomre Q Parameter"
        ]

        # Set up color palette
        palette = sns.color_palette("icefire", n_colors=10)

        # Create figure with GridSpec: 4 rows for plots only
        fig = plt.figure(figsize=figsize)
        gs = GridSpec(4, 2, figure=fig, hspace=0.3, top=0.93, bottom=0.05)

        # Create the 4x2 grid of plot axes
        axs = []
        for row in range(4):
            for col in range(2):
                ax = fig.add_subplot(gs[row, col])
                axs.append(ax)

        # Format parameter values for display
        M_str = f"{self.M.to(u.Msun).value:.2e} M$_\\odot$"
        Mdot_str = f"{self.Mdot.to(u.Msun/u.yr).value:.2f} M$_\\odot$/yr"
        alpha_str = f"{self.alpha:.1e}"
        mu_str = f"{self.mu:.2f}"
        Rstar_str = f"{(self.Rstar / self.gRad).decompose().value:.1f} R$_g$"
        gRad_str = f"{self.gRad.to(u.AU).value:.2e} m"

        # Build parameter text
        param_text = (
            f"M = {M_str}    "
            f"$\\dot{{M}}$ = {Mdot_str}    "
            f"$\\alpha$ = {alpha_str}    "
            f"$\\mu$ = {mu_str}    "
            f"R$_*$ = {Rstar_str}    "
            f"R$_g$ = {gRad_str}"
        )

        # Add centered text box with parameters above the plots
        bbox_props = dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.8, edgecolor='black')
        fig.text(
            0.5, 0.96, param_text,
            fontsize=10,
            verticalalignment='top',
            horizontalalignment='center',
            bbox=bbox_props
        )

        # Get power-law solutions if requested
        if include_powerlaw:
            df_powerlaw = self._get_powerlaw_df()

        for i, col_name in enumerate(cols):
            # Only add labels on the first plot for the legend
            sns.scatterplot(
                ax=axs[i], data=plot_data, x="R/Rg", y=col_name,
                s=5, color=palette[8], label=plot_label if i == 0 else None
            )

            if include_standard and self.model_ss is not None and col_name in self.model_ss.columns:
                sns.scatterplot(
                    ax=axs[i], data=self.model_ss, x="R/Rg", y=col_name,
                    s=5, color=palette[2], label="Shakura-Sunyaev" if i == 0 else None, alpha=0.7
                )

            if include_powerlaw and col_name in df_powerlaw.columns:
                sns.lineplot(
                    ax=axs[i], data=df_powerlaw, x="R/Rg", y=col_name, hue="zone",
                    palette=[palette[1], palette[2], palette[3]],
                    alpha=0.7, legend=(i == 0), linestyle="--", lw=1
                )

            axs[i].set_xscale("log")
            axs[i].set_yscale("log")
            axs[i].set_ylabel(ylabels[i], labelpad=-2)

        axs[0].legend(loc="upper right", fontsize=8)

        # Apply tight_layout with space for the parameter box at the top
        plt.tight_layout(rect=[0, 0, 1, 0.93])

        return fig, axs

    def _get_powerlaw_df(self):
        """
        Compute power-law solutions for all variables at each radius.
        Used for comparison plotting.
        """
        rows = []
        for zone in [1, 2, 3]:
            for r in self.rvals:
                f = (1 - np.sqrt(self.Rstar / r))
                fMdot = f * self.Mdot

                Tc = sd.T_central(alpha=self.alpha, M=self.M, fMdot=fMdot, r=r, zone=zone, scale="AGN")
                rhov = sd.rho(alpha=self.alpha, M=self.M, fMdot=fMdot, r=r, zone=zone, scale="AGN")
                P = sd.pressure(alpha=self.alpha, M=self.M, fMdot=fMdot, r=r, zone=zone, scale="AGN")
                H = sd.h(alpha=self.alpha, M=self.M, fMdot=fMdot, r=r, zone=zone, scale="AGN")
                Sig = sd.Sigma(alpha=self.alpha, M=self.M, fMdot=fMdot, r=r, zone=zone, scale="AGN")
                tau = sd.tau_es(alpha=self.alpha, M=self.M, fMdot=fMdot, r=r, zone=zone, scale="AGN")
                cs = np.sqrt(P / rhov).si
                visc = (self.alpha * cs * H).si

                rows.append({
                    "R/Rg": (r / self.gRad).value,
                    "Tc": Tc.si.value,
                    "rho": rhov.si.value,
                    "P": P.si.value,
                    "cs": cs.value,
                    "H": H.si.value,
                    "visc": visc.value,
                    "Sigma": Sig.si.value,
                    "tau": tau.si.value,
                    "zone": zone
                })

        return pd.DataFrame(rows)


    def spacehub_pretab(self, rvals=None):
        """
        Convert the disk model to SpaceHub's native unit system.

        SpaceHub unit system:
            - Mass: Solar masses (M☉)
            - Length: AU
            - Time: year/(2π), so that 1 year = 2π and G = 1
            - Temperature: Kelvin (unchanged)

        Parameters:
            rvals: Optional list of radii. If provided and model doesn't exist,
                   calls generate(rvals) first. Can be dimensionless (multiples of Rg)
                   or astropy Quantities with length units.
                   If None and model exists, uses existing model.
                   If None and model doesn't exist, uses default radii (10-10^7 Rg).

        Returns:
            DataFrame with columns in SpaceHub units:
                - R: AU
                - R/Rg: dimensionless (unchanged)
                - Tc: K (unchanged)
                - rho: M☉/AU³
                - P: M☉/(AU·T²) where T = yr/(2π)
                - cs: AU/T
                - H: AU
                - visc: AU²/T
                - Sigma: M☉/AU²
                - Q: dimensionless (unchanged)
                - grad_P: -dlnP/dlnr, dimensionless (unchanged)
                - zone: string (unchanged)

        The result is also stored in self.model_spacehub.
        """
        # Validate radial resolution before generating
        n_points = len(rvals) if rvals is not None else (len(self.model) if self.model is not None else 0)
        if n_points < 500:
            raise ValueError(
                f"spacehub_pretab requires at least 500 radial points, got {n_points}. "
                "Increase the number of radii passed via rvals."
            )

        # Generate model if needed
        if self.model is None or rvals is not None:
            self.generate(rvals=rvals)

        # Unit conversion constants (SI values)
        AU_m = 1.495978707e11       # meters per AU
        Msun_kg = 1.98847e30        # kg per solar mass
        year_s = 365.25636042 * 24 * 3600  # seconds per year
        T_unit = year_s / (2 * np.pi)      # SpaceHub time unit in seconds

        # Derived conversion factors (SI → SpaceHub)
        # To convert X_SI to X_SH: X_SH = X_SI * conv_factor
        conv_length = 1 / AU_m                          # m → AU
        conv_density = AU_m**3 / Msun_kg                # kg/m³ → M☉/AU³
        conv_pressure = AU_m * T_unit**2 / Msun_kg      # Pa → M☉/(AU·T²)
        conv_velocity = T_unit / AU_m                   # m/s → AU/T
        conv_viscosity = T_unit / AU_m**2               # m²/s → AU²/T
        conv_surface_density = AU_m**2 / Msun_kg        # kg/m² → M☉/AU²

        # Create converted DataFrame
        df = self.model.copy()

        # Apply conversions
        df["R"] = df["R"] * conv_length
        # R/Rg is dimensionless, unchanged
        # Tc is in Kelvin, unchanged
        df["rho"] = df["rho"] * conv_density
        df["P"] = df["P"] * conv_pressure
        df["cs"] = df["cs"] * conv_velocity
        df["H"] = df["H"] * conv_length
        df["visc"] = df["visc"] * conv_viscosity
        df["Sigma"] = df["Sigma"] * conv_surface_density
        # Q, grad_T, grad_Sigma, grad_P are dimensionless, unchanged
        # zone is a string label, unchanged

        # Reorder: 13 numeric base columns → zone → _raw columns (if any).
        # C++ reads by position and stops at grad_P (col 13), so _raw cols
        # placed after zone are invisible to disk-model.hpp.
        raw_cols  = [c for c in df.columns if c.endswith('_raw')]
        base_cols = [c for c in df.columns if c != 'zone' and not c.endswith('_raw')]
        if 'zone' in df.columns:
            df = df[base_cols + ['zone'] + raw_cols]
        elif raw_cols:
            df = df[base_cols + raw_cols]

        self.model_spacehub = df
        return df


if __name__ == "__main__":
    import argparse
    from astropy.constants import sigma_T

    parser = argparse.ArgumentParser(
        description="Generate AGN disk profile CSV in SpaceHub units."
    )
    parser.add_argument("--M_msun", type=float, required=True,
                        help="Central BH mass in solar masses")
    parser.add_argument("--alpha", type=float, required=True,
                        help="Shakura-Sunyaev viscosity parameter")
    parser.add_argument("--mdot_edd_frac", type=float, required=True,
                        help="Accretion rate as fraction of Eddington")
    parser.add_argument("--epsilon", type=float, default=0.1,
                        help="Radiative efficiency (default: 0.1)")
    parser.add_argument("--r_min_rg", type=float, default=10,
                        help="Min radius in Rg (default: 10)")
    parser.add_argument("--r_max_rg", type=float, default=1e9,
                        help="Max radius in Rg (default: 1e9)")
    parser.add_argument("--n_points", type=int, default=500,
                        help="Number of radial points (default: 500)")
    parser.add_argument("--output", type=str, required=True,
                        help="Output CSV file path")

    args = parser.parse_args()

    disk = AGNDisk()
    disk.M = args.M_msun * u.Msun
    disk.alpha = args.alpha

    # Compute Eddington accretion rate
    kappa = sigma_T / m_p
    L_Edd = 4 * np.pi * G * disk.M * c / kappa
    Mdot_Edd = (L_Edd / (args.epsilon * c**2)).to(u.kg / u.s)
    disk.Mdot = args.mdot_edd_frac * Mdot_Edd

    # Generate and save in SpaceHub units
    rvals = np.logspace(np.log10(args.r_min_rg), np.log10(args.r_max_rg), args.n_points)
    df = disk.spacehub_pretab(rvals=rvals)
    df.to_csv(args.output, index=False)

    print(f"Disk profile written to {args.output} ({len(df)} rows)")