from astropy import units as u
from astropy.constants import G, c
import numpy as np
from numpy import pi
from math import gamma as Gamma
from scipy.integrate import quad
from datetime import datetime, timezone
import itertools
import argparse
import h5py

def f(M, r, v, mf, gamma_bw, r_inf):
    """Distribution function in terms of r, v. Merritt 3.49"""
    phi0 = G*M / r_inf
    Energy = 0.5*v**2 - (G*M)/r
    Af = (3-gamma_bw)/8 * np.sqrt(2)/np.sqrt(pi**5)  * (Gamma(gamma_bw+1) / Gamma(gamma_bw - 1/2)) * (M/mf) * (phi0**(3/2) / (G*M)**3)
    return Af * (np.abs(Energy)/phi0)**(gamma_bw-3/2)

def _standard_integrand_EF(n, M, r, vf, v, mf, gamma_bw, r_inf):
    """Private helper method for functions E and F. f is evaluated at field-star speed vf."""
    return (vf/v)**n * f(M, r, vf, mf, gamma_bw, r_inf)

def E(n, M, r, v, mf, gamma_bw, r_inf):
    """ 'Standard' form integral used in diffusion coefficient calculation. Merritt 5.56a"""
    vu = v.unit
    v_esc = np.sqrt(2*G*M/r).to(vu).value
    gu = _standard_integrand_EF(n, M, r, v, v, mf, gamma_bw, r_inf).decompose().unit
    result, _ = quad(lambda vf: _standard_integrand_EF(n, M, r, vf*vu, v, mf, gamma_bw, r_inf).to_value(gu), v.to_value(vu), v_esc)
    return result * gu * vu

def F(n, M, r, v, mf, gamma_bw, r_inf):
    """ 'Standard' form integral used in diffusion coefficient calculation. Merritt 5.56b"""
    vu = v.unit
    gu = _standard_integrand_EF(n, M, r, v, v, mf, gamma_bw, r_inf).decompose().unit
    result, _ = quad(lambda vf: _standard_integrand_EF(n, M, r, vf*vu, v, mf, gamma_bw, r_inf).to_value(gu), 0, v.to_value(vu))
    return result * gu * vu

def sigma2(M, r, gamma_bw):
    """1D velocity dispersion of stars in the NSC"""
    return G * M / (r * (1+gamma_bw))

def lnLambda(M, r, m, mf, gamma_bw):
    """Coulomb logaritm ?Maybe Chandrasekhar?"""
    b_max = r/gamma_bw #n(r)/(dn/dr) using equation 3.48 in Merritt
    b_min = G * (m+mf) / sigma2(M, r, gamma_bw) #! possibly missing a factor of 1/3 to match the Chandrasekhar convention?
    return np.log(b_max/b_min)

#Diffusion Coefficients — source?
def DC_1par(M, r, v, m, mf, gamma_bw, r_inf):
    """1st order parallel diffusion coefficient. Corresponds to dynamical friction."""
    return -16 * pi**2 * G**2 * (mf + m) * mf * lnLambda(M, r, m, mf, gamma_bw) * F(2, M, r, v, mf, gamma_bw, r_inf)

def DC_2par(M, r, v, m, mf, gamma_bw, r_inf):
    """2nd order parallel diffusion coefficient."""
    const = 32 * pi**2 * G**2 * mf**2 / 3
    return const * lnLambda(M, r, m, mf, gamma_bw) * (F(4, M, r, v, mf, gamma_bw, r_inf) + E(1, M, r, v, mf, gamma_bw, r_inf)) * v

def DC_2perp(M, r, v, m, mf, gamma_bw, r_inf):
    """2nd order perpendicular diffusion coefficient."""
    const = 32 * pi**2 * G**2 * mf**2 / 3
    return const * lnLambda(M, r, m, mf, gamma_bw) * (3*F(2, M, r, v, mf, gamma_bw, r_inf) - F(4, M, r, v, mf, gamma_bw, r_inf) + 2*E(1, M, r, v, mf, gamma_bw, r_inf)) * v


#========RELAXATION TIME==============
#! notes say to fix $v=\sqrt{\frac{GM_{\bullet}}{r}}$ for validation plots
#! I added a negative sign to the Dynamical Friction relaxation time because otherwise it is negative, and that did not make sense to me.
#! more common form would add a factor of 1/3 somewhere, but we are choosing to ignore that


def t_rel_DF(M, r, v, m, mf, gamma_bw, r_inf):
    """first order (dynamical friction) relaxation time"""
    dc = DC_1par(M, r, v, m, mf, gamma_bw, r_inf)
    sigma = np.sqrt(sigma2(M, r, gamma_bw))
    return -sigma/dc

def t_rel_2par(M, r, v, m, mf, gamma_bw, r_inf):
    """parallel second order relaxation time"""
    dc = DC_2par(M, r, v, m, mf, gamma_bw, r_inf)
    return sigma2(M, r, gamma_bw)/dc

def t_rel_2perp(M, r, v, m, mf, gamma_bw, r_inf):
    """perpdendicular second order relaxation time"""
    dc = DC_2perp(M, r, v, m, mf, gamma_bw, r_inf)
    return sigma2(M, r, gamma_bw)/dc



def _progress(slab_i, n_slabs, frac, width=40):
    """In-place two-line status: current slab + a bar tied to frac in [0,1] (ANSI, refreshes rather than scrolls)."""
    bar = "#" * int(round(width * frac)) + "-" * (width - int(round(width * frac)))
    print(f"\rComputing slab [{slab_i}/{n_slabs}].\033[K\n\rProgress: [{bar}] {frac*100:5.1f}%\033[K\033[F", end="", flush=True)


def compute(out_path, rv_table_len, M_exps=[5,6,7,8], m_star_vals = [0.3], gamma_vals=[7/4], r_inf_scaled_vals = [0.5, 1], Rmin = 10, Rmax = 1e5):
    """Pre-tabulate the Merritt standard integrals F2, F4, E1 over an (r, v) grid for every (M, m_star, gamma, r_inf) combo and dump one slab each into a single HDF5 file (SpaceHub natural units).

    out_path           where to write the .h5
    rv_table_len       grid size; tables are square (rv_table_len x rv_table_len)
    M_exps             central BH masses as log10(M / Msun)
    m_star_vals        field-star masses [Msun] (plain floats, not Quantities)
    gamma_vals         Bahcall-Wolf density slopes [dimensionless]
    r_inf_scaled_vals  influence-radius scale factors (r_inf = val * pc * sqrt(M / 1e6 Msun))
    Rmin, Rmax         inner/outer radial range [units of Rg]
    """
    try:
        m_star_vals = [mf * u.Msun for mf in m_star_vals]

        # --- SpaceHub natural units: G=1, length=AU, mass=Msun, time=yr/(2pi) ---
        T_unit = (1 * u.yr) / (2 * np.pi)   # SpaceHub time unit
        F_nat = T_unit**2 / u.AU**5         # F2/F4/E1 decompose to time^2 / length^5
        combos = list(itertools.product(M_exps, m_star_vals, gamma_vals, r_inf_scaled_vals))
        v_norm = np.linspace(0.01, 1.0, rv_table_len)   # one normalized x=v/v_esc axis for all slabs

        with h5py.File(out_path, "w") as h5:
            # ---- root metadata ----
            h5.attrs["description"] = "Merritt standard integrals F_n, E_n for a Bahcall-Wolf NSC"
            h5.attrs["created"] = datetime.now(timezone.utc).isoformat()
            h5.attrs["generated_by"] = "compute_DCs.py(v0.1)"
            h5.attrs["unit_system"] = "SpaceHub natural: G=1, length=AU, mass=Msun, time=yr/(2pi)"
            h5.attrs["length_unit"] = "AU"
            h5.attrs["mass_unit"] = "Msun"
            h5.attrs["time_unit"] = "yr/(2pi)"
            h5.attrs["F_units"] = "(yr/(2pi))^2 / AU^5"
            h5.attrs["v_norm_def"] = "x = v / v_esc, v_esc = sqrt(2 G M / r) with G=1"
            h5.attrs["Rmin_Rg"] = Rmin
            h5.attrs["Rmax_Rg"] = Rmax
            h5.attrs["rv_table_len"] = rv_table_len
            h5.attrs["n_slabs"] = len(combos)

            for slab_i, (M_exp, m_star, gamma, r_inf_scaled) in enumerate(combos):
                # runs once per combination of (M_exp, m_star, gamma_bw, r_inf_scaled)
                M = 10**M_exp * u.Msun
                r_inf = r_inf_scaled * u.pc * (M / 1e6 / u.Msun) ** (0.5)
                n_inf = (3 - gamma) * M / (4 * pi * m_star * r_inf**3)
                Rg = 2 * G * M / c / c
                # r_vals log-spaced between Rmin and Rmax (in Rg), stored in AU
                r_vals = np.geomspace(Rmin, Rmax, rv_table_len) * Rg.to(u.pc)
                r_AU = r_vals.to_value(u.AU)

                F2arr = np.empty((rv_table_len, rv_table_len)) #creates empty/unitialized 2d np.ndarray w/ shape rv_table_len x rv_table_len
                F4arr = np.empty((rv_table_len, rv_table_len))
                E1arr = np.empty((rv_table_len, rv_table_len))

                for i, r in enumerate(r_vals):
                    v_esc = np.sqrt(2 * G * M / r).to(u.km / u.s)
                    for j, x in enumerate(v_norm):
                        v = x * v_esc   # physical velocity at this (r, x)
                        F2arr[i, j] = F(2, M, r, v, m_star, gamma, r_inf).to_value(F_nat)
                        F4arr[i, j] = F(4, M, r, v, m_star, gamma, r_inf).to_value(F_nat)
                        E1arr[i, j] = E(1, M, r, v, m_star, gamma, r_inf).to_value(F_nat)
                    _progress(slab_i + 1, len(combos), (i + 1) / rv_table_len)

                # ---- write this slab ----
                g = h5.create_group(f"slab_{slab_i:04d}")
                g.attrs["M_exp"] = M_exp
                g.attrs["M_Msun"] = M.to_value(u.Msun)
                g.attrs["m_star_Msun"] = m_star.to_value(u.Msun)
                g.attrs["gamma"] = gamma
                g.attrs["r_inf_pc"] = r_inf.to_value(u.pc)
                g.attrs["r_inf_scaled"] = r_inf_scaled
                g.attrs["n_inf_per_pc3"] = n_inf.to_value(u.pc**-3)
                g.attrs["n_r"] = rv_table_len
                g.attrs["n_v"] = rv_table_len

                d = g.create_dataset("r", data=r_AU);      d.attrs["units"] = "AU"
                d = g.create_dataset("v_norm", data=v_norm); d.attrs["units"] = "dimensionless (v / v_esc)"
                for name, arr in (("F2", F2arr), ("F4", F4arr), ("E1", E1arr)):
                    d = g.create_dataset(name, data=arr)
                    d.attrs["units"] = "(yr/(2pi))^2 / AU^5"


        print("\n")
        print(f"Wrote {len(combos)} slabs to {out_path}")
    except Exception as e:
        raise Exception(f"Compute failed, exception occurred: {e.with_traceback(None)}")
        

    return 0


def _parse_args():
    p = argparse.ArgumentParser(description="Pre-tabulate Merritt F2/F4/E1 slabs to a single HDF5 file.")
    p.add_argument("--out-path", default="dc_tables.h5", help="output HDF5 file")
    p.add_argument("--rv-table-len", type=int, default=25, help="grid size; tables are square")
    p.add_argument("--M-exps", type=int, nargs="+", default=[5, 6, 7, 8], help="central BH masses as log10(M/Msun)")
    p.add_argument("--m-star-vals", type=float, nargs="+", default=[0.3], help="field-star masses [Msun]")
    p.add_argument("--gamma-vals", type=float, nargs="+", default=[7/4], help="Bahcall-Wolf density slopes")
    p.add_argument("--r-inf-scaled-vals", type=float, nargs="+", default=[0.5, 1.0], help="influence-radius scale factors")
    p.add_argument("--Rmin", type=float, default=10, help="inner radius [Rg]")
    p.add_argument("--Rmax", type=float, default=1e5, help="outer radius [Rg]")
    return p.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    compute(
        out_path=args.out_path,
        rv_table_len=args.rv_table_len,
        M_exps=args.M_exps,
        m_star_vals=args.m_star_vals,
        gamma_vals=args.gamma_vals,
        r_inf_scaled_vals=args.r_inf_scaled_vals,
        Rmin=args.Rmin,
        Rmax=args.Rmax,
    )


"""
====EXAMPLE COMMAND TO RUN FROM TERMINAL=====
cd disktab/diffusion_coeffs

python compute_DCs.py \
  --out-path test_script.h5 \
  --rv-table-len 11 \
  --M-exps 5 6 7 \
  --m-star-vals 0.3 1.0 \
  --gamma-vals 1.75 \
  --r-inf-scaled-vals 1.0 \
  --Rmin 10 \
  --Rmax 1e5


"""