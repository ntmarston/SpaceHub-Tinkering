
from astropy import units as u
from astropy.constants import G, c

def calc(scaling, alpha, M, fM, r, pindices):
    """Helper function to calculate disk properties"""
    return scaling * (alpha ** pindices[0]) * (M**pindices[1]) * (fM ** pindices[2]) * (r ** pindices[3])

# ============================================================================
# DISK HALF THICKNESS
# ============================================================================
def h(alpha, M, fMdot, r, zone, scale="AGN"):
    """Disk half thickness"""
    scale = scale.lower().strip()
    if scale == "agn":
        scaling_zones = [1.59e12 * u.m, 4.45e8*u.m, 2.09e8*u.m]
        #scaling_zones = [2.09e8*u.m, 4.45e8*u.m, 1.59e12 * u.m]
        runits = (G * M / c**2)**(-1)
        Munits = (1e8 * u.Msun)**(-1)
        fmunits = (1e23 * u.kg / u.s)**(-1)
    elif scale == "galactic":
        scaling_zones = [1.59e3 * u.m, 17.7*u.m, 9.31*u.m]
        runits = (G * M / c**2)**(-1)
        Munits = (1 * u.Msun)**(-1)
        fmunits = (1e14 * u.kg / u.s)**(-1)
    else:
        raise ValueError("scale must be 'AGN' or 'galactic'")

    if zone == 1:
        pindices = [0, 0, 1, 0]
    elif zone == 2:
        pindices = [-1/10, 7/10, 1/5, 21/20]
    elif zone == 3:
        pindices = [-1/10, 3/4, 3/20, 9/8]
    else:
        raise ValueError("zone must be 1, 2, or 3")

    return calc(scaling_zones[zone-1], alpha, M*Munits, fMdot*fmunits, r*runits, pindices=pindices)

# ============================================================================
# KEPLER FREQUENCY
# ============================================================================
def Omega(alpha, M, fMdot, r, zone, scale="AGN"):
    """Kepler frequency (same for all zones)"""
    scale = scale.lower().strip()
    if scale == "agn":
        scaling = 2.03e-3 / u.s
        runits = (G * M / c**2)**(-1)
        Munits = (1e8 * u.Msun)**(-1)
        fmunits = (1e23 * u.kg / u.s)**(-1)
    elif scale == "galactic":
        scaling = 2.03e5 / u.s
        runits = (G * M / c**2)**(-1)
        Munits = (1 * u.Msun)**(-1)
        fmunits = (1e14 * u.kg / u.s)**(-1)
    else:
        raise ValueError("scale must be 'AGN' or 'galactic'")

    pindices = [0, -1, 0, -3/2]
    return calc(scaling, alpha, M*Munits, fMdot*fmunits, r*runits, pindices=pindices)

# ============================================================================
# ORBITAL PERIOD
# ============================================================================
def period(alpha, M, fMdot, r, zone, scale="AGN"):
    """Orbital period 2π/Ω (same for all zones)"""
    scale = scale.lower().strip()
    if scale == "agn":
        scaling = 3100 * u.s
        runits = (G * M / c**2)**(-1)
        Munits = (1e8 * u.Msun)**(-1)
        fmunits = (1e23 * u.kg / u.s)**(-1)
    elif scale == "galactic":
        scaling = 3.1e-5 * u.s
        runits = (G * M / c**2)**(-1)
        Munits = (1 * u.Msun)**(-1)
        fmunits = (1e14 * u.kg / u.s)**(-1)
    else:
        raise ValueError("scale must be 'AGN' or 'galactic'")

    pindices = [0, 1, 0, 3/2]
    return calc(scaling, alpha, M*Munits, fMdot*fmunits, r*runits, pindices=pindices)

# ============================================================================
# RADIAL TIME
# ============================================================================
def radial_time(alpha, M, fMdot, r, zone, scale="AGN"):
    """Radial time r/vr"""
    scale = scale.lower().strip()
    if scale == "agn":
        scaling_zones = [12.7 * u.s, 1.63e8*u.s, 2.47e8*u.s]
        runits = (G * M / c**2)**(-1)
        Munits = (1e8 * u.Msun)**(-1)
        fmunits = (1e23 * u.kg / u.s)**(-1)
    elif scale == "galactic":
        scaling_zones = [1.27e-5 * u.s, 0.103*u.s, 0.124*u.s]
        runits = (G * M / c**2)**(-1)
        Munits = (1 * u.Msun)**(-1)
        fmunits = (1e14 * u.kg / u.s)**(-1)
    else:
        raise ValueError("scale must be 'AGN' or 'galactic'")

    if zone == 1:
        pindices = [-1, 3, -2, 7/2]
    elif zone == 2:
        pindices = [-4/5, 8/5, -2/5, 7/5]
    elif zone == 3:
        pindices = [-4/5, 3/2, -3/10, 5/4]
    else:
        raise ValueError("zone must be 1, 2, or 3")

    return calc(scaling_zones[zone-1], alpha, M*Munits, fMdot*fmunits, r*runits, pindices=pindices)

# ============================================================================
# RADIAL VELOCITY
# ============================================================================
def vr(alpha, M, fMdot, r, zone, scale="AGN"):
    """Radial velocity"""
    scale = scale.lower().strip()
    if scale == "agn":
        scaling_zones = [1.16e10 * u.m/u.s, 0.908e3*u.m/u.s, 0.598e3*u.m/u.s]
        runits = (G * M / c**2)**(-1)
        Munits = (1e8 * u.Msun)**(-1)
        fmunits = (1e23 * u.kg / u.s)**(-1)
    elif scale == "galactic":
        scaling_zones = [1.16e8 * u.m/u.s, 14.4e3*u.m/u.s, 11.9e3*u.m/u.s]
        runits = (G * M / c**2)**(-1)
        Munits = (1 * u.Msun)**(-1)
        fmunits = (1e14 * u.kg / u.s)**(-1)
    else:
        raise ValueError("scale must be 'AGN' or 'galactic'")

    if zone == 1:
        pindices = [1, -2, 2, -5/2]
    elif zone == 2:
        pindices = [4/5, -3/5, 2/5, -2/5]
    elif zone == 3:
        pindices = [4/5, -1/2, 3/10, -1/4]
    else:
        raise ValueError("zone must be 1, 2, or 3")

    return calc(scaling_zones[zone-1], alpha, M*Munits, fMdot*fmunits, r*runits, pindices=pindices)

# ============================================================================
# MASS DENSITY
# ============================================================================
def rho(alpha, M, fMdot, r, zone, scale="AGN"):
    """Mass density"""
    scale = scale.lower().strip()
    if scale == "agn":
        scaling_zones = [2.91e-12 * u.kg/u.m**3, 0.133*u.kg/u.m**3, 0.432*u.kg/u.m**3]
        runits = (G * M / c**2)**(-1)
        Munits = (1e8 * u.Msun)**(-1)
        fmunits = (1e23 * u.kg / u.s)**(-1)
    elif scale == "galactic":
        scaling_zones = [0.0291 * u.kg/u.m**3, 2.11e4*u.kg/u.m**3, 4.85e4*u.kg/u.m**3]
        runits = (G * M / c**2)**(-1)
        Munits = (1 * u.Msun)**(-1)
        fmunits = (1e14 * u.kg / u.s)**(-1)
    else:
        raise ValueError("scale must be 'AGN' or 'galactic'")

    if zone == 1:
        pindices = [-1, 1, -2, 3/2]
    elif zone == 2:
        pindices = [-7/10, -11/10, 2/5, -33/20]
    elif zone == 3:
        pindices = [-7/10, -5/4, 11/20, -15/8]
    else:
        raise ValueError("zone must be 1, 2, or 3")

    return calc(scaling_zones[zone-1], alpha, M*Munits, fMdot*fmunits, r*runits, pindices=pindices)

# ============================================================================
# SURFACE DENSITY
# ============================================================================
def Sigma(alpha, M, fMdot, r, zone, scale="AGN"):
    """Surface density"""
    scale = scale.lower().strip()
    if scale == "agn":
        scaling_zones = [9.27 * u.kg/u.m**2, 1.19e8*u.kg/u.m**2, 1.80e8*u.kg/u.m**2]
        runits = (G * M / c**2)**(-1)
        Munits = (1e8 * u.Msun)**(-1)
        fmunits = (1e23 * u.kg / u.s)**(-1)
    elif scale == "galactic":
        scaling_zones = [92.7 * u.kg/u.m**2, 7.49e5*u.kg/u.m**2, 9.03e5*u.kg/u.m**2]
        runits = (G * M / c**2)**(-1)
        Munits = (1 * u.Msun)**(-1)
        fmunits = (1e14 * u.kg / u.s)**(-1)
    else:
        raise ValueError("scale must be 'AGN' or 'galactic'")

    if zone == 1:
        pindices = [-1, 1, -1, 3/2]
    elif zone == 2:
        pindices = [-4/5, -2/5, 3/5, -3/5]
    elif zone == 3:
        pindices = [-4/5, -1/2, 7/10, -3/4]
    else:
        raise ValueError("zone must be 1, 2, or 3")

    return calc(scaling_zones[zone-1], alpha, M*Munits, fMdot*fmunits, r*runits, pindices=pindices)

# ============================================================================
# TOTAL PRESSURE
# ============================================================================
def pressure(alpha, M, fMdot, r, zone, scale="AGN"):
    """Total pressure"""
    scale = scale.lower().strip()
    if scale == "agn":
        scaling_zones = [1.52e7 * u.Pa, 5.44e10*u.Pa, 3.87e10*u.Pa]
        runits = (G * M / c**2)**(-1)
        Munits = (1e8 * u.Msun)**(-1)
        fmunits = (1e23 * u.kg / u.s)**(-1)
    elif scale == "galactic":
        scaling_zones = [1.52e15 * u.Pa, 1.37e17*u.Pa, 8.67e16*u.Pa]
        runits = (G * M / c**2)**(-1)
        Munits = (1 * u.Msun)**(-1)
        fmunits = (1e14 * u.kg / u.s)**(-1)
    else:
        raise ValueError("scale must be 'AGN' or 'galactic'")

    if zone == 1:
        pindices = [-1, -1, 0, -3/2]
    elif zone == 2:
        pindices = [-9/10, -17/10, 4/5, -51/20]
    elif zone == 3:
        pindices = [-9/10, -7/4, 17/20, -21/8]
    else:
        raise ValueError("zone must be 1, 2, or 3")

    return calc(scaling_zones[zone-1], alpha, M*Munits, fMdot*fmunits, r*runits, pindices=pindices)

# ============================================================================
# PRESSURE RATIO
# ============================================================================
def pressure_ratio(alpha, M, fMdot, r, zone, scale="AGN"):
    """Pressure ratio pgas/prad (dimensionless)"""
    scale = scale.lower().strip()
    if scale == "agn":
        scaling_zones = [1.31e-9, 2.80e-4, 0.0855]
        runits = (G * M / c**2)**(-1)
        Munits = (1e8 * u.Msun)**(-1)
        fmunits = (1e23 * u.kg / u.s)**(-1)
    elif scale == "galactic":
        scaling_zones = [1.31e-5, 0.0111, 1.21]
        runits = (G * M / c**2)**(-1)
        Munits = (1 * u.Msun)**(-1)
        fmunits = (1e14 * u.kg / u.s)**(-1)
    else:
        raise ValueError("scale must be 'AGN' or 'galactic'")

    if zone == 1:
        pindices = [-1/4, 7/4, -2, 21/8]
    elif zone == 2:
        pindices = [-1/10, 7/10, -4/5, 21/20]
    elif zone == 3:
        pindices = [-1/10, 1/4, -7/20, 3/8]
    else:
        raise ValueError("zone must be 1, 2, or 3")

    return calc(scaling_zones[zone-1], alpha, M*Munits, fMdot*fmunits, r*runits, pindices=pindices)

# ============================================================================
# CENTRAL TEMPERATURE
# ============================================================================
def T_central(alpha, M, fMdot, r, zone, scale="AGN"):
    """Central temperature"""
    scale = scale.lower().strip()
    if scale == "agn":
        scaling_zones = [4.96e5 * u.K, 2.96e7*u.K, 6.51e6*u.K]
        runits = (G * M / c**2)**(-1)
        Munits = (1e8 * u.Msun)**(-1)
        fmunits = (1e23 * u.kg / u.s)**(-1)
    elif scale == "galactic":
        scaling_zones = [4.96e7 * u.K, 4.70e8*u.K, 1.30e8*u.K]
        runits = (G * M / c**2)**(-1)
        Munits = (1 * u.Msun)**(-1)
        fmunits = (1e14 * u.kg / u.s)**(-1)
    else:
        raise ValueError("scale must be 'AGN' or 'galactic'")

    if zone == 1:
        pindices = [-1/4, -1/4, 0, -3/8]
    elif zone == 2:
        pindices = [-1/5, -3/5, 2/5, -9/10]
    elif zone == 3:
        pindices = [-1/5, -1/2, 3/10, -3/4]
    else:
        raise ValueError("zone must be 1, 2, or 3")

    return calc(scaling_zones[zone-1], alpha, M*Munits, fMdot*fmunits, r*runits, pindices=pindices)

# ============================================================================
# SURFACE TEMPERATURE
# ============================================================================
def T_surface(alpha, M, fMdot, r, zone, scale="AGN"):
    """Surface temperature (zone 1&2: Modified BB, zone 3: BB)"""
    scale = scale.lower().strip()
    if scale == "agn":
        scaling_zones = [5.17e8 * u.K, 2.21e6*u.K, 9.65e5*u.K]
        runits = (G * M / c**2)**(-1)
        Munits = (1e8 * u.Msun)**(-1)
        fmunits = (1e23 * u.kg / u.s)**(-1)
    elif scale == "galactic":
        scaling_zones = [4.00e9 * u.K, 1.99e8*u.K, 5.43e7*u.K]
        runits = (G * M / c**2)**(-1)
        Munits = (1 * u.Msun)**(-1)
        fmunits = (1e14 * u.kg / u.s)**(-1)
    else:
        raise ValueError("scale must be 'AGN' or 'galactic'")

    if zone == 1:
        pindices = [2/9, -10/9, 8/9, -5/3]
    elif zone == 2:
        pindices = [7/45, -29/45, 16/45, -29/30]
    elif zone == 3:
        pindices = [0, -1/2, 1/4, -3/4]
    else:
        raise ValueError("zone must be 1, 2, or 3")

    return calc(scaling_zones[zone-1], alpha, M*Munits, fMdot*fmunits, r*runits, pindices=pindices)

# ============================================================================
# SURFACE FLUX
# ============================================================================
def Q(alpha, M, fMdot, r, zone, scale="AGN"):
    """Surface flux (same for all zones)"""
    scale = scale.lower().strip()
    if scale == "agn":
        scaling = 4.92e16 * u.W/u.m**2
        runits = (G * M / c**2)**(-1)
        Munits = (1e8 * u.Msun)**(-1)
        fmunits = (1e23 * u.kg / u.s)**(-1)
    elif scale == "galactic":
        scaling = 4.92e23 * u.W/u.m**2
        runits = (G * M / c**2)**(-1)
        Munits = (1 * u.Msun)**(-1)
        fmunits = (1e14 * u.kg / u.s)**(-1)
    else:
        raise ValueError("scale must be 'AGN' or 'galactic'")

    pindices = [0, -2, 1, -3]
    return calc(scaling, alpha, M*Munits, fMdot*fmunits, r*runits, pindices=pindices)

# ============================================================================
# OPTICAL DEPTH (SCATTERING)
# ============================================================================
def tau_es(alpha, M, fMdot, r, zone, scale="AGN"):
    """Optical depth (scattering) - dimensionless"""
    scale = scale.lower().strip()
    if scale == "agn":
        scaling_zones = [0.185, 2.37e6, 3.60e6]
        runits = (G * M / c**2)**(-1)
        Munits = (1e8 * u.Msun)**(-1)
        fmunits = (1e23 * u.kg / u.s)**(-1)
    elif scale == "galactic":
        scaling_zones = [1.85, 1.50e4, 1.81e4]
        runits = (G * M / c**2)**(-1)
        Munits = (1 * u.Msun)**(-1)
        fmunits = (1e14 * u.kg / u.s)**(-1)
    else:
        raise ValueError("scale must be 'AGN' or 'galactic'")

    if zone == 1:
        pindices = [-1, 1, -1, 3/2]
    elif zone == 2:
        pindices = [-4/5, -2/5, 3/5, -3/5]
    elif zone == 3:
        pindices = [-4/5, -1/2, 7/10, -3/4]
    else:
        raise ValueError("zone must be 1, 2, or 3")

    return calc(scaling_zones[zone-1], alpha, M*Munits, fMdot*fmunits, r*runits, pindices=pindices)

# ============================================================================
# TRUE OPTICAL DEPTH
# ============================================================================
def tau_star(alpha, M, fMdot, r, zone, scale="AGN"):
    """True optical depth τ* - dimensionless"""
    scale = scale.lower().strip()
    if scale == "agn":
        scaling_zones = [4.32e-7, 1.02e-7, 354]
        runits = (G * M / c**2)**(-1)
        Munits = (1e8 * u.Msun)**(-1)
        fmunits = (1e23 * u.kg / u.s)**(-1)
    elif scale == "galactic":
        scaling_zones = [1.36e-4, 0.204, 5.61]
        runits = (G * M / c**2)**(-1)
        Munits = (1 * u.Msun)**(-1)
        fmunits = (1e14 * u.kg / u.s)**(-1)
    else:
        raise ValueError("scale must be 'AGN' or 'galactic'")

    if zone == 1:
        pindices = [-17/16, 31/16, -2, 93/32]
    elif zone == 2:
        pindices = [-4/5, -9/10, 1/10, 3/20]
    elif zone == 3:
        pindices = [-4/5, 0, 1/5, 0]
    else:
        raise ValueError("zone must be 1, 2, or 3")

    return calc(scaling_zones[zone-1], alpha, M*Munits, fMdot*fmunits, r*runits, pindices=pindices)

# ============================================================================
# OPTICAL DEPTH RATIO
# ============================================================================
def tau_ratio(alpha, M, fMdot, r, zone, scale="AGN"):
    """Optical depth ratio τff/τes - dimensionless"""
    scale = scale.lower().strip()
    if scale == "agn":
        scaling_zones = [5.40e-12, 1.50e-7, 9.82e-5]
        runits = (G * M / c**2)**(-1)
        Munits = (1e8 * u.Msun)**(-1)
        fmunits = (1e23 * u.kg / u.s)**(-1)
    elif scale == "galactic":
        scaling_zones = [5.40e-9, 1.50e-6, 3.11e-4]
        runits = (G * M / c**2)**(-1)
        Munits = (1 * u.Msun)**(-1)
        fmunits = (1e14 * u.kg / u.s)**(-1)
    else:
        raise ValueError("scale must be 'AGN' or 'galactic'")

    if zone == 1:
        pindices = [-1/8, 15/8, -2, 45/16]
    elif zone == 2:
        pindices = [0, 1, -1, 3/2]
    elif zone == 3:
        pindices = [0, 1/2, -1/2, 3/4]
    else:
        raise ValueError("zone must be 1, 2, or 3")

    return calc(scaling_zones[zone-1], alpha, M*Munits, fMdot*fmunits, r*runits, pindices=pindices)

# ============================================================================
# MAXIMUM MAGNETIC FIELD
# ============================================================================
def B_max(alpha, M, fMdot, r, zone, scale="AGN"):
    """Maximum magnetic field"""
    scale = scale.lower().strip()
    if scale == "agn":
        scaling_zones = [6.18 * u.T, 370*u.T, 312*u.T]
        runits = (G * M / c**2)**(-1)
        Munits = (1e8 * u.Msun)**(-1)
        fmunits = (1e23 * u.kg / u.s)**(-1)
    elif scale == "galactic":
        scaling_zones = [6.18e4 * u.T, 5.86e5*u.T, 4.67e5*u.T]
        runits = (G * M / c**2)**(-1)
        Munits = (1 * u.Msun)**(-1)
        fmunits = (1e14 * u.kg / u.s)**(-1)
    else:
        raise ValueError("scale must be 'AGN' or 'galactic'")

    if zone == 1:
        pindices = [0, -1/2, 0, -3/4]
    elif zone == 2:
        pindices = [1/20, -17/20, 2/5, -51/40]
    elif zone == 3:
        pindices = [-9/20, -7/8, 17/40, -21/16]
    else:
        raise ValueError("zone must be 1, 2, or 3")

    return calc(scaling_zones[zone-1], alpha, M*Munits, fMdot*fmunits, r*runits, pindices=pindices)
