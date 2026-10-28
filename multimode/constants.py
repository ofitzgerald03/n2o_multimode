"""
Physical constants (CODATA 2018) and atomic masses (NIST).
"""

# CODATA 2018
E_CHARGE = 1.602176634e-19  # Elementary charge [C] (also J/eV conversion factor)
M_ELECTRON = 9.1093837015e-31  # Electron mass [kg]
AMU = 1.66053906660e-27  # Atomic mass unit [kg]
G0 = 9.80665  # Standard gravitational acceleration [m/s^2] (IAU)

# Earth (IERS 2010 / EGM2008), for the LEO orbital period
MU_EARTH = 3.986004418e14  # Gravitational parameter [m^3/s^2]
R_EARTH = 6378.137e3  # Equatorial radius [m]

# Atomic masses (NIST Atomic Weights and Isotopic Compositions)
# These are the mass values most commonly needed at module level.
M_XE_AMU = 131.293
M_AR_AMU = 39.948
M_XE = M_XE_AMU * AMU
M_AR = M_AR_AMU * AMU
