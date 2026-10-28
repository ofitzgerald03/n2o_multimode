"""
Thrust efficiency eta_T(Te, R_m), Petro and Sedwick (2017) Eqs. (16)-(18),
and their argon ionisation-cost fit, Eqs. (9)-(10).
"""

import numpy as np


def thrust_efficiency(Te_eV, delta_phi_sh_eV, Ec_eV, R_m):
    """Thrust efficiency, Eq. (18) of Petro and Sedwick (2017):

        eta_T = |dphi_sh| / [(1 + 1/R_m)(2 Te + |dphi_sh| + Ec)
                             - (Te / R_m) ln R_m]

    Energies in eV; R_m is the magnetic mirror ratio.
    """
    Te = np.asarray(Te_eV, dtype=float)
    dphi = np.asarray(delta_phi_sh_eV, dtype=float)
    Ec = np.asarray(Ec_eV, dtype=float)

    numerator = dphi
    denominator = (1.0 + 1.0 / R_m) * (2.0 * Te + dphi + Ec) - (Te / R_m) * np.log(R_m)
    return numerator / denominator


def delta_phi_mirror(delta_phi_sh_eV, Te_eV, R_m):
    """Sheath drop with a magnetic mirror, Eq. (16):
    |dphi_mirror| = |dphi_sh| - Te ln R_m."""
    return delta_phi_sh_eV - Te_eV * np.log(R_m)


def thrust_efficiency_from_eq17(Te_eV, delta_phi_sh_eV, Ec_eV, R_m):
    """Thrust efficiency in the form of Eq. (17), with A_wall = A_0 / R_m:

    eta_T = |dphi_sh| / [(1/R_m)(2 Te + |dphi_mirror| + Ec)
                         + (2 Te + |dphi_sh| + Ec)]
    """
    Te = np.asarray(Te_eV, dtype=float)
    dphi = np.asarray(delta_phi_sh_eV, dtype=float)
    Ec = np.asarray(Ec_eV, dtype=float)
    dphi_m = delta_phi_mirror(dphi, Te, R_m)

    num = dphi
    den = (1.0 / R_m) * (2.0 * Te + dphi_m + Ec) + (2.0 * Te + dphi + Ec)
    return num / den


# Coefficients of Eqs. (9) and (10)
PS_ARGON_LOW = dict(
    c0=72.78,
    c1=2.673,
    c2=-3.195,
    c3=-0.01939,
    c4=5.606,
    c5=-25.86,
)
PS_ARGON_HIGH = dict(
    k0=-42.42,
    k1=-0.935,
    k2=60.58,
    k3=579.3,
    k4=-4.316,
)


def PS_Ar_curve_fit(Te):
    """Argon ionisation cost [eV], Eqs. (9)-(10) of Petro and Sedwick (2017).

    Te <= 3 eV:  Ec = c0 exp(c1 Te^c2) + c3 Te^c4 + c5
    Te > 3 eV:   Ec = k0 Te / (Te - k1) + k2 Te / (Te + k1)
                      + (Te - k3) / (Te - k4)^2
    """
    Te = np.asarray(Te, dtype=float)
    out = np.zeros_like(Te)

    low = Te <= 3.0
    high = Te > 3.0

    c = PS_ARGON_LOW
    out[low] = (
        c["c0"] * np.exp(c["c1"] * Te[low] ** c["c2"])
        + c["c3"] * Te[low] ** c["c4"]
        + c["c5"]
    )

    k = PS_ARGON_HIGH
    out[high] = (
        k["k0"] * Te[high] / (Te[high] - k["k1"])
        + k["k2"] * Te[high] / (Te[high] + k["k1"])
        + (Te[high] - k["k3"]) / (Te[high] - k["k4"]) ** 2
    )

    return out
