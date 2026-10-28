"""
Component mass models: tanks, power electronics, added power system, feed
system and structure margin.
"""

# Tank mass:
#
#     m_tank = max(P_MEOP V / (g0 PF), m_min) + m_fixed
#     V      = m_prop * residual / rho_min * ullage
#
# with PF the performance factor of the tank type (metal or COPV), rho_min
# the propellant density at the maximum storage temperature, m_min the
# minimum tank mass and m_fixed the fixed mass per tank. An installed empty
# tank keeps m_min + m_fixed.

from ..constants import G0
from .parameters import value, power_margin_w

_PROP_ROWS = {
    "n2o": ("p_meop_n2o", "rho_min_n2o"),
    "propene": ("p_meop_propene", "rho_min_propene"),
    "xe": ("p_meop_xe", "rho_min_xe"),
}

# Tank type per propellant; None takes the 'chem_tank_technology' row
_DEFAULT_TECH = {
    "n2o": None,
    "propene": None,
    "xe": "copv",
}


def copv_pf():
    """COPV performance factor [m] that reproduces the 4.5 % Xe tank
    fraction: alpha = P residual ullage / (rho g0 PF)."""
    ru = value("residual_factor") * value("ullage_factor")
    return value("p_meop_xe") / (value("rho_min_xe") / ru * G0 * value("alpha_xe_copv"))


def performance_factor(technology):
    """Performance factor [m] of a tank type ('metal' | 'copv')."""
    if technology == "metal":
        return value("pf_metal")
    if technology == "copv":
        return copv_pf()
    raise ValueError(f"Unknown tank technology {technology!r}")


def tank_technology(propellant, overrides=None):
    """Tank type of a propellant; overrides, e.g. {'n2o': 'copv'}, take
    precedence."""
    if overrides and propellant in overrides:
        return overrides[propellant]
    tech = _DEFAULT_TECH[propellant]
    if tech is None:
        tech = value("chem_tank_technology")
    return tech


def tank_mass(m_prop, propellant, technology, installed=True):
    """Tank mass [kg]; an installed empty tank keeps its minimum mass."""
    if m_prop <= 0.0:
        return (value("m_min_tank") + value("m_fixed_tank")) if installed else 0.0
    p_key, rho_key = _PROP_ROWS[propellant]
    V = m_prop * value("residual_factor") / value(rho_key) * value("ullage_factor")
    m_pv = value(p_key) * V / (G0 * performance_factor(technology))
    return max(m_pv, value("m_min_tank")) + value("m_fixed_tank")


def floor_crossover(propellant, technology="metal"):
    """Propellant load [kg] at which the pressure-volume term reaches the
    minimum tank mass."""
    p_key, rho_key = _PROP_ROWS[propellant]
    ru = value("residual_factor") * value("ullage_factor")
    return (
        value("m_min_tank")
        * G0
        * performance_factor(technology)
        * value(rho_key)
        / (value(p_key) * ru)
    )


def electronics_mass(p_dc_rated_w, kind):
    """Power-electronics mass [kg] on rated DC power; kind 'ppu', 'rf' or
    None."""
    if kind is None or p_dc_rated_w <= 0.0:
        return 0.0
    sm = value("sm_ppu") if kind == "ppu" else value("sm_rf")
    return sm * p_dc_rated_w / 1e3


def eps_increment(p_dc_peak_w, bus_class, arr_degr=0.0):
    """Added power-system mass [kg]:
    sm_eps_inc max(0, P_dc,peak - margin) / (1 - arr_degr).

    arr_degr : solar-array power loss over the transfer, in [0, 1); the
    array is sized for the end of the transfer.
    """
    if not 0.0 <= arr_degr < 1.0:
        raise ValueError(f"arr_degr must lie in [0, 1); got {arr_degr!r}")
    return (
        value("sm_eps_inc")
        * max(0.0, p_dc_peak_w - power_margin_w(bus_class))
        / 1e3
        / (1.0 - arr_degr)
    )


def transfer_array_derate(mission, ep_transfer_share):
    """Array power loss for a mission: arr_degr_transfer times the electric
    share of the transfer; zero without a transfer class."""
    if not any(c.kind == "transfer" for c in mission.classes):
        return 0.0
    return value("arr_degr_transfer") * ep_transfer_share


def eps_mass_cap(bus_mass):
    """Power-system limit as a mass [kg]: eps_cap_frac x bus mass."""
    return value("eps_cap_frac") * bus_mass


def feed_mass(m_prop_total):
    """Feed-system mass [kg]."""
    return value("feed_frac") * m_prop_total


def with_structure(dry_subtotal):
    """Dry mass with the structure margin."""
    return dry_subtotal * (1.0 + value("struct_frac"))
