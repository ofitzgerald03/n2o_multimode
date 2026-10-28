"""
Exhaust quantities: Bohm velocity, exit velocity, sheath drop and specific
impulse, for single-species and multi-species propellants.
"""

# Multi-species closures (Baalrud and Hegna 2011):
#     case 'A': each ion species leaves at its own Bohm velocity;
#     case 'B': all species share the system Bohm velocity.

import numpy as np

from ..constants import E_CHARGE, M_ELECTRON, G0
from .cross_sections import (
    compute_rate_constants_vs_Te,
    compute_ion_fractions,
)


# Single species


def bohm_velocity_individual(Te_eV, m_kg):
    """Bohm velocity v_b = sqrt(e Te / m_i) [m/s]; Te in eV, m in kg."""
    return np.sqrt(E_CHARGE * Te_eV / m_kg)


def exit_velocity_PS(Te_eV, m_kg):
    """Exit velocity, Eq. (7) of Petro and Sedwick (2017):
    u_e = sqrt(e Te (1 + ln(m_i / (2 pi m_e))) / m_i)."""
    log_term = 1.0 + np.log(m_kg / (2.0 * np.pi * M_ELECTRON))
    return np.sqrt(E_CHARGE * Te_eV * log_term / m_kg)


def sheath_potential_PS(Te_eV, m_kg):
    """Sheath potential drop [eV], Eq. (6) of Petro and Sedwick (2017):
    |dphi_sh| = (Te / 2) (1 + ln(m_i / (2 pi m_e)))."""
    log_term = 1.0 + np.log(m_kg / (2.0 * np.pi * M_ELECTRON))
    return 0.5 * Te_eV * log_term


def isp_single_species(Te_eV, m_kg):
    """Specific impulse u_e / g0 [s]."""
    return exit_velocity_PS(Te_eV, m_kg) / G0


# Multi-species


def bohm_velocity_system(Te_eV, fractions, masses_kg):
    """System Bohm velocity, Baalrud and Hegna (2011) Eq. (5):
    v_b,sys = sqrt(sum_i f_i e Te / m_i)."""
    sum_term = np.zeros_like(Te_eV, dtype=float)
    for sp in fractions:
        sum_term += fractions[sp] / masses_kg[sp]
    return np.sqrt(E_CHARGE * Te_eV * sum_term)


def effective_isp_closure_A(Te_eV, fractions, masses_kg):
    """Specific impulse with species-resolved Bohm velocities:
    Isp = sum_i(f_i m_i v_b,i u_e,i) / (g0 sum_i(f_i m_i v_b,i))."""
    thrust_num = np.zeros_like(Te_eV, dtype=float)
    mdot_num = np.zeros_like(Te_eV, dtype=float)

    for sp in fractions:
        m = masses_kg[sp]
        f = fractions[sp]
        vb = bohm_velocity_individual(Te_eV, m)
        ue = exit_velocity_PS(Te_eV, m)
        thrust_num += f * m * vb * ue
        mdot_num += f * m * vb

    return thrust_num / (G0 * mdot_num)


def effective_isp_closure_B(Te_eV, fractions, masses_kg):
    """Specific impulse with a common Bohm velocity:
    Isp = sum_i(f_i m_i u_e,i) / (g0 sum_i(f_i m_i))."""
    thrust_num = np.zeros_like(Te_eV, dtype=float)
    mdot_num = np.zeros_like(Te_eV, dtype=float)

    for sp in fractions:
        m = masses_kg[sp]
        f = fractions[sp]
        ue = exit_velocity_PS(Te_eV, m)
        thrust_num += f * m * ue
        mdot_num += f * m

    return thrust_num / (G0 * mdot_num)


def effective_sheath_potential_closure_A(Te_eV, fractions, masses_kg):
    """Particle-flux-weighted sheath drop, species-resolved Bohm velocities:
    sum_i(f_i v_b,i |dphi_i|) / sum_i(f_i v_b,i)."""
    num = np.zeros_like(Te_eV, dtype=float)
    den = np.zeros_like(Te_eV, dtype=float)

    for sp in fractions:
        m = masses_kg[sp]
        f = fractions[sp]
        vb = bohm_velocity_individual(Te_eV, m)
        phi = sheath_potential_PS(Te_eV, m)
        num += f * vb * phi
        den += f * vb

    return num / den


def effective_sheath_potential_closure_B(Te_eV, fractions, masses_kg):
    """Fraction-weighted sheath drop with a common Bohm velocity:
    sum_i f_i |dphi_i|."""
    result = np.zeros_like(Te_eV, dtype=float)
    for sp in fractions:
        m = masses_kg[sp]
        f = fractions[sp]
        result += f * sheath_potential_PS(Te_eV, m)
    return result


def compute_fractions_for_propellant(ion_channels, ion_energies, Te_array):
    """Rate constants and ion fractions of a multi-species propellant.

    Returns ({species: f_i(Te)}, {species: kappa_i(Te) [m^3/s]}).
    """
    kappas = {}
    for sp, data in ion_channels.items():
        kappas[sp] = compute_rate_constants_vs_Te(ion_energies, data["sigma"], Te_array)

    fractions = compute_ion_fractions(kappas)
    return fractions, kappas


def _fraction_weighted_mass(fractions, masses_kg):
    """Mean ion mass sum_i f_i m_i [kg]."""
    Te_like = next(iter(fractions.values()))
    m_avg = np.zeros_like(Te_like, dtype=float)
    for sp in fractions:
        m_avg += fractions[sp] * masses_kg[sp]
    return m_avg


def _closure_A_flux_weighted(Te_eV, fractions, masses_kg):
    """Flux-weighted Bohm velocity and ion mass, species-resolved:
    v_b,eff = sum_i f_i v_b,i;  m_eff = sum_i(f_i v_b,i m_i) / v_b,eff."""
    Te_like = next(iter(fractions.values()))
    v_b_eff = np.zeros_like(Te_like, dtype=float)
    mflux_eff = np.zeros_like(Te_like, dtype=float)
    for sp in fractions:
        v_bi = bohm_velocity_individual(Te_eV, masses_kg[sp])
        v_b_eff += fractions[sp] * v_bi
        mflux_eff += fractions[sp] * v_bi * masses_kg[sp]
    m_eff = np.where(v_b_eff > 0, mflux_eff / v_b_eff, 0.0)
    return v_b_eff, m_eff


def _closure_B_system_bohm_and_mass(Te_eV, fractions, masses_kg):
    """System Bohm velocity and mean ion mass, common closure."""
    v_b_sys = bohm_velocity_system(Te_eV, fractions, masses_kg)
    return v_b_sys, _fraction_weighted_mass(fractions, masses_kg)


def compute_per_propellant_quantities(
    Te_eV,
    propellant_label,  # 'single' or 'multi'
    single_species_mass_kg=None,
    fractions=None,
    species_masses_kg=None,
    case="A",
):
    """Per-particle quantities of one propellant.

    'single': pass single_species_mass_kg. 'multi': pass fractions,
    species_masses_kg and case ('A' or 'B').

    Returns a dict with 'I_sp' [s], 'u_e' [m/s], 'v_b_flux' [m/s],
    'phi_sh' [eV] and 'm_i_eff' [kg].
    """
    Te = np.asarray(Te_eV, dtype=float)

    if propellant_label == "single":
        if single_species_mass_kg is None:
            raise ValueError("single_species_mass_kg required for 'single'")
        m = single_species_mass_kg
        v_b = bohm_velocity_individual(Te, m)
        u_e = exit_velocity_PS(Te, m)
        Isp = isp_single_species(Te, m)
        phi = sheath_potential_PS(Te, m)
        return {
            "I_sp": Isp,
            "u_e": u_e,
            "v_b_flux": v_b,
            "phi_sh": phi,
            "m_i_eff": np.full_like(Te, m),
        }

    elif propellant_label == "multi":
        if fractions is None or species_masses_kg is None:
            raise ValueError("fractions and species_masses_kg required for 'multi'")
        if case == "A":
            Isp = effective_isp_closure_A(Te, fractions, species_masses_kg)
            phi = effective_sheath_potential_closure_A(Te, fractions, species_masses_kg)
            v_b_flux, m_eff = _closure_A_flux_weighted(Te, fractions, species_masses_kg)
        elif case == "B":
            Isp = effective_isp_closure_B(Te, fractions, species_masses_kg)
            phi = effective_sheath_potential_closure_B(Te, fractions, species_masses_kg)
            v_b_flux, m_eff = _closure_B_system_bohm_and_mass(
                Te, fractions, species_masses_kg
            )
        else:
            raise ValueError(f"Unknown case: {case!r}")

        u_e = Isp * G0
        return {
            "I_sp": Isp,
            "u_e": u_e,
            "v_b_flux": v_b_flux,
            "phi_sh": phi,
            "m_i_eff": m_eff,
        }

    else:
        raise ValueError(f"Unknown propellant_label: {propellant_label!r}")
