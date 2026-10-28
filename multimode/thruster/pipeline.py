"""
Thruster-model pipeline: ionisation cost, sheath drop, ion fractions, thrust
efficiency and the composed thrust, mass flow and power of each propellant
against the xenon anchor.
"""

# Multi-species propellants carry both Bohm closures: suffix _A for
# species-resolved Bohm velocities, _B for the common system Bohm velocity.
# Single-species propellants use the _A keys.

import numpy as np

from ..config import PARAMS
from .propellants import ALL_PROPELLANTS, Propellant
from .exhaust import (
    sheath_potential_PS,
    effective_sheath_potential_closure_A,
    effective_sheath_potential_closure_B,
    compute_fractions_for_propellant,
)
from .cross_sections import compute_ion_fractions_species_resolved
from .ionisation_cost import solve_Ec
from .thrust_efficiency import thrust_efficiency
from .exhaust import compute_per_propellant_quantities
from .anchoring import (
    AnchorReference,
    framework_extensive_arrays,
    compute_scale_factors,
    compose_extensive_curves,
    at_anchor_operating_point,
)


def run_pipeline(
    params: dict = None,
    propellants: tuple = ("Xe", "N2O"),
    anchor_propellant: str = "Xe",
) -> dict:
    """Run the thruster model for the given propellants.

    params : PARAMS-style dict (default config.PARAMS).
    propellants : labels to compute; must include the anchor propellant.
    anchor_propellant : single-species propellant whose measured operating
        point sets the absolute scale.

    Returns {'Te', 'params', 'anchor', <label>: {...}} where each propellant
    entry holds 'Ec', 'K_by_name', 'kappa_ion', 'channels', 'phi_sh_A/B',
    'fractions', 'q_A/B' (per-particle bundles), 'eta_T_A/B' (at
    Rm_nominal), 'eta_T_A/B_at_Rm_ref', 'fw_A/B' (framework arrays),
    'ext_A/B' (composed curves) and 'op_point_A/B' (at the anchor). The _B
    entries are None for single-species propellants.
    """
    if params is None:
        params = PARAMS

    Te = np.arange(
        params["Te_min"], params["Te_max"] + params["Te_step"] / 2, params["Te_step"]
    )

    # Ar takes the Xe value, H2O the N2O value
    eta_u_by_label = {
        "Xe": params.get("eta_u_xe", 1.0),
        "Ar": params.get("eta_u_xe", 1.0),
        "N2O": params.get("eta_u_n2o", 1.0),
        "H2O": params.get("eta_u_n2o", 1.0),
    }

    Rm_nominal = params["Rm_nominal"]
    Rm_ref = params["Rm_ref"]
    eta_RF = params["eta_RF"]
    eta_D = params["eta_D"]

    results = {
        "Te": Te,
        "params": dict(params),
    }

    # Ionisation cost, sheath drop and ion fractions
    for label in propellants:
        prop = ALL_PROPELLANTS[label]
        prop_results = _compute_intensive(prop, Te)
        results[label] = prop_results

    # Per-particle quantities
    for label in propellants:
        prop = ALL_PROPELLANTS[label]
        results[label].update(_compute_per_particle(prop, Te, results[label]))

    # Thrust efficiency at Rm_nominal and Rm_ref
    for label in propellants:
        Ec = results[label]["Ec"]
        phi_A = results[label]["phi_sh_A"]
        phi_B = results[label]["phi_sh_B"]

        results[label]["eta_T_A"] = thrust_efficiency(Te, phi_A, Ec, Rm_nominal)
        results[label]["eta_T_A_at_Rm_ref"] = thrust_efficiency(
            Te,
            phi_A,
            Ec,
            Rm_ref,
        )
        if phi_B is not None:
            results[label]["eta_T_B"] = thrust_efficiency(
                Te,
                phi_B,
                Ec,
                Rm_nominal,
            )
            results[label]["eta_T_B_at_Rm_ref"] = thrust_efficiency(
                Te,
                phi_B,
                Ec,
                Rm_ref,
            )
        else:
            results[label]["eta_T_B"] = None
            results[label]["eta_T_B_at_Rm_ref"] = None

    # Composition against the anchor
    anchor_ref = _build_anchor_reference(params)

    if results[anchor_propellant]["q_B"] is not None:
        raise NotImplementedError(
            f"Anchor propellant {anchor_propellant!r} must be single-species."
        )
    fw_xe = framework_extensive_arrays(
        results[anchor_propellant]["q_A"],
        results[anchor_propellant]["eta_T_A_at_Rm_ref"],
        eta_RF,
        eta_D,
    )
    scale_factors = compute_scale_factors(fw_xe, Te, anchor_ref)
    i_anchor = scale_factors["i_anchor"]

    results["anchor"] = {
        "reference": anchor_ref,
        "anchor_propellant": anchor_propellant,
        "i_anchor": i_anchor,
        "scale_factors": scale_factors,
    }

    for label in propellants:
        bundle = results[label]
        eta_u_target = eta_u_by_label.get(label, 1.0)

        fw_A = framework_extensive_arrays(
            bundle["q_A"],
            bundle["eta_T_A_at_Rm_ref"],
            eta_RF,
            eta_D,
        )
        composed_A = compose_extensive_curves(fw_A, scale_factors, eta_u_target)
        op_point_A = at_anchor_operating_point(
            composed_A,
            fw_xe,
            fw_A,
            anchor_ref,
            eta_u_target,
            i_anchor,
        )
        composed_A["eta_tot"] = bundle["eta_T_A_at_Rm_ref"] * eta_RF * eta_D

        bundle["fw_A"] = fw_A
        bundle["ext_A"] = composed_A
        bundle["op_point_A"] = op_point_A

        if bundle["q_B"] is not None and bundle["eta_T_B_at_Rm_ref"] is not None:
            fw_B = framework_extensive_arrays(
                bundle["q_B"],
                bundle["eta_T_B_at_Rm_ref"],
                eta_RF,
                eta_D,
            )
            composed_B = compose_extensive_curves(
                fw_B,
                scale_factors,
                eta_u_target,
            )
            op_point_B = at_anchor_operating_point(
                composed_B,
                fw_xe,
                fw_B,
                anchor_ref,
                eta_u_target,
                i_anchor,
            )
            composed_B["eta_tot"] = bundle["eta_T_B_at_Rm_ref"] * eta_RF * eta_D

            bundle["fw_B"] = fw_B
            bundle["ext_B"] = composed_B
            bundle["op_point_B"] = op_point_B
        else:
            bundle["fw_B"] = None
            bundle["ext_B"] = None
            bundle["op_point_B"] = None

    return results


def _build_anchor_reference(params: dict) -> AnchorReference:
    """AnchorReference from a PARAMS-style dict."""
    return AnchorReference(
        T_meas=params["T_ref"],
        mdot_meas=params["mdot_ref"],
        Pin_meas=params["Pin_ref"],
        Isp_meas=params["Isp_ref"],
        eta_u_meas=params["eta_u_xe_ref"],
        Te_ref=params["Te_ref"],
        Rm_ref=params["Rm_ref"],
    )


def _compute_intensive(prop: Propellant, Te: np.ndarray) -> dict:
    """Ionisation cost, channel rate coefficients, sheath drop and ion
    fractions of one propellant."""
    if prop.channel_builder is None:
        raise RuntimeError(
            f"Propellant {prop.label!r} has no channel_builder. "
            "thruster.ionisation_cost has not registered the channel builders."
        )
    channels, ion_channels = prop.channel_builder()
    Ec, K_by_name, kappa_ion = solve_Ec(channels, ion_channels, Te)

    common = {
        "Ec": Ec,
        "K_by_name": K_by_name,
        "kappa_ion": kappa_ion,
        "channels": channels,
    }

    if prop.is_multi_species:
        # Ion fractions from the balance n_i ~ kappa_i / v_b,i: kappa ratio
        # with a common Bohm velocity (B), sqrt(m)-weighted with
        # species-resolved velocities (A).
        fractions_B, kappas = compute_fractions_for_propellant(
            prop.ion_channels,
            prop.ion_energies,
            Te,
        )
        masses_amu = {sp: prop.ion_channels[sp]["mass_amu"] for sp in prop.ion_channels}
        fractions_A = compute_ion_fractions_species_resolved(kappas, masses_amu)
        phi_A = effective_sheath_potential_closure_A(
            Te,
            fractions_A,
            prop.species_masses_kg,
        )
        phi_B = effective_sheath_potential_closure_B(
            Te,
            fractions_B,
            prop.species_masses_kg,
        )
        return {
            **common,
            "phi_sh_A": phi_A,
            "phi_sh_B": phi_B,
            "fractions": fractions_B,
            "fractions_A": fractions_A,
        }
    else:
        phi = sheath_potential_PS(Te, prop.mass_kg)
        return {
            **common,
            "phi_sh_A": phi,
            "phi_sh_B": None,
            "fractions": None,
        }


def _compute_per_particle(prop: Propellant, Te: np.ndarray, intensive: dict) -> dict:
    """Per-particle bundles q_A and q_B (None for a single species)."""
    if prop.is_multi_species:
        q_A = compute_per_propellant_quantities(
            Te,
            "multi",
            fractions=intensive["fractions_A"],
            species_masses_kg=prop.species_masses_kg,
            case="A",
        )
        q_B = compute_per_propellant_quantities(
            Te,
            "multi",
            fractions=intensive["fractions"],
            species_masses_kg=prop.species_masses_kg,
            case="B",
        )
        return {"q_A": q_A, "q_B": q_B}
    else:
        q = compute_per_propellant_quantities(
            Te,
            "single",
            single_species_mass_kg=prop.mass_kg,
        )
        return {"q_A": q, "q_B": None}


def baseline_at_anchor(params: dict) -> dict:
    """Xe and N2O at eta_u = 1: thrust, ion mass flow and absorbed power at
    Te_ref and over the Te grid, and the specific impulse u_e / g0.

    Returns {'Te_anchor_eV', 'Te_grid', 'Xe': {...}, 'N2O': {...}} with
    'T_N', 'm_dot_ion', 'P_in_W' (at the anchor) and 'T_N_arr',
    'm_dot_ion_arr', 'P_in_arr', 'I_sp_eta_u_1_arr' (over the grid).
    """
    p = dict(params)
    p["eta_u_xe"] = 1.0
    p["eta_u_n2o"] = 1.0

    results = run_pipeline(
        params=p,
        propellants=("Xe", "N2O"),
        anchor_propellant="Xe",
    )
    Te_grid = results["Te"]
    i = int(np.argmin(np.abs(Te_grid - p["Te_ref"])))

    out = {
        "Te_anchor_eV": float(Te_grid[i]),
        "Te_grid": np.asarray(Te_grid, dtype=float),
    }
    for label in ("Xe", "N2O"):
        if label == "N2O":
            ext = results[label]["ext_A" if p["n2o_case"] == "A" else "ext_B"]
            q = results[label]["q_A" if p["n2o_case"] == "A" else "q_B"]
        else:
            ext = results[label]["ext_A"]
            q = results[label]["q_A"]
        out[label] = {
            "T_N": float(ext["T"][i]),
            "m_dot_ion": float(ext["m_dot_ion"][i]),
            "P_in_W": float(ext["P_in"][i]),
            "T_N_arr": np.asarray(ext["T"], dtype=float),
            "m_dot_ion_arr": np.asarray(ext["m_dot_ion"], dtype=float),
            "P_in_arr": np.asarray(ext["P_in"], dtype=float),
            "I_sp_eta_u_1_arr": np.asarray(q["I_sp"], dtype=float),
        }
    return out
