"""
Ionisation cost E_c(T_e), Eq. (8) of Petro and Sedwick (2017), and the
energy-loss channel lists of each propellant.
"""

# A channel is a dict with keys
#     name      : display label
#     category  : 'elastic' | 'ionisation' | 'electronic' | 'rotation' |
#                 'vibration' | 'dissociation' | 'radiation' | 'attachment'
#     E_j       : energy loss per collision [eV], or a string where it
#                 depends on the electron energy
#     compute_K : callable, Te array -> K_j = E_j kappa_j [eV m^3/s]

import os
import re

import numpy as np

from .cross_sections import (
    compute_K_constant_Ej,
    compute_K_elastic,
    compute_K_Ee_captured,
    compute_rate_constants_vs_Te,
    parse_lxcat_file,
    sigma_m2_to_1e16cm2,
    load_vib_csv,
)
from .propellant_data import n2o as _n2o
from .propellant_data import h2o as _h2o


def solve_Ec(channels, ionisation_channels, Te_array):
    """E_c = sum_j K_j / kappa_ion, Eq. (8) of Petro and Sedwick (2017).

    Returns (Ec [eV], {name: K_j}, kappa_ion [m^3/s]).
    """
    K_by_name = {}
    numerator = np.zeros_like(Te_array)
    for ch in channels:
        K_j = ch["compute_K"](Te_array)
        K_by_name[ch["name"]] = K_j
        numerator += K_j

    kappa_ion = np.zeros_like(Te_array)
    for ch in ionisation_channels:
        kappa_ion += compute_rate_constants_vs_Te(ch["E_grid"], ch["sigma"], Te_array)

    Ec = np.where(kappa_ion > 0, numerator / kappa_ion, np.nan)
    return Ec, K_by_name, kappa_ion


def fractional_loss_by_category(K_by_name, channels, Te, merge_elec_diss_rad=True):
    """Fraction of the energy loss in each channel category,
    chi_j = K_category / sum K, as in Fig. 5 of Petro and Sedwick (2017).

    With merge_elec_diss_rad, 'electronic', 'dissociation' and 'radiation'
    are merged into 'electronic excitation'. Returns {category: chi(Te)}.
    """
    categories = {}
    for ch in channels:
        cat = ch["category"]
        if cat not in categories:
            categories[cat] = np.zeros_like(Te)
        categories[cat] += K_by_name[ch["name"]]

    if merge_elec_diss_rad:
        merged = np.zeros_like(Te)
        for cat in ["electronic", "dissociation", "radiation"]:
            if cat in categories:
                merged += categories.pop(cat)
        if merged.max() > 0:
            categories["electronic excitation"] = merged

    total = sum(categories.values())
    return {
        cat: np.where(total > 0, arr / total, 0.0) for cat, arr in categories.items()
    }


# N2O


def build_n2o_channels(vib_dir: str = "."):
    """N2O channel list: five ionisation channels, momentum-transfer
    elastic, three vibrational modes (from the CSVs in `vib_dir`), three
    swarm-derived electronic channels and O- attachment.

    Returns (channels, ionisation_channels).
    """
    channels = []
    ionisation_channels = []

    # Ionisation, constant E_j
    for species, data in _n2o.N2O_ION_CHANNELS.items():
        Ej = _n2o.N2O_ION_Ej[species]
        channels.append(
            {
                "name": f"Ionisation: {species}",
                "category": "ionisation",
                "E_j": Ej,
                "compute_K": lambda Te, _s=data["sigma"], _Ej=Ej: compute_K_constant_Ej(
                    _n2o.N2O_ION_ENERGIES, _s, _Ej, Te
                ),
            }
        )
        ionisation_channels.append(
            {
                "name": f"Ionisation: {species}",
                "sigma": data["sigma"],
                "E_grid": _n2o.N2O_ION_ENERGIES,
            }
        )

    # Elastic, E_j = 2 E_e m_e / m_N2O
    channels.append(
        {
            "name": "Elastic (MT)",
            "category": "elastic",
            "E_j": "2*E_e*(m_e/m_N2O)",
            "compute_K": lambda Te: compute_K_elastic(
                _n2o.N2O_ELASTIC_MT_ENERGIES, _n2o.N2O_ELASTIC_MTCS, _n2o.M_N2O, Te
            ),
        }
    )

    # Vibrational excitation, constant E_j
    vib_files = {
        "nu2": "vibration_nu2.csv",
        "nu3_2nu2": "vibration_nu3_2nu2.csv",
        "nu1_2nu3": "vibration_nu1_2nu3.csv",
    }
    vib_labels = {
        "nu2": "Vibration: ν₂ (bend)",
        "nu3_2nu2": "Vibration: ν₃/2ν₂",
        "nu1_2nu3": "Vibration: ν₁/2ν₃",
    }
    for mode, filename in vib_files.items():
        E_vib, sigma_vib = load_vib_csv(os.path.join(vib_dir, filename))
        Ej = _n2o.N2O_VIB_Ej[mode]
        channels.append(
            {
                "name": vib_labels[mode],
                "category": "vibration",
                "E_j": Ej,
                "compute_K": lambda Te, _E=E_vib, _s=sigma_vib, _Ej=Ej: (
                    compute_K_constant_Ej(_E, _s, _Ej, Te)
                ),
            }
        )

    # Electronic excitation, E_j = threshold
    elec_data = {
        f"swarm {d['threshold']} eV": {
            "sigma": d["sigma"],
            "Ej": d["threshold"],
            "E": d["energies"],
        }
        for d in _n2o.N2O_ELEC_SWARM_SPLIT.values()
    }
    for state, d in elec_data.items():
        channels.append(
            {
                "name": f"Electronic: {state}",
                "category": "electronic",
                "E_j": d["Ej"],
                "compute_K": lambda Te, _s=d["sigma"], _Ej=d["Ej"], _E=d["E"]: (
                    compute_K_constant_Ej(_E, _s, _Ej, Te)
                ),
            }
        )

    # Dissociative attachment, the electron energy is lost
    for species, data in _n2o.N2O_NEG_ION_CHANNELS.items():
        channels.append(
            {
                "name": f"Attachment: {species}",
                "category": "attachment",
                "E_j": "E_e (captured)",
                "compute_K": lambda Te, _E=data["energies"], _s=data["sigma"]: (
                    compute_K_Ee_captured(_E, _s, Te)
                ),
            }
        )

    return channels, ionisation_channels


# H2O, Petro and Sedwick (2017) Table 1


def build_h2o_channels_full():
    """H2O channel list of Petro and Sedwick (2017) Table 1: rotation (4),
    vibration (2), elastic, attachment (3), dissociation (4), electronic
    excitation (2), radiation (7) and H2O+ ionisation.

    Returns (channels, ionisation_channels).
    """
    channels = []
    ionisation_channels = []
    Ej = _h2o.H2O_Ej

    # Ionisation: H2O+ only
    h2o_plus = _h2o.H2O_ION_CHANNELS["H₂O⁺"]
    channels.append(
        {
            "name": "Ionisation: H₂O⁺",
            "category": "ionisation",
            "E_j": Ej["ionisation"],
            "compute_K": lambda Te: compute_K_constant_Ej(
                _h2o.H2O_ION_ENERGIES, h2o_plus["sigma"], Ej["ionisation"], Te
            ),
        }
    )
    ionisation_channels.append(
        {
            "name": "Ionisation: H₂O⁺",
            "sigma": h2o_plus["sigma"],
            "E_grid": _h2o.H2O_ION_ENERGIES,
        }
    )

    channels.append(
        {
            "name": "Elastic",
            "category": "elastic",
            "E_j": "2*E_e*(m_e/m_H2O)",
            "compute_K": lambda Te: compute_K_elastic(
                _h2o.H2O_ELASTIC_ENERGIES, _h2o.H2O_ELASTIC_ICS, _h2o.M_H2O, Te
            ),
        }
    )

    # Rotation 0->0 has the elastic energy loss
    channels.append(
        {
            "name": "Rotation: 0→0 (elastic)",
            "category": "rotation",
            "E_j": "2*E_e*(m_e/m_H2O)",
            "compute_K": lambda Te: compute_K_elastic(
                _h2o.H2O_ROT_ENERGIES, _h2o.H2O_ROT_0to0, _h2o.M_H2O, Te
            ),
        }
    )

    rot_data = [
        ("Rotation: 0→1", _h2o.H2O_ROT_0to1, Ej["rot_0to1"]),
        ("Rotation: 0→2", _h2o.H2O_ROT_0to2, Ej["rot_0to2"]),
        ("Rotation: 0→3", _h2o.H2O_ROT_0to3, Ej["rot_0to3"]),
    ]
    for name, sigma, ej in rot_data:
        channels.append(
            {
                "name": name,
                "category": "rotation",
                "E_j": ej,
                "compute_K": lambda Te, _s=sigma, _Ej=ej: compute_K_constant_Ej(
                    _h2o.H2O_ROT_ENERGIES, _s, _Ej, Te
                ),
            }
        )

    channels.append(
        {
            "name": "Vibration: stretch (100+001)",
            "category": "vibration",
            "E_j": Ej["vib_stretch"],
            "compute_K": lambda Te: compute_K_constant_Ej(
                _h2o.H2O_VIB_STRETCH_E,
                _h2o.H2O_VIB_STRETCH_SIGMA,
                Ej["vib_stretch"],
                Te,
            ),
        }
    )
    channels.append(
        {
            "name": "Vibration: bend (010)",
            "category": "vibration",
            "E_j": Ej["vib_bend"],
            "compute_K": lambda Te: compute_K_constant_Ej(
                _h2o.H2O_VIB_BEND_E, _h2o.H2O_VIB_BEND_SIGMA, Ej["vib_bend"], Te
            ),
        }
    )

    for species, data in _h2o.H2O_NEG_ION_CHANNELS.items():
        channels.append(
            {
                "name": f"Attachment: {species}",
                "category": "attachment",
                "E_j": "E_e (captured)",
                "compute_K": lambda Te, _E=data["energies"], _s=data["sigma"]: (
                    compute_K_Ee_captured(_E, _s, Te)
                ),
            }
        )

    # Dissociation into neutrals; OH(X)+H in 1e-16 cm^2
    channels.append(
        {
            "name": "Dissociation: OH(X)+H",
            "category": "dissociation",
            "E_j": Ej["diss_OHX_H"],
            "compute_K": lambda Te: compute_K_constant_Ej(
                _h2o.H2O_DISS_OHX_E, _h2o.H2O_DISS_OHX_SIGMA, Ej["diss_OHX_H"], Te
            ),
        }
    )
    # The three O + H2 channels use the O(1S) cross section (Petro and
    # Sedwick 2017, Table 1 footnote a); 1e-18 -> 1e-16 cm^2.
    O1S_sigma_converted = _h2o.H2O_DISS_O1S_SIGMA * 0.01
    for diss_label, ej_key in [
        ("Dissociation: O(³P)+H₂", "diss_O3P_H2"),
        ("Dissociation: O(¹D)+H₂", "diss_O1D_H2"),
        ("Dissociation: O(¹S)+H₂", "diss_O1S_H2"),
    ]:
        ej = Ej[ej_key]
        channels.append(
            {
                "name": diss_label,
                "category": "dissociation",
                "E_j": ej,
                "compute_K": lambda Te, _s=O1S_sigma_converted, _Ej=ej: (
                    compute_K_constant_Ej(_h2o.H2O_DISS_O1S_E, _s, _Ej, Te)
                ),
            }
        )

    # Electronic excitation, 1e-19 -> 1e-16 cm^2
    channels.append(
        {
            "name": "Electronic: ã³B₁",
            "category": "electronic",
            "E_j": Ej["elec_a3B1"],
            "compute_K": lambda Te: compute_K_constant_Ej(
                _h2o.H2O_ELEC_ENERGIES, _h2o.H2O_SIGMA_a3B1 * 0.001, Ej["elec_a3B1"], Te
            ),
        }
    )
    channels.append(
        {
            "name": "Electronic: Ã¹B₁",
            "category": "electronic",
            "E_j": Ej["elec_A1B1"],
            "compute_K": lambda Te: compute_K_constant_Ej(
                _h2o.H2O_ELEC_ENERGIES, _h2o.H2O_SIGMA_B1A1 * 0.001, Ej["elec_A1B1"], Te
            ),
        }
    )

    # Radiation; 1e-18 cm^2 (* 0.01) or 1e-19 cm^2 (* 0.001)
    channels.append(
        {
            "name": "Radiation: OH(A−X)",
            "category": "radiation",
            "E_j": Ej["rad_OHAX"],
            "compute_K": lambda Te: compute_K_constant_Ej(
                _h2o.H2O_RAD_OHAX_E, _h2o.H2O_RAD_OHAX_SIGMA * 0.01, Ej["rad_OHAX"], Te
            ),
        }
    )
    channels.append(
        {
            "name": "Radiation: H Lyman α",
            "category": "radiation",
            "E_j": Ej["rad_Lyman_a"],
            "compute_K": lambda Te: compute_K_constant_Ej(
                _h2o.H2O_RAD_LYMAN_E,
                _h2o.H2O_RAD_LYMAN_SIGMA * 0.01,
                Ej["rad_Lyman_a"],
                Te,
            ),
        }
    )
    channels.append(
        {
            "name": "Radiation: H Balmer α",
            "category": "radiation",
            "E_j": Ej["rad_Balmer_a"],
            "compute_K": lambda Te: compute_K_constant_Ej(
                _h2o.H2O_RAD_BALMER_A_E,
                _h2o.H2O_RAD_BALMER_A_SIGMA * 0.01,
                Ej["rad_Balmer_a"],
                Te,
            ),
        }
    )
    channels.append(
        {
            "name": "Radiation: H Balmer β",
            "category": "radiation",
            "E_j": Ej["rad_Balmer_b"],
            "compute_K": lambda Te: compute_K_constant_Ej(
                _h2o.H2O_RAD_T20_E,
                _h2o.H2O_RAD_BALMER_B_SIGMA * 0.001,
                Ej["rad_Balmer_b"],
                Te,
            ),
        }
    )
    channels.append(
        {
            "name": "Radiation: O 844.7 nm",
            "category": "radiation",
            "E_j": Ej["rad_O8447"],
            "compute_K": lambda Te: compute_K_constant_Ej(
                _h2o.H2O_RAD_T20_E,
                _h2o.H2O_RAD_O8447_SIGMA * 0.001,
                Ej["rad_O8447"],
                Te,
            ),
        }
    )
    channels.append(
        {
            "name": "Radiation: O 777.4 nm",
            "category": "radiation",
            "E_j": Ej["rad_O7774"],
            "compute_K": lambda Te: compute_K_constant_Ej(
                _h2o.H2O_RAD_T20_E,
                _h2o.H2O_RAD_O7774_SIGMA * 0.001,
                Ej["rad_O7774"],
                Te,
            ),
        }
    )
    channels.append(
        {
            "name": "Radiation: O 130.4 nm",
            "category": "radiation",
            "E_j": Ej["rad_O1304"],
            "compute_K": lambda Te: compute_K_constant_Ej(
                _h2o.H2O_RAD_LYMAN_E,
                _h2o.H2O_RAD_O1304_SIGMA * 0.001,
                Ej["rad_O1304"],
                Te,
            ),
        }
    )

    return channels, ionisation_channels


# Noble gases, from LXCat files


def build_noble_gas_channels_from_lxcat(
    lxcat_path: str,
    m_gas_kg: float,
    species_label: str,
):
    """Channel list of a noble gas from an LXCat file.

    Returns (channels, ionisation_channels).
    """
    blocks = parse_lxcat_file(lxcat_path)
    return build_noble_gas_channels_from_blocks(blocks, m_gas_kg, species_label)


def build_noble_gas_channels_from_blocks(blocks, m_gas_kg, species_label):
    """Channel list of a noble gas from parsed LXCat blocks."""
    channels = []
    ionisation_channels = []

    for b in blocks:
        bt = b["type"]
        E_grid = b["E"]
        sigma_1e16 = sigma_m2_to_1e16cm2(b["sigma"])

        if bt == "ELASTIC":
            channels.append(
                {
                    "name": f"{species_label} Elastic (MT)",
                    "category": "elastic",
                    "E_j": f"2*E_e*m_e/m_{species_label}",
                    "compute_K": (
                        lambda E_g=E_grid, s_g=sigma_1e16, m=m_gas_kg: (
                            lambda Te: compute_K_elastic(E_g, s_g, m, Te)
                        )
                    )(),
                }
            )

        elif bt == "EXCITATION":
            ej = b["threshold"]
            state_label = b["comment"].rstrip(".") if b["comment"] else f"{ej:.3f} eV"
            state_label = re.sub(r"^Excitation\s+", "", state_label)
            channels.append(
                {
                    "name": f"{species_label} Exc: {state_label} ({ej:.2f} eV)",
                    "category": "electronic",
                    "E_j": ej,
                    "compute_K": (
                        lambda E_g=E_grid, s_g=sigma_1e16, _ej=ej: (
                            lambda Te: compute_K_constant_Ej(E_g, s_g, _ej, Te)
                        )
                    )(),
                }
            )

        elif bt == "IONIZATION":
            ej = b["threshold"]
            channels.append(
                {
                    "name": f"{species_label}+ ({ej:.2f} eV)",
                    "category": "ionisation",
                    "E_j": ej,
                    "compute_K": (
                        lambda E_g=E_grid, s_g=sigma_1e16, _ej=ej: (
                            lambda Te: compute_K_constant_Ej(E_g, s_g, _ej, Te)
                        )
                    )(),
                }
            )
            ionisation_channels.append(
                {
                    "name": f"{species_label}+",
                    "sigma": sigma_1e16,
                    "E_grid": E_grid,
                }
            )

    return channels, ionisation_channels


# Attach the channel builders to the propellant registry
from . import propellants  # noqa: E402

propellants.attach_channel_builders()
