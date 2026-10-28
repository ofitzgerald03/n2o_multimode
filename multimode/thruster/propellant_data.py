"""
Electron-impact cross sections of N2O and H2O, read from the CSV files in
data/cross_sections/ at import time. One namespace per species (`n2o`,
`h2o`). Energies in eV; cross sections in 1e-16 cm^2 unless the attribute
name says otherwise.
"""

from types import SimpleNamespace

import numpy as np
import pandas as pd

from ..constants import AMU
from ..config import PATHS


_XS_DIR = PATHS["cross_sections"]


def _read(species: str, filename: str) -> pd.DataFrame:
    """Read one cross-section CSV ('#' lines are header comments)."""
    return pd.read_csv(
        _XS_DIR / species / filename, comment="#", float_precision="round_trip"
    )


def _cols(df: pd.DataFrame, *names):
    """Named columns as float arrays (empty cells -> NaN)."""
    out = tuple(df[n].to_numpy(dtype=float) for n in names)
    return out[0] if len(out) == 1 else out


def _read_state_table(species: str, filename: str):
    """(energies, sigma) of a two-column state table."""
    df = _read(species, filename)
    return _cols(df, "E_eV", "sigma_1e16_cm2")


def _process_elec_state(E_grid, raw, threshold_eV):
    """Electronic-excitation cross section: zero below threshold, linear
    interpolation across gaps at and above it."""
    sigma = raw.copy()
    sigma[E_grid < threshold_eV] = 0.0
    above = E_grid >= threshold_eV
    known = above & ~np.isnan(sigma)
    if np.sum(known) >= 2:
        sigma[above] = np.interp(E_grid[above], E_grid[known], sigma[known])
    sigma = np.where(np.isfinite(sigma), sigma, 0.0)
    return sigma


def _ion_threshold_onset(E_tab, raw, Ej, step_eV=0.01):
    """Partial-ionisation cross sections with a linear threshold onset.

    `raw` maps channel -> sigma tabulated on E_tab; `Ej` maps channel ->
    appearance energy. Each channel rises linearly from zero at E_j to its
    first non-zero tabulated value, in steps of `step_eV` on a grid shared
    by all channels; between tabulated points a channel takes the log-log
    interpolant of its tabulated values. Returns (E_grid, {channel: sigma}).
    """
    first = {}
    ramps = []
    for sp, s in raw.items():
        i1 = int(np.flatnonzero(s > 0)[0])
        first[sp] = (E_tab[i1], s[i1])
        n = max(1, int(np.ceil((E_tab[i1] - Ej[sp]) / step_eV)))
        ramps.append(np.round(np.linspace(Ej[sp], E_tab[i1], n + 1)[:-1], 9))
    E = np.unique(np.concatenate([E_tab] + ramps))
    out = {}
    for sp, s in raw.items():
        E1, s1 = first[sp]
        tab = E_tab >= E1
        sig = np.exp(np.interp(np.log(E), np.log(E_tab[tab]), np.log(s[tab])))
        on = (E >= Ej[sp]) & (E < E1)
        sig[on] = s1 * (E[on] - Ej[sp]) / (E1 - Ej[sp])
        sig[E < Ej[sp]] = 0.0
        sig[np.isin(E, E_tab[tab])] = s[tab]
        out[sp] = sig
    return E, out


# N2O: Song et al. (2019) Tables XII, XIV, XV and Fig. 28


def _load_n2o() -> SimpleNamespace:
    ns = SimpleNamespace()

    # Ionisation appearance energies E_j [eV], NIST Chemistry WebBook
    ns.N2O_ION_Ej = {
        "N₂O⁺": 12.889,  # ionisation energy (Kimura, Katsumata et al.)
        "NO⁺": 15.3,  # electron impact (Olivier, Locht and Momigny 1982)
        "N₂⁺": 17.3,  # electron impact (Olivier, Locht and Momigny 1982)
        "N⁺": 19.49,  # photoelectron (Berkowitz and Eland 1977)
        "O⁺": 15.31,  # photoionisation (Dibeler et al. 1967)
    }

    # Partial ionisation (Table XIV), with threshold onset
    df = _read("n2o", "ionisation.csv")
    _keys = (
        ("N₂O⁺", "N2O+"),
        ("NO⁺", "NO+"),
        ("N₂⁺", "N2+"),
        ("N⁺", "N+"),
        ("O⁺", "O+"),
    )
    ns.N2O_ION_ENERGIES, _sig = _ion_threshold_onset(
        _cols(df, "E_eV"), {sp: _cols(df, col) for sp, col in _keys}, ns.N2O_ION_Ej
    )
    ns.N2O_SIGMA_N2Oplus = _sig["N₂O⁺"]
    ns.N2O_SIGMA_NOplus = _sig["NO⁺"]
    ns.N2O_SIGMA_N2plus = _sig["N₂⁺"]
    ns.N2O_SIGMA_Nplus = _sig["N⁺"]
    ns.N2O_SIGMA_Oplus = _sig["O⁺"]

    ns.N2O_ION_CHANNELS = {
        "N₂O⁺": {"sigma": ns.N2O_SIGMA_N2Oplus, "mass_amu": 44.0},
        "NO⁺": {"sigma": ns.N2O_SIGMA_NOplus, "mass_amu": 30.0},
        "N₂⁺": {"sigma": ns.N2O_SIGMA_N2plus, "mass_amu": 28.0},
        "N⁺": {"sigma": ns.N2O_SIGMA_Nplus, "mass_amu": 14.0},
        "O⁺": {"sigma": ns.N2O_SIGMA_Oplus, "mass_amu": 16.0},
    }
    ns.N2O_MASSES_AMU = {sp: d["mass_amu"] for sp, d in ns.N2O_ION_CHANNELS.items()}
    ns.N2O_MASSES_KG = {sp: m * AMU for sp, m in ns.N2O_MASSES_AMU.items()}

    # Dissociative attachment (Table XV)
    df = _read("n2o", "attachment_Ominus.csv")
    ns.N2O_DEA_Ominus_ENERGIES, ns.N2O_DEA_Ominus_SIGMA = _cols(
        df, "E_eV", "sigma_1e18_cm2"
    )
    ns.N2O_NEG_ION_CHANNELS = {
        "O⁻": {
            "energies": ns.N2O_DEA_Ominus_ENERGIES,
            "sigma": ns.N2O_DEA_Ominus_SIGMA * 0.01,  # 1e-18 -> 1e-16 cm^2
            "mass_amu": 16.0,
        },
    }

    # Elastic momentum transfer (Table XII)
    df = _read("n2o", "elastic_momentum_transfer.csv")
    ns.N2O_ELASTIC_MT_ENERGIES, ns.N2O_ELASTIC_MTCS = _cols(
        df, "E_eV", "sigma_m_1e16_cm2"
    )

    # Electronic excitation: swarm-derived components (Fig. 28)
    # Three channels, each on its own grid; E_j = threshold.
    ns.N2O_ELEC_SWARM_SPLIT = {}
    for key, fname, thr in (
        ("Qex1", "electronic_swarm_Qex1_5p6eV.csv", 5.6),
        ("Qex2", "electronic_swarm_Qex2_8p5eV.csv", 8.5),
        ("Qex3", "electronic_swarm_Qex3_9p6eV.csv", 9.6),
    ):
        E_sw, raw_sw = _read_state_table("n2o", fname)
        ns.N2O_ELEC_SWARM_SPLIT[key] = {
            "energies": E_sw,
            "sigma": _process_elec_state(E_sw, raw_sw, thr),
            "threshold": thr,
        }

    # Vibrational E_j [eV], Bosi (2014) Table 4
    ns.N2O_VIB_Ej = {
        "nu2": 0.073,  # bend
        "nu3_2nu2": 0.159,  # NO stretch / bend overtone
        "nu1_2nu3": 0.276,  # NN stretch / NO stretch overtone
    }

    ns.M_N2O_AMU = 44.013
    ns.M_N2O = ns.M_N2O_AMU * AMU
    return ns


# H2O: Itikawa and Mason (2005); Thorn et al. (2007) Table 4.
# Energy losses E_j from Petro and Sedwick (2017) Table 1.


def _load_h2o() -> SimpleNamespace:
    ns = SimpleNamespace()

    # Partial ionisation (Itikawa and Mason 2005, Table 11)
    df = _read("h2o", "ionisation.csv")
    ns.H2O_ION_ENERGIES = _cols(df, "E_eV")
    ns.H2O_SIGMA_H2Oplus = _cols(df, "H2O+")
    ns.H2O_SIGMA_OHplus = _cols(df, "OH+")
    ns.H2O_SIGMA_Hplus = _cols(df, "H+")
    ns.H2O_SIGMA_Oplus = _cols(df, "O+")
    ns.H2O_SIGMA_H2plus = _cols(df, "H2+")
    ns.H2O_SIGMA_Oplusplus = _cols(df, "O++")

    ns.H2O_ION_CHANNELS = {
        "H₂O⁺": {"sigma": ns.H2O_SIGMA_H2Oplus, "mass_amu": 18.0},
        "OH⁺": {"sigma": ns.H2O_SIGMA_OHplus, "mass_amu": 17.0},
        "H⁺": {"sigma": ns.H2O_SIGMA_Hplus, "mass_amu": 1.0},
        "O⁺": {"sigma": ns.H2O_SIGMA_Oplus, "mass_amu": 16.0},
        "H₂⁺": {"sigma": ns.H2O_SIGMA_H2plus, "mass_amu": 2.0},
        "O⁺⁺": {"sigma": ns.H2O_SIGMA_Oplusplus, "mass_amu": 16.0},
    }
    ns.H2O_MASSES_AMU = {sp: d["mass_amu"] for sp, d in ns.H2O_ION_CHANNELS.items()}
    ns.H2O_MASSES_KG = {sp: m * AMU for sp, m in ns.H2O_MASSES_AMU.items()}

    # Dissociative attachment (Itikawa and Mason 2005, Tables 13-15)
    dea = {}
    for key, fname in (
        ("Hminus", "attachment_Hminus.csv"),
        ("Ominus", "attachment_Ominus.csv"),
        ("OHminus", "attachment_OHminus.csv"),
    ):
        df = _read("h2o", fname)
        dea[key] = _cols(df, "E_eV", "sigma_1e18_cm2")
    ns.H2O_DEA_Hminus_ENERGIES, ns.H2O_DEA_Hminus_SIGMA = dea["Hminus"]
    ns.H2O_DEA_Ominus_ENERGIES, ns.H2O_DEA_Ominus_SIGMA = dea["Ominus"]
    ns.H2O_DEA_OHminus_ENERGIES, ns.H2O_DEA_OHminus_SIGMA = dea["OHminus"]
    ns.H2O_NEG_ION_CHANNELS = {
        "H⁻": {
            "energies": ns.H2O_DEA_Hminus_ENERGIES,
            "sigma": ns.H2O_DEA_Hminus_SIGMA * 0.01,  # 1e-18 -> 1e-16
            "mass_amu": 1.0,
        },
        "O⁻": {
            "energies": ns.H2O_DEA_Ominus_ENERGIES,
            "sigma": ns.H2O_DEA_Ominus_SIGMA * 0.01,
            "mass_amu": 16.0,
        },
        "OH⁻": {
            "energies": ns.H2O_DEA_OHminus_ENERGIES,
            "sigma": ns.H2O_DEA_OHminus_SIGMA * 0.01,
            "mass_amu": 17.0,
        },
    }

    # Elastic scattering (Itikawa and Mason 2005, Table 4)
    df = _read("h2o", "elastic.csv")
    ns.H2O_ELASTIC_ENERGIES, ns.H2O_ELASTIC_ICS = _cols(df, "E_eV", "sigma_1e16_cm2")

    # Rotational excitation (Itikawa and Mason 2005, Tables 7-8)
    # Empty cells below threshold are zero.
    df = _read("h2o", "rotation.csv")
    ns.H2O_ROT_ENERGIES = _cols(df, "E_eV")
    ns.H2O_ROT_0to0 = _cols(df, "rot_0to0")
    for trans in ("0to1", "0to2", "0to3"):
        raw = _cols(df, f"rot_{trans}")
        setattr(ns, f"H2O_ROT_{trans}", np.where(np.isfinite(raw), raw, 0.0))

    # Vibrational excitation (Itikawa and Mason 2005, Table 9)
    df = _read("h2o", "vibration_stretch.csv")
    ns.H2O_VIB_STRETCH_E, ns.H2O_VIB_STRETCH_SIGMA = _cols(df, "E_eV", "sigma_1e16_cm2")
    df = _read("h2o", "vibration_bend.csv")
    ns.H2O_VIB_BEND_E, ns.H2O_VIB_BEND_SIGMA = _cols(df, "E_eV", "sigma_1e16_cm2")

    # Electronic excitation (Thorn et al. 2007, Table 4; 1e-19 cm^2)
    df = _read("h2o", "electronic.csv")
    ns.H2O_ELEC_ENERGIES = _cols(df, "E_eV")
    ns.H2O_SIGMA_a3B1 = _cols(df, "a3B1")
    ns.H2O_SIGMA_B1A1 = _cols(df, "B1A1")

    # Dissociation into neutrals (Itikawa and Mason 2005, Tables 23-24)
    df = _read("h2o", "dissociation_OHX.csv")
    ns.H2O_DISS_OHX_E, ns.H2O_DISS_OHX_SIGMA = _cols(df, "E_eV", "sigma_1e16_cm2")
    df = _read("h2o", "dissociation_O1S.csv")
    ns.H2O_DISS_O1S_E, ns.H2O_DISS_O1S_SIGMA = _cols(df, "E_eV", "sigma_1e18_cm2")

    # Radiation (Itikawa and Mason 2005, Tables 17, 19, 20, 22)
    df = _read("h2o", "radiation_OHAX.csv")
    ns.H2O_RAD_OHAX_E, ns.H2O_RAD_OHAX_SIGMA = _cols(df, "E_eV", "sigma_1e18_cm2")
    df = _read("h2o", "radiation_Lyman_grid.csv")
    ns.H2O_RAD_LYMAN_E = _cols(df, "E_eV")
    ns.H2O_RAD_LYMAN_SIGMA = _cols(df, "Lyman_a_1e18_cm2")
    ns.H2O_RAD_O1304_SIGMA = _cols(df, "O1304_1e19_cm2")
    df = _read("h2o", "radiation_Balmer_a.csv")
    ns.H2O_RAD_BALMER_A_E, ns.H2O_RAD_BALMER_A_SIGMA = _cols(
        df, "E_eV", "sigma_1e18_cm2"
    )
    df = _read("h2o", "radiation_T20_grid.csv")
    ns.H2O_RAD_T20_E = _cols(df, "E_eV")
    ns.H2O_RAD_BALMER_B_SIGMA = _cols(df, "Balmer_b")
    ns.H2O_RAD_O8447_SIGMA = _cols(df, "O8447")
    ns.H2O_RAD_O7774_SIGMA = _cols(df, "O7774")

    # Energy loss per channel E_j [eV], Petro and Sedwick (2017) Table 1
    # The 5.1 eV dissociation energy is subtracted from the states that
    # dissociate (Petro and Sedwick 2017, p. 4).
    ns.H2O_Ej = {
        "rot_0to0": "2*E_e*(m_e/m_H2O)",  # elastic form
        "rot_0to1": 0.0049,
        "rot_0to2": 0.0128,
        "rot_0to3": 0.0265,
        "vib_stretch": 0.45,
        "vib_bend": 0.20,
        "diss_OHX_H": 5.1,
        "diss_O3P_H2": 5.0,
        "diss_O1D_H2": 7.0,
        "diss_O1S_H2": 9.2,
        "elec_a3B1": 7.14 - 5.1,
        "elec_A1B1": 7.49 - 5.1,
        "rad_OHAX": 9.15 - 5.1,
        "rad_Lyman_a": 15.3 - 5.1,
        "rad_Balmer_a": 17.2 - 5.1,
        "rad_Balmer_b": 18.6 - 5.1,
        "rad_O8447": 14.6,
        "rad_O7774": 15.7,
        "rad_O1304": 14.7,
        "ionisation": 12.6,
    }

    ns.M_H2O_AMU = 18.015
    ns.M_H2O = ns.M_H2O_AMU * AMU
    return ns


n2o = _load_n2o()
h2o = _load_h2o()
