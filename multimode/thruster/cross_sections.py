"""
Cross-section processing: Maxwellian rate constants, energy-loss rate
coefficients K_j, ion fractions, detachment cross sections and LXCat
parsing. Cross sections in 1e-16 cm^2, energies in eV.
"""

import re
from typing import Callable

import numpy as np
from scipy import interpolate

from ..constants import E_CHARGE, M_ELECTRON


def compute_rate_constant(
    energies_eV: np.ndarray,
    sigma_1e16cm2: np.ndarray,
    Te_eV: float,
) -> float:
    """Maxwellian rate constant kappa(Te) [m^3/s], Eq. (1) of Petro and
    Sedwick (2017):

        kappa = sqrt(8 e / (pi m_e Te^3)) * int Ee sigma(Ee) exp(-Ee/Te) dEe

    with sigma interpolated log-log over its non-zero points.
    """
    sigma_m2 = sigma_1e16cm2 * 1e-20  # 1e-16 cm^2 = 1e-20 m^2

    nonzero = sigma_m2 > 0
    if not np.any(nonzero):
        return 0.0

    e_valid = energies_eV[nonzero]
    s_valid = sigma_m2[nonzero]

    log_interp = interpolate.interp1d(
        np.log(e_valid),
        np.log(s_valid),
        kind="linear",
        bounds_error=False,
        fill_value=-np.inf,
    )

    def sigma_func(Ee):
        result = np.zeros_like(Ee, dtype=float)
        valid = (Ee >= e_valid[0]) & (Ee <= e_valid[-1])
        if np.any(valid):
            result[valid] = np.exp(log_interp(np.log(Ee[valid])))
        return result

    E_max = energies_eV[-1]
    E_grid = np.linspace(0, E_max, 10000)
    sigma_values = sigma_func(E_grid)
    integrand = E_grid * sigma_values * np.exp(-E_grid / Te_eV)
    integral_value = np.trapezoid(integrand, E_grid)

    prefactor = np.sqrt(8.0 * E_CHARGE / (np.pi * M_ELECTRON)) * Te_eV ** (-1.5)
    return prefactor * integral_value


def compute_rate_constants_vs_Te(
    energies_eV: np.ndarray,
    sigma_1e16cm2: np.ndarray,
    Te_array: np.ndarray,
) -> np.ndarray:
    """compute_rate_constant over a Te grid."""
    return np.array(
        [compute_rate_constant(energies_eV, sigma_1e16cm2, Te) for Te in Te_array]
    )


def compute_ion_fractions(kappa_dict: dict) -> dict:
    """Ion fractions with a common system Bohm velocity:
    n_i / n = kappa_i / sum_j kappa_j. Returns {species: f_i(Te)}."""
    total = sum(kappa_dict.values())
    return {
        sp: np.where(total > 0, kappa / total, 0.0) for sp, kappa in kappa_dict.items()
    }


def compute_ion_fractions_species_resolved(kappa_dict: dict, mass_dict: dict) -> dict:
    """Ion fractions with species-resolved Bohm velocities, Eq. (4) of Petro
    and Sedwick (2017): n_i / n = kappa_i sqrt(m_i) / sum_j kappa_j sqrt(m_j).
    Returns {species: f_i(Te)}."""
    weighted = {sp: kappa_dict[sp] * np.sqrt(mass_dict[sp]) for sp in kappa_dict}
    total = sum(weighted.values())
    return {sp: np.where(total > 0, w / total, 0.0) for sp, w in weighted.items()}


# Energy-loss rate coefficients K_j = E_j kappa_j, Eq. (8) of Petro and
# Sedwick (2017)


def compute_K_constant_Ej(
    E_grid: np.ndarray,
    sigma_1e16cm2: np.ndarray,
    E_loss_eV: float,
    Te_array: np.ndarray,
) -> np.ndarray:
    """K_j = E_j kappa_j for a constant E_j [eV m^3/s]."""
    kappa = compute_rate_constants_vs_Te(E_grid, sigma_1e16cm2, Te_array)
    return E_loss_eV * kappa


def _build_loglog_interpolant(E_grid, sigma_m2):
    """Log-log interpolant of sigma [m^2] against E [eV]; returns
    (sigma_func, E_max)."""
    nonzero = (sigma_m2 > 0) & (E_grid > 0)
    if np.sum(nonzero) < 2:
        return lambda E: np.zeros_like(E, dtype=float), float(np.max(E_grid))

    log_interp = interpolate.interp1d(
        np.log(E_grid[nonzero]),
        np.log(sigma_m2[nonzero]),
        kind="linear",
        bounds_error=False,
        fill_value=-np.inf,
    )
    E_max = float(np.max(E_grid))

    def sigma_func(E_fine):
        result = np.zeros_like(E_fine, dtype=float)
        valid = (E_fine >= E_grid[nonzero][0]) & (E_fine <= E_grid[nonzero][-1])
        if np.any(valid):
            result[valid] = np.exp(log_interp(np.log(E_fine[valid])))
        return result

    return sigma_func, E_max


def compute_K_elastic(
    E_grid: np.ndarray,
    sigma_1e16cm2: np.ndarray,
    m_gas_kg: float,
    Te_array: np.ndarray,
    n_points: int = 10000,
) -> np.ndarray:
    """K for elastic scattering, E_j = 2 Ee m_e / m_gas [eV m^3/s]:

    K = (2 m_e / m_gas) sqrt(8 e / (pi m_e)) Te^(-3/2)
        * int Ee^2 sigma(Ee) exp(-Ee/Te) dEe
    """
    sigma_m2 = np.asarray(sigma_1e16cm2) * 1e-20
    sigma_func, E_max = _build_loglog_interpolant(E_grid, sigma_m2)

    E_fine = np.linspace(1e-6, E_max, n_points)
    sigma_fine = sigma_func(E_fine)

    mass_ratio = M_ELECTRON / m_gas_kg
    base_prefactor = 2.0 * mass_ratio * np.sqrt(8.0 * E_CHARGE / (np.pi * M_ELECTRON))

    result = np.zeros(len(Te_array))
    for k, Te in enumerate(Te_array):
        integrand = E_fine**2 * sigma_fine * np.exp(-E_fine / Te)
        result[k] = base_prefactor * Te ** (-1.5) * np.trapezoid(integrand, E_fine)

    return result


def compute_K_Ee_captured(
    E_grid: np.ndarray,
    sigma_1e16cm2: np.ndarray,
    Te_array: np.ndarray,
    n_points: int = 10000,
) -> np.ndarray:
    """K for dissociative attachment, E_j = Ee [eV m^3/s]:

    K = sqrt(8 e / (pi m_e)) Te^(-3/2) int Ee^2 sigma(Ee) exp(-Ee/Te) dEe
    """
    sigma_m2 = np.asarray(sigma_1e16cm2) * 1e-20
    sigma_func, E_max = _build_loglog_interpolant(E_grid, sigma_m2)

    E_fine = np.linspace(1e-6, E_max, n_points)
    sigma_fine = sigma_func(E_fine)
    base_prefactor = np.sqrt(8.0 * E_CHARGE / (np.pi * M_ELECTRON))

    result = np.zeros(len(Te_array))
    for k, Te in enumerate(Te_array):
        integrand = E_fine**2 * sigma_fine * np.exp(-E_fine / Te)
        result[k] = base_prefactor * Te ** (-1.5) * np.trapezoid(integrand, E_fine)

    return result


# Electron-impact detachment


def compute_detach_rate_constant(
    sigma_func: Callable,
    Te_eV: float,
    E_max: float = 1000.0,
    n_points: int = 10000,
) -> float:
    """Maxwellian rate constant [m^3/s] of a callable sigma(E) in
    1e-16 cm^2."""
    E_grid = np.linspace(0, E_max, n_points)
    sigma_m2 = sigma_func(E_grid) * 1e-20
    integrand = E_grid * sigma_m2 * np.exp(-E_grid / Te_eV)
    integral_value = np.trapezoid(integrand, E_grid)
    prefactor = np.sqrt(8.0 * E_CHARGE / (np.pi * M_ELECTRON)) * Te_eV ** (-1.5)
    return prefactor * integral_value


def compute_detach_rate_constants_vs_Te(
    sigma_func: Callable,
    Te_array: np.ndarray,
) -> np.ndarray:
    return np.array([compute_detach_rate_constant(sigma_func, Te) for Te in Te_array])


def sigma_detach_Hminus(E):
    """H- detachment cross section [1e-16 cm^2], Janev et al. (2003)
    Eq. (58); threshold 0.754 eV."""
    E = np.asarray(E, dtype=float)
    sigma_cm2 = np.zeros_like(E)
    valid = E > 0.754
    Ev = E[valid]
    euler = np.e
    sigma_cm2[valid] = (
        (2.06 / Ev)
        * np.log(euler + 2.335e-4 * Ev)
        * np.exp(-13.75 / (Ev**0.868))
        * (1.0 - (0.754 / Ev) ** 2)
        * 1e-13
    )
    sigma_cm2 = np.maximum(sigma_cm2, 0.0)
    return sigma_cm2 * 1e16


def sigma_detach_Ominus(E):
    """O- detachment cross section [1e-16 cm^2], Vejby-Christensen et al.
    (1996) Eq. (16): sigma = p pi R^2 (1 - E_th/E), R = 8 a0, E_th = 3.4 eV,
    p = 0.12."""
    E = np.asarray(E, dtype=float)
    a0_cm = 5.29177210903e-9
    R = 8.0 * a0_cm
    Eth = 3.4
    p = 0.12
    sigma_cm2 = np.zeros_like(E, dtype=float)
    valid = E > Eth
    sigma_cm2[valid] = p * np.pi * R**2 * (1.0 - Eth / E[valid])
    return sigma_cm2 * 1e16


# LXCat-format parsing


def parse_lxcat_file(path: str) -> list:
    """Blocks of an LXCat text file, each a dict with 'type', 'species',
    'threshold', 'process', 'comment', 'E' [eV] and 'sigma' [m^2]."""
    with open(path, "r") as f:
        text = f.read()

    type_re = re.compile(
        r"^(ELASTIC|EFFECTIVE|EXCITATION|IONIZATION|ATTACHMENT)\s*$",
        re.MULTILINE,
    )
    matches = list(type_re.finditer(text))

    blocks = []
    for i, m in enumerate(matches):
        start = m.start()
        end = matches[i + 1].start() if i + 1 < len(matches) else len(text)
        blocks.append(_parse_lxcat_block(text[start:end]))
    return blocks


def _parse_lxcat_block(raw: str) -> dict:
    lines = [ln.rstrip("\r\n") for ln in raw.splitlines()]
    block = {
        "type": lines[0].strip(),
        "species": lines[1].strip(),
    }
    try:
        block["threshold"] = float(lines[2].strip())
    except ValueError:
        block["threshold"] = None

    block["process"] = ""
    block["comment"] = ""
    for ln in lines[3:]:
        if ln.startswith("PROCESS:"):
            block["process"] = ln[len("PROCESS:") :].strip()
        elif ln.startswith("COMMENT:"):
            block["comment"] = ln[len("COMMENT:") :].strip()

    dash_re = re.compile(r"^-{5,}\s*$")
    dash_idx = [i for i, ln in enumerate(lines) if dash_re.match(ln)]
    data_lines = lines[dash_idx[0] + 1 : dash_idx[1]] if len(dash_idx) >= 2 else []

    E_list, sigma_list = [], []
    for ln in data_lines:
        ln = ln.strip()
        if not ln:
            continue
        parts = ln.split()
        if len(parts) >= 2:
            try:
                E_list.append(float(parts[0]))
                sigma_list.append(float(parts[1]))
            except ValueError:
                continue

    block["E"] = np.array(E_list)
    block["sigma"] = np.array(sigma_list)
    return block


def sigma_m2_to_1e16cm2(sigma_m2):
    """m^2 -> 1e-16 cm^2."""
    return sigma_m2 / 1e-20


# Digitised vibrational cross sections


def load_vib_csv(filepath: str):
    """(energies [eV], sigma [1e-16 cm^2]) of a two-column CSV, sorted by
    energy with repeated energies removed."""
    data = np.loadtxt(filepath, delimiter=",", comments="#")
    energies = data[:, 0]
    sigma = data[:, 1]

    order = np.argsort(energies)
    energies = energies[order]
    sigma = sigma[order]

    _, unique_idx = np.unique(energies, return_index=True)
    return energies[unique_idx], sigma[unique_idx]
