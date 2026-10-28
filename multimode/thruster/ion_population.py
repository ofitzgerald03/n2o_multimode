"""
Positive-ion rate constants and fractions, attachment rate constants and
the negative-ion-to-neutral ratio of a propellant.
"""

import numpy as np

from .cross_sections import (
    compute_rate_constants_vs_Te,
    compute_ion_fractions,
    compute_detach_rate_constants_vs_Te,
    sigma_detach_Hminus,
    sigma_detach_Ominus,
)


def analyse_ion_populations(
    ion_channels,
    ion_energies,
    Te_array,
    neg_ion_channels=None,
):
    """Ion rate constants and fractions of one propellant.

    ion_channels : {species: {'sigma', 'mass_amu'}} on ion_energies [eV].
    neg_ion_channels : {species: {'energies', 'sigma', 'mass_amu'}},
        attachment cross sections in 1e-16 cm^2, each on its own grid.

    Returns a dict with 'Te', 'kappa_pos' {species: kappa_i [m^3/s]},
    'kappa_total', 'fractions' {species: f_i} and 'kappa_neg'
    {species: kappa_attach}.
    """
    kappa_pos = {}
    for sp, data in ion_channels.items():
        kappa_pos[sp] = compute_rate_constants_vs_Te(
            ion_energies, data["sigma"], Te_array
        )

    fractions = compute_ion_fractions(kappa_pos)
    kappa_total = sum(kappa_pos.values())

    kappa_neg = {}
    if neg_ion_channels:
        for sp, data in neg_ion_channels.items():
            kappa_neg[sp] = compute_rate_constants_vs_Te(
                data["energies"], data["sigma"], Te_array
            )

    return {
        "Te": Te_array,
        "kappa_pos": kappa_pos,
        "kappa_total": kappa_total,
        "fractions": fractions,
        "kappa_neg": kappa_neg,
    }


# Detachment cross sections by negative ion; no model for OH-
DETACH_MODELS = {
    "H⁻": sigma_detach_Hminus,
    "O⁻": sigma_detach_Ominus,
}


def compute_neg_ion_ratios(kappa_neg_dict, Te_array):
    """Negative-ion-to-neutral ratio, Eq. (3) of Petro and Sedwick (2017),
    summed over species: n_- / n_n = sum_j kappa_attach,j / kappa_detach,j.
    Species without a detachment model are skipped.

    Returns ({species: ratio(Te)}, total ratio(Te)).
    """
    individual_ratios = {}
    for sp, kappa_attach in kappa_neg_dict.items():
        if sp not in DETACH_MODELS:
            continue
        kappa_detach = compute_detach_rate_constants_vs_Te(
            DETACH_MODELS[sp],
            Te_array,
        )
        individual_ratios[sp] = np.where(
            kappa_detach > 0,
            kappa_attach / kappa_detach,
            0.0,
        )

    total_ratio = (
        sum(individual_ratios.values())
        if individual_ratios
        else np.zeros_like(Te_array)
    )
    return individual_ratios, total_ratio
