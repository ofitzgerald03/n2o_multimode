"""
Propellant registry: Xe and Ar (single species), N2O and H2O (multi-species).
"""

from dataclasses import dataclass
from typing import Callable, Optional

from ..constants import M_XE, M_AR
from .propellant_data import n2o as _n2o_data
from .propellant_data import h2o as _h2o_data


@dataclass
class Propellant:
    """One propellant.

    mass_kg : ion mass of a single-species propellant [kg].
    species_masses_kg : {ion: mass [kg]} of a multi-species propellant.
    channel_builder : callable returning (channels, ionisation_channels).
    ion_channels : {ion: {'sigma': array, 'mass_amu': float}} on ion_energies.
    neg_ion_channels : {ion: {'energies', 'sigma', 'mass_amu'}}, attachment.
    """

    label: str
    mass_kg: Optional[float] = None
    species_masses_kg: Optional[dict] = None
    channel_builder: Optional[Callable] = None
    ion_channels: Optional[dict] = None
    ion_energies: Optional[object] = None
    neg_ion_channels: Optional[dict] = None

    @property
    def is_multi_species(self) -> bool:
        return self.species_masses_kg is not None


XE = Propellant(label="Xe", mass_kg=M_XE)

AR = Propellant(label="Ar", mass_kg=M_AR)

N2O = Propellant(
    label="N2O",
    species_masses_kg=_n2o_data.N2O_MASSES_KG,
    ion_channels=_n2o_data.N2O_ION_CHANNELS,
    ion_energies=_n2o_data.N2O_ION_ENERGIES,
    neg_ion_channels=_n2o_data.N2O_NEG_ION_CHANNELS,
)

H2O = Propellant(
    label="H2O",
    species_masses_kg=_h2o_data.H2O_MASSES_KG,
    ion_channels=_h2o_data.H2O_ION_CHANNELS,
    ion_energies=_h2o_data.H2O_ION_ENERGIES,
    neg_ion_channels=_h2o_data.H2O_NEG_ION_CHANNELS,
)

ALL_PROPELLANTS = {p.label: p for p in (XE, AR, N2O, H2O)}


def attach_channel_builders():
    """Set each propellant's channel builder. Called from ionisation_cost
    at import, which avoids a circular import."""
    from .ionisation_cost import (
        build_n2o_channels,
        build_h2o_channels_full,
        build_noble_gas_channels_from_lxcat,
    )
    from ..config import PATHS

    N2O.channel_builder = lambda: build_n2o_channels(
        vib_dir=str(PATHS["cross_sections"] / "n2o")
    )
    H2O.channel_builder = lambda: build_h2o_channels_full()
    XE.channel_builder = lambda: build_noble_gas_channels_from_lxcat(
        str(PATHS["xe_lxcat"]),
        XE.mass_kg,
        XE.label,
    )
    AR.channel_builder = lambda: build_noble_gas_channels_from_lxcat(
        str(PATHS["ar_lxcat"]),
        AR.mass_kg,
        AR.label,
    )
