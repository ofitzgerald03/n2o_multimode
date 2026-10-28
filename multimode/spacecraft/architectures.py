"""
Thruster units and the five architectures: single-electric, bipropellant
multimode, hybrid, single-chemical and monopropellant multimode.
"""

# Each unit is a fixed head with given mass, Isp, thrust, mass flow and DC
# power; more thrust means more units. Electric DC power includes the
# conversion chain: DC-to-RF for the electrodeless unit, PPU for the Hall.

from dataclasses import dataclass
from typing import Optional, Tuple

from ..constants import G0
from .parameters import value


@dataclass(frozen=True)
class ThrusterUnit:
    """A thruster unit.

    p_dc        : DC bus power per unit [W], conversion chain included.
    electronics : 'rf' (DC-to-RF chain), 'ppu' or None (chemical).
    """

    key: str
    head_mass: float  # [kg]
    isp: float  # [s]
    thrust: float  # [N] per unit
    mdot: float  # [kg/s] per unit
    p_dc: float  # [W] per unit at the DC bus
    electronics: Optional[str]


def ecr_unit(eta_u=None):
    """N2O electrodeless unit. eta_u None: the register values (eta_u =
    0.61). Otherwise each quantity is the register value times the thruster
    model's ratio op(eta_u) / op(0.61): thrust and absorbed power held,
    mass flow ~ 1/eta_u, Isp ~ eta_u.
    """
    if eta_u is None:
        return ThrusterUnit(
            key="ecr",
            head_mass=value("m_head_ecr"),
            isp=value("isp_ecr"),
            thrust=value("t_unit_ecr"),
            mdot=value("mdot_unit_ecr"),
            p_dc=value("p_abs_ecr") / value("eta_dcrf"),  # absorbed to DC
            electronics="rf",
        )
    # Imported here: the thruster model runs only when eta_u is given
    from ..thruster.ecr_unit import ecr_unit_at, eta_u_parity

    op = ecr_unit_at(eta_u)
    op0 = ecr_unit_at(eta_u_parity())
    return ThrusterUnit(
        key="ecr",
        head_mass=value("m_head_ecr"),
        isp=value("isp_ecr") * (op.isp / op0.isp),
        thrust=value("t_unit_ecr") * (op.thrust / op0.thrust),
        mdot=value("mdot_unit_ecr") * (op.mdot / op0.mdot),
        p_dc=value("p_abs_ecr") * (op.p_abs / op0.p_abs) / value("eta_dcrf"),
        electronics="rf",
    )


def hall_unit():
    """Xe Hall unit, BHT-200 class."""
    isp = value("isp_hall")
    thrust = value("t_unit_hall")
    return ThrusterUnit(
        key="hall",
        head_mass=value("m_head_hall"),
        isp=isp,
        thrust=thrust,
        mdot=thrust / (isp * G0),
        p_dc=value("p_disch_hall") / value("eta_ppu"),  # discharge to DC
        electronics="ppu",
    )


def b20_unit():
    """Dawn Aerospace B20 N2O/propene unit."""
    isp = value("isp_chem")
    thrust = value("t_b20_nominal")
    return ThrusterUnit(
        key="b20",
        head_mass=value("m_head_b20"),
        isp=isp,
        thrust=thrust,
        mdot=thrust / (isp * G0),
        p_dc=0.0,  # valve and igniter power negligible
        electronics=None,
    )


def mono_unit():
    """N2O catalytic monopropellant unit; p_dc is the catalyst heater."""
    isp = value("isp_mono")
    thrust = value("t_mono_nominal")
    return ThrusterUnit(
        key="mono",
        head_mass=value("m_head_mono"),
        isp=isp,
        thrust=thrust,
        mdot=thrust / (isp * G0),
        p_dc=value("p_heater_mono"),
        electronics=None,
    )


@dataclass(frozen=True)
class Architecture:
    """One architecture.

    tanks     : installed tanks ('n2o', 'propene', 'xe'); an installed tank
                carries its minimum mass even when empty.
    chem_kind : 'biprop' (N2O/propene) | 'mono' (N2O) | None.
    ep_propellant : tank the electric unit draws from ('n2o' | 'xe').
    n_chem_units  : number of chemical units.
    """

    key: str
    display: str
    ep_unit: Optional[ThrusterUnit]
    ep_propellant: Optional[str]
    chem_unit: Optional[ThrusterUnit]
    chem_kind: Optional[str]
    tanks: Tuple[str, ...]
    n_chem_units: int = 1

    @property
    def has_ep(self):
        return self.ep_unit is not None

    @property
    def has_chem(self):
        return self.chem_unit is not None


_ARCH_CACHE = None

ARCHITECTURE_ORDER = (
    "single_elec",
    "multimode",
    "hybrid",
    "single_chem",
    "mono_multimode",
)
CONTENDERS = ARCHITECTURE_ORDER


def get_architectures():
    """Architectures built from the register, cached."""
    global _ARCH_CACHE
    if _ARCH_CACHE is not None:
        return _ARCH_CACHE
    ecr, hall, b20, mono = ecr_unit(), hall_unit(), b20_unit(), mono_unit()
    _ARCH_CACHE = {
        "single_elec": Architecture(
            key="single_elec",
            display="Single-electric (Xe Hall)",
            ep_unit=hall,
            ep_propellant="xe",
            chem_unit=None,
            chem_kind=None,
            tanks=("xe",),
            n_chem_units=0,
        ),
        "multimode": Architecture(
            key="multimode",
            display="Multimode (shared N2O: ECR + biprop)",
            ep_unit=ecr,
            ep_propellant="n2o",
            chem_unit=b20,
            chem_kind="biprop",
            tanks=("n2o", "propene"),
        ),
        "hybrid": Architecture(
            key="hybrid",
            display="Hybrid (Xe Hall + N2O/propene biprop)",
            ep_unit=hall,
            ep_propellant="xe",
            chem_unit=b20,
            chem_kind="biprop",
            tanks=("xe", "n2o", "propene"),
        ),
        "single_chem": Architecture(
            key="single_chem",
            display="Single-chemical (N2O/propene biprop)",
            ep_unit=None,
            ep_propellant=None,
            chem_unit=b20,
            chem_kind="biprop",
            tanks=("n2o", "propene"),
        ),
        "mono_multimode": Architecture(
            key="mono_multimode",
            display="Monoprop-multimode (shared N2O: ECR + monoprop)",
            ep_unit=ecr,
            ep_propellant="n2o",
            chem_unit=mono,
            chem_kind="mono",
            tanks=("n2o",),
        ),
    }
    return _ARCH_CACHE


def get_architecture(key):
    return get_architectures()[key]


def invalidate_cache():
    """Drop the cached architectures, so the next call rebuilds them from
    the current register (needed around parameters.override)."""
    global _ARCH_CACHE
    _ARCH_CACHE = None


def with_ecr_eta_u(arch_key, eta_u):
    """Copy of an electrodeless architecture with its electric unit at
    eta_u, for closure.close(..., arch=...)."""
    import dataclasses

    arch = get_architecture(arch_key)
    if arch.ep_unit is None or arch.ep_unit.key != "ecr":
        raise ValueError(
            f"{arch_key} does not carry the ECR head; eta_u^N2O does not apply"
        )
    return dataclasses.replace(arch, ep_unit=ecr_unit(eta_u))
