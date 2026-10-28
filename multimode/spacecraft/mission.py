"""
The two missions as impulsive-equivalent delta-v budgets per manoeuvre
class: LEO (station keeping, collision avoidance and its return, disposal)
and GEO (transfer from GTO, station keeping, graveyard).
"""

# Each class carries an impulsive delta-v and a loss ratio applied to the
# share flown electrically. The assignment rank sets the order in which the
# electric mode takes classes; the execution rank sets the order of the
# burns in the rocket-equation chain.

from dataclasses import dataclass
from typing import Optional, Tuple

from .parameters import value, YEAR_S


@dataclass(frozen=True)
class ManoeuvreClass:
    """One manoeuvre class of a mission budget.

    dv_impulsive    : impulsive-equivalent delta-v [m/s].
    ep_loss_ratio   : electric delta-v = dv_impulsive x ratio.
    divisible       : may be split between the two modes.
    assignment_rank : order in which the electric mode takes classes
                      (0 first).
    execution_rank  : order of the burns (0 first, at the highest mass).
    window_s        : schedule window [s]; None for collision avoidance
                      (the warning time) and the return burn.
    dv_per_event / n_events : per-event delta-v and number of events.
    w_duty          : fraction of time the electric mode can thrust.
    """

    key: str
    kind: str  # 'sk' | 'ca' | 'ca_restore' |
    # 'transfer' | 'eol'
    dv_impulsive: float
    ep_loss_ratio: float
    divisible: bool
    assignment_rank: int
    execution_rank: int
    window_s: Optional[float]
    dv_per_event: float
    n_events: int
    w_duty: float


@dataclass(frozen=True)
class Mission:
    """A mission; t_warn_s is the collision-avoidance warning time [s]
    (None without a collision-avoidance class)."""

    key: str
    display: str
    bus_mass: float  # [kg]
    bus_class: int  # 50 | 150 | 300 (power-margin row)
    duration_s: float
    classes: Tuple[ManoeuvreClass, ...]
    t_warn_s: Optional[float] = None

    @property
    def total_dv_impulsive(self):
        return sum(c.dv_impulsive for c in self.classes)

    def get_class(self, key):
        for c in self.classes:
            if c.key == key:
                return c
        raise KeyError(f"Mission {self.key!r} has no class {key!r}")


def ca_dv_per_event(t_warn_s):
    """Collision-avoidance delta-v per event, budget_1d x (1 d / t_warn):
    miss distance scales with delta-v times lead time."""
    return value("dv_ca_budget_1d") * value("t_warn_ref_ca") / t_warn_s


def leo_mission(t_warn_s=None):
    """LEO mission: 500 km, 150 kg bus, 5 years; station keeping,
    collision avoidance and its return, disposal.

    t_warn_s : warning time [s] (default: the register value, 1 d).
    """
    if t_warn_s is None:
        t_warn_s = value("t_warn")
    life = value("mission_life_leo") * YEAR_S
    dv_sk = value("dv_sk_leo_rate") * value("mission_life_leo")
    w = value("w_duty_leo")
    dv_ca_ev = ca_dv_per_event(t_warn_s)
    n_ca = int(value("n_ca_events"))
    dv_rs_ev = value("dv_ca_restore")
    classes = (
        ManoeuvreClass(
            key="sk",
            kind="sk",
            dv_impulsive=dv_sk,
            ep_loss_ratio=1.0,
            divisible=True,
            assignment_rank=1,
            execution_rank=1,
            window_s=life,
            dv_per_event=dv_sk,
            n_events=1,
            w_duty=w,
        ),
        ManoeuvreClass(
            key="ca",
            kind="ca",
            dv_impulsive=dv_ca_ev * n_ca,
            ep_loss_ratio=1.0,
            divisible=False,
            assignment_rank=3,
            execution_rank=2,
            window_s=None,
            dv_per_event=dv_ca_ev,
            n_events=n_ca,
            w_duty=w,
        ),
        # Return after each avoidance burn: not time-critical, so no window
        # and no feasibility verdict
        ManoeuvreClass(
            key="ca_restore",
            kind="ca_restore",
            dv_impulsive=dv_rs_ev * n_ca,
            ep_loss_ratio=1.0,
            divisible=True,
            assignment_rank=2,
            execution_rank=3,
            window_s=None,
            dv_per_event=dv_rs_ev,
            n_events=n_ca,
            w_duty=w,
        ),
        ManoeuvreClass(
            key="eol",
            kind="eol",
            dv_impulsive=value("dv_eol_leo"),
            ep_loss_ratio=value("eol_leo_ep_loss_ratio"),
            divisible=False,
            assignment_rank=0,
            execution_rank=4,
            window_s=value("eol_leo_window") * YEAR_S,
            dv_per_event=value("dv_eol_leo"),
            n_events=1,
            w_duty=w,
        ),
    )
    return Mission(
        key="leo",
        display="LEO mission — 500 km, SK + CA + disposal (5 yr)",
        bus_mass=value("bus_mass_leo"),
        bus_class=int(value("bus_mass_leo")),
        duration_s=life,
        classes=classes,
        t_warn_s=t_warn_s,
    )


def geo_mission(t_transfer_max_yr=None):
    """GEO mission: transfer from GTO, 15 years of station keeping,
    graveyard; 300 kg bus.

    t_transfer_max_yr : transfer window [yr] (default: the register value);
    it sizes the electric unit count and sets the transfer verdict.
    """
    if t_transfer_max_yr is None:
        t_transfer_max_yr = value("t_transfer_max")
    life = value("mission_life_geo") * YEAR_S
    dv_sk = value("dv_geo_sk_rate") * value("mission_life_geo")
    classes = (
        ManoeuvreClass(
            key="transfer",
            kind="transfer",
            dv_impulsive=value("dv_transfer_chem"),
            # electric transfer budget = chemical budget / planning efficiency
            ep_loss_ratio=1.0 / value("eta_planning_transfer"),
            divisible=True,
            assignment_rank=2,
            execution_rank=0,
            window_s=t_transfer_max_yr * YEAR_S,
            dv_per_event=value("dv_transfer_chem"),
            n_events=1,
            w_duty=value("w_duty_transfer"),
        ),
        ManoeuvreClass(
            key="sk",
            kind="sk",
            dv_impulsive=dv_sk,
            ep_loss_ratio=1.0,
            divisible=True,
            assignment_rank=1,
            execution_rank=1,
            window_s=life,
            dv_per_event=dv_sk,
            n_events=1,
            w_duty=value("w_duty_geo"),
        ),
        ManoeuvreClass(
            key="eol",
            kind="eol",
            dv_impulsive=value("dv_eol_geo"),
            ep_loss_ratio=1.0,
            divisible=False,
            assignment_rank=0,
            execution_rank=2,
            window_s=value("eol_geo_window") * YEAR_S,
            dv_per_event=value("dv_eol_geo"),
            n_events=1,
            w_duty=value("w_duty_geo"),
        ),
    )
    return Mission(
        key="geo",
        display="GEO mission — GTO to GEO transfer, SK, graveyard (15 yr)",
        bus_mass=value("bus_mass_geo"),
        bus_class=int(value("bus_mass_geo")),
        duration_s=life,
        classes=classes,
    )
