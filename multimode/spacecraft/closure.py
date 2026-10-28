"""
Propulsion wet-mass closure at fixed bus mass. For one architecture,
mission and allocation the wet mass is the fixed point of
"""

#     m_wet -> dry(m_wet) + loaded propellant(m_wet),
#
# with the propellant from the rocket equation, burned in execution order
# (transfer, station keeping, collision avoidance, return, disposal) from
# m0 = bus + m_wet; within a split class the chemical share burns first.
#
# Unit counts are integers resolved in the loop: one chemical unit; electric
# units sized so that each scheduled class fits its window,
# n = max ceil(m_prop / (mdot_unit window w)), at least one. Collision
# avoidance is not a unit driver; it is a feasibility verdict.
#
# Electronics are sized on rated DC power, the added power system on the
# peak DC power in use (the larger of the two modes).

import math
from dataclasses import dataclass
from typing import Dict, Optional


from ..constants import G0
from .allocation import transfer_ep_share
from .architectures import get_architecture
from .chemical import split_biprop
from .mass_models import (
    tank_mass,
    tank_technology,
    electronics_mass,
    eps_increment,
    transfer_array_derate,
    feed_mass,
    with_structure,
)
from .parameters import value


@dataclass
class ClosureResult:
    """Result of one closure; masses in kg.

    diverged : no fixed point (the map does not contract); detected by
               growth beyond 1e5 x bus mass or by a run of non-decreasing
               increments.
    ca_charge : for electric collision avoidance, {'factor': charged /
               impulsive delta-v, 'drift_ok': bool}; None otherwise.
               drift_ok False: the required drift is unreachable at any
               delta-v; the maximum-drift burn is loaded instead.
    arr_degr : array power loss applied to the added power system.
    """

    arch_key: str
    mission_key: str
    wet: float
    dry: float
    prop: float
    tanks: float
    eps_inc: float
    heads: float
    electronics: float
    feed: float
    p_dc_peak: float  # [W]
    n_units: Dict[str, int]  # {'ep': n, 'chem': n}
    per_class_prop: Dict[str, Dict[str, float]]  # class -> mode -> kg
    prop_by_species: Dict[str, float]  # 'n2o'/'propene'/'xe' -> kg
    tanks_by_species: Dict[str, float]
    n_iter: int
    converged: bool
    diverged: bool = False
    ca_charge: Optional[dict] = None
    arr_degr: float = 0.0


def sizing_window_s(cls, mission):
    """Thrusting time available to size the electric unit count on a class
    [s]; None for collision avoidance and its return."""
    if cls.kind in ("ca", "ca_restore"):
        return None
    if cls.kind == "sk":
        return value("sk_duty_max") * mission.duration_s * cls.w_duty
    return cls.window_s * cls.w_duty


def _ca_centroid_charge(dv_ev_imp, m_run, ve, mdot_total, t_warn, w_duty=1.0):
    """Per-event electric avoidance delta-v after the finite-burn charge.

    A finite tangential burn acts as the same impulse at its centroid. With
    thrusting in sunlight only, a burn of thrust-on time c dv_c takes
    c dv_c / w of wall-clock time, and the charged delta-v solves

        dv_c (t_warn - c dv_c / (2 w)) = dv_imp t_warn,

    c = m / (ve mdot). Smaller root, with c_eff = c / w:

        dv_c = [t_warn - sqrt(t_warn^2 - 2 c_eff dv_imp t_warn)] / c_eff.

    Without a real root the drift is unreachable before conjunction; the
    maximum-drift burn dv_c = w t_warn / c is returned.

    Returns (dv_ev_charged, drift_ok).
    """
    c_eff = m_run / (ve * mdot_total) / w_duty  # s per (m/s)
    disc = t_warn * t_warn - 2.0 * c_eff * dv_ev_imp * t_warn
    if disc >= 0.0:
        return (t_warn - math.sqrt(disc)) / c_eff, True
    return t_warn / c_eff, False


def _prop_masses(m0, arch, mission, allocation, n_ep_prev=1):
    """Propellant by class and mode, burned in execution order from m0.

    Returns (per_class {class: {mode: kg}}, totals {mode: kg}, ca_charge).
    Electric collision avoidance is charged the finite-burn factor with the
    unit count n_ep_prev of the previous iterate.
    """
    per_class = {c.key: {"chem": 0.0, "ep": 0.0} for c in mission.classes}
    totals = {"chem": 0.0, "ep": 0.0}
    ca_charge = None
    m_run = m0
    for cls in sorted(mission.classes, key=lambda c: c.execution_rank):
        for mode in ("chem", "ep"):
            dv = allocation[cls.key][mode]
            if dv <= 0.0:
                continue
            isp = arch.chem_unit.isp if mode == "chem" else arch.ep_unit.isp
            if cls.kind == "ca" and mode == "ep":
                # finite-burn charge per event
                dv_ev = dv / cls.n_events
                dv_ev_c, drift_ok = _ca_centroid_charge(
                    dv_ev,
                    m_run,
                    isp * G0,
                    max(n_ep_prev, 1) * arch.ep_unit.mdot,
                    mission.t_warn_s,
                    w_duty=cls.w_duty,
                )
                ca_charge = {"factor": dv_ev_c / dv_ev, "drift_ok": drift_ok}
                dv = dv_ev_c * cls.n_events
            mp = m_run * (1.0 - math.exp(-dv / (isp * G0)))
            per_class[cls.key][mode] = mp
            totals[mode] += mp
            m_run -= mp
    return per_class, totals, ca_charge


def _species_split(arch, totals):
    """Propellant by tank [kg] from the mode totals."""
    species = {"n2o": 0.0, "propene": 0.0, "xe": 0.0}
    if arch.has_chem:
        if arch.chem_kind == "biprop":
            ox, fuel = split_biprop(totals["chem"])
            species["n2o"] += ox
            species["propene"] += fuel
        else:  # 'mono': all N2O
            species["n2o"] += totals["chem"]
    if arch.has_ep:
        species[arch.ep_propellant] += totals["ep"]
    return species


def _ep_unit_count(arch, mission, per_class):
    """Electric unit count: the most any scheduled class needs, at least 1."""
    if not arch.has_ep:
        return 0
    n = 1
    for cls in mission.classes:
        mp = per_class[cls.key]["ep"]
        window = sizing_window_s(cls, mission)
        if mp <= 0.0 or window is None:
            continue
        n = max(n, math.ceil(mp / (arch.ep_unit.mdot * window)))
    return n


def close(
    arch_key,
    mission,
    allocation,
    *,
    arch=None,
    tank_tech_overrides=None,
    redundancy_plus_one=False,
    tol=1e-10,
    max_iter=4000,
):
    """Wet-mass closure of one design.

    arch : architecture to use instead of the registry entry (e.g. with
        the electric unit at another eta_u).
    tank_tech_overrides : tank type per propellant, e.g. {'n2o': 'copv'}.
    redundancy_plus_one : one extra unit of each installed mode.

    Returns a ClosureResult; wet is nan if the iteration neither converges
    nor diverges within max_iter.
    """
    arch = get_architecture(arch_key) if arch is None else arch
    bus = mission.bus_mass

    # Array power loss from the electric share of the transfer
    s_transfer = transfer_ep_share(mission, allocation)
    arr_degr = transfer_array_derate(mission, s_transfer)

    m_wet_prop = 10.0  # initial guess
    converged = False
    diverged = False
    n_ep = 1 if arch.has_ep else 0  # lagged into the CA charge
    # Divergence: a run of non-decreasing positive increments (a map with
    # slope >= 1), or growth beyond 1e5 x bus mass
    d_prev = None
    n_nondecr = 0
    _DIV_STREAK = 40
    for it in range(1, max_iter + 1):
        m0 = bus + m_wet_prop

        per_class, totals, ca_charge = _prop_masses(
            m0, arch, mission, allocation, n_ep_prev=n_ep
        )
        m_prop_tot = totals["chem"] + totals["ep"]
        species = _species_split(arch, totals)

        # Loaded propellant includes the unusable residual; tank volume
        # already carries it
        residual = value("residual_factor")
        m_loaded = m_prop_tot * residual
        species_loaded = {p: m * residual for p, m in species.items()}

        # Tanks of the installed set
        tanks_by = {}
        for prop in arch.tanks:
            tech = tank_technology(prop, tank_tech_overrides)
            tanks_by[prop] = tank_mass(species[prop], prop, tech, installed=True)
        tanks = sum(tanks_by.values())

        # Unit counts
        n_ep = _ep_unit_count(arch, mission, per_class)
        n_chem = arch.n_chem_units
        if redundancy_plus_one:
            n_ep = n_ep + 1 if arch.has_ep else 0
            n_chem = n_chem + 1 if arch.has_chem else 0

        heads = 0.0
        electronics = 0.0
        p_ep_peak = 0.0
        p_chem_peak = 0.0
        if arch.has_ep:
            heads += n_ep * arch.ep_unit.head_mass
            # electronics on rated power; added power system on used power
            electronics += electronics_mass(
                n_ep * arch.ep_unit.p_dc, arch.ep_unit.electronics
            )
            if totals["ep"] > 0.0:
                p_ep_peak = n_ep * arch.ep_unit.p_dc
        if arch.has_chem:
            heads += n_chem * arch.chem_unit.head_mass
            # catalyst heater of the monopropellant unit
            electronics += electronics_mass(
                n_chem * arch.chem_unit.p_dc, arch.chem_unit.electronics
            )
            if totals["chem"] > 0.0:
                p_chem_peak = n_chem * arch.chem_unit.p_dc
        # the two modes never fire together
        p_dc_peak = max(p_ep_peak, p_chem_peak)

        eps = eps_increment(p_dc_peak, mission.bus_class, arr_degr)
        feed = feed_mass(m_loaded)
        dry = with_structure(heads + electronics + eps + tanks + feed)

        new_wet = dry + m_loaded
        if new_wet > 1e5 * bus:  # slope > 1: no fixed point
            m_wet_prop = new_wet
            diverged = True
            break
        d = new_wet - m_wet_prop
        if abs(d) <= tol * max(1.0, new_wet):
            m_wet_prop = new_wet
            converged = True
            break
        if d_prev is not None and d > 0.0 and d >= d_prev:
            n_nondecr += 1
        else:
            n_nondecr = 0
        d_prev = d
        m_wet_prop = new_wet
        if n_nondecr >= _DIV_STREAK:
            diverged = True
            break

    # no wet mass without convergence or divergence
    wet_out = m_wet_prop if (converged or diverged) else float("nan")

    return ClosureResult(
        arch_key=arch_key,
        mission_key=mission.key,
        wet=wet_out,
        dry=dry,
        prop=m_loaded,
        tanks=tanks,
        eps_inc=eps,
        heads=heads,
        electronics=electronics,
        feed=feed,
        p_dc_peak=p_dc_peak,
        n_units={"ep": n_ep, "chem": n_chem},
        per_class_prop=per_class,
        prop_by_species=species_loaded,
        tanks_by_species=tanks_by,
        n_iter=it,
        converged=converged,
        diverged=diverged,
        ca_charge=ca_charge,
        arr_degr=arr_degr,
    )


def expected_idle_offsets():
    """Mass of the idle hardware of the bipropellant multimode [kg]:
    (electrodeless head and rated RF chain, at f_EP = 0; B20 head and
    empty propene tank, at f_EP = 1), with the structure margin."""
    arch_mm = get_architecture("multimode")
    ecr = arch_mm.ep_unit
    off0 = with_structure(ecr.head_mass + electronics_mass(ecr.p_dc, ecr.electronics))
    off1 = with_structure(
        arch_mm.chem_unit.head_mass
        + tank_mass(0.0, "propene", tank_technology("propene"), installed=True)
    )
    return off0, off1
