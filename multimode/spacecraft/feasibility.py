"""
Feasibility of a closed design: schedule checks per manoeuvre class, the
power-system limit and the bus mass-fraction limit, and the closed-form
transfer check of the GEO mission.
"""

#     station keeping   burn time <= sk_duty_max x duration x w_duty
#     collision avoid.  per-event wall-clock burn (t_on / w) <= 0.5 t_warn,
#                       and the required drift reachable
#     disposal, transfer  wall-clock burn (t_on / w) <= the class window
#     power             added power-system mass <= eps_cap_frac x bus mass
#     mass fraction     propulsion wet mass <= envelope_frac x bus mass
#
# A design failing any check is excluded, with the reason named.

import math
from dataclasses import dataclass
from typing import List

from ..constants import G0
from .architectures import get_architecture, ecr_unit, hall_unit
from .mass_models import eps_increment, eps_mass_cap, electronics_mass
from .parameters import value, YEAR_S, dv_transfer_ep


@dataclass(frozen=True)
class Verdict:
    """One (class, mode) feasibility verdict."""

    class_key: str
    mode: str
    t_burn_s: float  # [s]
    t_limit_s: float  # [s]
    feasible: bool
    reason: str


@dataclass
class FeasibilityReport:
    arch_key: str
    mission_key: str
    verdicts: List[Verdict]
    eps_inc: float
    eps_cap: float
    power_ok: bool
    # bus mass-fraction limit
    wet: float = float("nan")
    envelope_cap: float = float("inf")
    envelope_ok: bool = True

    @property
    def excluded(self):
        return (
            (not self.power_ok)
            or (not self.envelope_ok)
            or any(not v.feasible for v in self.verdicts)
        )

    @property
    def reasons(self):
        r = [
            f"{v.class_key}/{v.mode}: {v.reason}"
            for v in self.verdicts
            if not v.feasible
        ]
        if not self.power_ok:
            r.append(
                f"power: EPS increment {self.eps_inc:.1f} kg exceeds "
                f"the power-system limit {self.eps_cap:.1f} kg"
            )
        if not self.envelope_ok:
            r.append(
                f"envelope: propulsion wet mass {self.wet:.1f} kg "
                f"exceeds envelope_frac x bus = "
                f"{self.envelope_cap:.1f} kg"
            )
        return r


def _mode_unit_and_count(arch, mode, result):
    if mode == "ep":
        return arch.ep_unit, result.n_units["ep"]
    return arch.chem_unit, result.n_units["chem"]


def assess(
    arch_key, mission, allocation, result, t_warn=None, arch=None, envelope_frac=None
):
    """Feasibility of a converged closure result.

    t_warn, arch and envelope_frac replace the mission warning time, the
    registry architecture and the register mass-fraction limit.
    """
    arch = get_architecture(arch_key) if arch is None else arch
    if t_warn is None:
        t_warn = mission.t_warn_s if mission.t_warn_s is not None else value("t_warn")

    verdicts = []
    for cls in mission.classes:
        for mode in ("chem", "ep"):
            mp = result.per_class_prop[cls.key][mode]
            if mp <= 0.0:
                continue
            if cls.kind == "ca_restore":
                continue
            unit, n = _mode_unit_and_count(arch, mode, result)
            mdot = n * unit.mdot

            if cls.kind == "ca":
                # per-event wall-clock burn of the charged delta-v
                t_burn = (mp / cls.n_events) / mdot / cls.w_duty
                t_limit = value("ca_burn_twarn_frac") * t_warn
                drift_ok = result.ca_charge is None or result.ca_charge["drift_ok"]
                ok = drift_ok and t_burn <= t_limit
                if ok:
                    reason = "ok"
                elif not drift_ok:
                    reason = (
                        "required drift unreachable before "
                        "conjunction at any Delta-v "
                        f"(burn ~ t_warn = {t_warn / 3600:.1f} h)"
                    )
                else:
                    reason = (
                        f"charged per-event wall-clock burn "
                        f"{t_burn / 3600:.2f} h exceeds "
                        f"0.5 x t_warn = {t_limit / 3600:.2f} h"
                    )
                verdicts.append(Verdict(cls.key, mode, t_burn, t_limit, ok, reason))
            elif cls.kind == "sk":
                t_burn = mp / mdot
                t_limit = value("sk_duty_max") * mission.duration_s * cls.w_duty
                ok = t_burn <= t_limit
                duty = t_burn / (mission.duration_s * cls.w_duty)
                verdicts.append(
                    Verdict(
                        cls.key,
                        mode,
                        t_burn,
                        t_limit,
                        ok,
                        "ok"
                        if ok
                        else f"SK duty {duty:.2f} exceeds {value('sk_duty_max')}",
                    )
                )
            else:  # transfer, disposal
                t_wall = (mp / mdot) / cls.w_duty
                t_limit = cls.window_s
                ok = t_wall <= t_limit
                verdicts.append(
                    Verdict(
                        cls.key,
                        mode,
                        t_wall,
                        t_limit,
                        ok,
                        "ok"
                        if ok
                        else f"wall-clock {t_wall / YEAR_S:.2f} yr exceeds the "
                        f"{cls.kind} window {t_limit / YEAR_S:.2f} yr",
                    )
                )

    cap = eps_mass_cap(mission.bus_mass)
    if envelope_frac is None:
        envelope_frac = value("envelope_frac")
    env_cap = envelope_frac * mission.bus_mass
    return FeasibilityReport(
        arch_key=arch_key,
        mission_key=mission.key,
        verdicts=verdicts,
        eps_inc=result.eps_inc,
        eps_cap=cap,
        power_ok=result.eps_inc <= cap,
        wet=result.wet,
        envelope_cap=env_cap,
        envelope_ok=result.wet <= env_cap,
    )


# Closed-form transfer check, GEO mission


def transfer_asymmetry_check(
    m_avg=None, t_max_yr=None, w_duty=None, bus_class=None, bus_mass=None, arr_degr=0.0
):
    """Hall and electrodeless fleets that fly the electric GTO to GEO
    transfer within the window, sized on thrust at an average mass:
    T_req = m_avg dv_ep / (t_max w), and the chemical transfer propellant.

    arr_degr : array power loss applied to the added power system.

    Returns {'T_req_N', 'hall', 'ecr', 'chem', 'eps_cap'}; each fleet with
    its unit count, DC power, added power-system mass, electronics mass
    and power verdict.
    """
    m_avg = value("geo_m_avg") if m_avg is None else m_avg
    t_max = (value("t_transfer_max") if t_max_yr is None else t_max_yr) * YEAR_S
    w = value("w_duty_transfer") if w_duty is None else w_duty
    bus_class = int(value("bus_mass_geo")) if bus_class is None else bus_class
    bus_mass = value("bus_mass_geo") if bus_mass is None else bus_mass

    a_req = dv_transfer_ep() / (t_max * w)
    T_req = a_req * m_avg
    cap = eps_mass_cap(bus_mass)

    out = {"T_req_N": T_req}
    for name, unit in (("hall", hall_unit()), ("ecr", ecr_unit())):
        n = math.ceil(T_req / unit.thrust)
        p_dc = n * unit.p_dc
        eps = eps_increment(p_dc, bus_class, arr_degr)
        out[name] = {
            "n_units": n,
            "p_dc_w": p_dc,
            "eps_inc": eps,
            "electronics": electronics_mass(p_dc, unit.electronics),
            "power_ok": eps <= cap,
        }
    m_chem = value("geo_chem_stack") * (
        1.0 - math.exp(-value("dv_transfer_chem") / (value("isp_chem") * G0))
    )
    out["chem"] = {"m_prop": m_chem, "stack": value("geo_chem_stack")}
    out["eps_cap"] = cap
    return out
