"""
Spacecraft-model checks: reference designs of both missions, closed-form and
analytic cross-checks, sweep structure and the sensitivity mechanics.
Tolerances are given with each target.
"""

import math

from ..spacecraft.allocation import nominal_allocation, sweep_allocation
from ..spacecraft.architectures import ARCHITECTURE_ORDER, get_architecture
from ..spacecraft.closure import close, expected_idle_offsets
from ..spacecraft.feasibility import assess, transfer_asymmetry_check
from ..spacecraft.mass_models import with_structure
from ..spacecraft.mission import leo_mission, geo_mission
from ..spacecraft.parameters import value, dv_transfer_ep
from ..spacecraft.sweeps import (
    sweep_f_ep,
    sweep_eta_u,
    sweep_t_warn,
    sweep_t_transfer_max,
    tornado,
    _close_arch_at_eta_u,
    mono_corner_grid,
    T_WARN_SET_S,
)
from .common import section
from .registry import task
from .spacecraft_results import (
    nominal_designs_leo,
    nominal_designs_geo,
    multimode_all_ep_transfer,
)


# Reference values (regression targets)

# LEO mission, nominal allocation, metal tanks:
# (wet, dry, loaded propellant, tanks, added power system) [kg]
_LEO_TARGETS = {
    "single_elec": (7.04, 5.72, 1.31, 0.80, 1.44),
    "multimode": (6.54, 5.46, 1.09, 1.60, 0.00),
    "hybrid": (9.90, 8.49, 1.41, 2.40, 1.44),
    "single_chem": (8.67, 3.67, 5.00, 1.99, 0.00),
    "mono_multimode": (5.41, 4.27, 1.14, 0.80, 0.00),
}
# Wet-mass differences, LEO [kg]
_LEO_DELTA_TARGETS = {"mm_minus_se": -0.50, "hy_minus_mm": +3.36, "mmm_minus_mm": -1.13}
# LEO wet mass with a composite N2O tank [kg]
_LEO_COPV_TARGETS = {"multimode": 6.54, "hybrid": 9.90, "single_chem": 8.19}
_ATOL = 0.02  # kg

# Closed-form transfer check, GEO: thrust [mN], and per fleet the unit
# count, DC power [kW], added power system and electronics [kg]
_TRANSFER_CHECK_TARGETS = {
    "T_req_mN": 77.976,
    "hall": {"n": 6, "p_kw": 1.2903, "eps": 13.63, "elec": 8.00},
    "ecr": {"n": 80, "p_kw": 7.5446, "eps": 91.81, "elec": 93.55},
    "chem_prop": 211.6,
}
# Idle hardware of the multimode at f_EP = 0 and 1 [kg]; tolerance 20 %
_IDLE_OFFSET_F0 = 2.603
_IDLE_OFFSET_F1 = 1.794
_IDLE_OFFSET_TOL = 0.20

# GEO mission, nominal allocation: ((wet, dry, propellant, tanks, added
# power system) [kg], electric unit count)
_GEO_TARGETS = {
    "single_elec": ((185.80, 68.17, 117.63, 5.34, 23.18), 8),
    "multimode": ((747.89, 183.05, 564.84, 105.15, 0.00), 2),
    "hybrid": ((753.56, 178.76, 574.80, 102.54, 0.19), 1),
    "single_chem": ((1731.13, 412.06, 1319.07, 243.76, 0.00), 0),
    "mono_multimode": ((13220.02, 3342.29, 9877.73, 1999.54, 14.00), 14),
}
# GEO multimode with an all-electric transfer: 'diverged' or
# 'converged_power_excluded'
_GEO_ALL_EP_TRANSFER_STATE = "diverged"

# Single-electric collision-avoidance verdicts at 2 h, 6 h, 1 d, 3 d and
# the finite-burn factor at 1 d
_CA_VERDICT_TARGETS = (False, True, True, True)
_CA_FACTOR_1D = 1.0113

# Electrodeless unit against eta_u: eta_u -> (Isp [s], mdot [kg/s])
_ECR_ETA_U_TARGETS = {
    0.20: (592.7868852459016, 1.686650e-07),
    0.85: (2519.3442622950818, 3.9685882352941176e-08),
}
# Thruster model at parity against the register rows (relative)
_ECR_PARITY_RTOL = 2e-3
# eta_u at which the electric Isp equals the bipropellant Isp
_ECR_ISP_CHEM_ETA_U = 0.0914242

# eta_u thresholds of the multimode (tolerance 1e-3)
_MIN_VIABLE_TARGETS = {
    "leo": {"own_chem": 0.12656, "single_chem": 0.21481},
    "geo": {"own_chem": 0.09534, "single_chem": 0.09628},
}
# Multimode wet mass at single eta_u [kg]
_ETA_U_SPOT_TARGETS = {
    "leo": {0.05: 22.05, 0.61: 6.54, 0.85: 6.23},
    "geo": {0.07: 2862.07, 0.61: 747.89, 0.85: 718.98},
}
# Multimode at f_EP = 0 [kg]
_ETA_U_F0_REFS = {"leo": 11.38, "geo": 1748.67}

# Monopropellant multimode: LEO eta_u thresholds, LEO all-chemical
# reference [kg], and the eta_u points at which the GEO closure diverges
_MONO_MIN_VIABLE_LEO = {"own_chem": 0.07299, "single_chem": 0.16584}
_MONO_F0_REF_LEO = 15.1955
_MONO_GEO_ETA_U_DIVERGED = [
    0.05,
    0.07,
    0.09,
    0.11,
    0.13,
    0.15,
    0.17,
    0.19,
    0.21,
    0.23,
    0.25,
]

# GEO designs excluded by the mass-fraction limit at each limit value
_GEO_ENVELOPE_TARGETS = {
    0.35: {"single_elec": True, "multimode": True, "hybrid": True, "single_chem": True},
    0.50: {"single_elec": True, "multimode": True, "hybrid": True, "single_chem": True},
    1.00: {
        "single_elec": False,
        "multimode": True,
        "hybrid": True,
        "single_chem": True,
    },
}
# GEO wet mass with a composite N2O tank, and multimode - hybrid [kg]
_GEO_COPV_TARGETS = {
    "single_elec": 185.80,
    "multimode": 498.61,
    "hybrid": 510.01,
    "single_chem": 869.23,
}
_GEO_COPV_GAP = -11.40

# f_EP points (21-point grid) at which the GEO closure diverges
_F_EP_DIVERGENCE_TARGETS = {
    ("mono_multimode", "geo"): (
        0.0,
        0.05,
        0.1,
        0.15,
        0.2,
        0.25,
        0.8,
        0.85,
        0.9,
        0.95,
        1.0,
    ),
    ("multimode", "geo"): (0.9, 0.95, 1.0),
}

# Warning-time sweep: verdicts and lightest feasible design at 2 h, 6 h,
# 1 d, 3 d; single-electric limit [s] and tolerance; multimode at 6 h [kg]
_T_WARN_VERDICT_TARGETS = {
    "single_elec": (False, True, True, True),
    "multimode": (True, True, True, True),
    "hybrid": (True, True, True, True),
    "single_chem": (True, True, True, True),
    "mono_multimode": (True, True, True, True),
}
_T_WARN_OPTIMAL_TARGETS = (
    "mono_multimode",
    "mono_multimode",
    "mono_multimode",
    "mono_multimode",
)
_T_WARN_FLIP_S = 5.838 * 3600.0
_T_WARN_FLIP_ATOL_S = 0.5 * 3600.0
_T_WARN_MM_6H_WET = 6.845

# Transfer-window sweep. Closed-form check at 3 and 24 months: (units,
# DC power [kW], added power system [kg], electronics [kg], power verdict)
_WINDOW_CHECK_TARGETS = {
    0.25: {
        "hall": (12, 2.5806, 29.76, 16.00, True),
        "ecr": (160, 15.0892, 186.12, 187.11, False),
    },
    2.00: {
        "hall": (2, 0.4301, 2.88, 2.67, True),
        "ecr": (20, 1.8862, 21.08, 23.39, True),
    },
}
# Single-electric Hall units over the windows (None: no fixed point)
_WINDOW_SE_N_HALL = (None, None, 318, 18, 8, 4, 2)
# Time-optimal multimode and hybrid wet masses [kg]; monopropellant
# multimode wet mass [kg] and units
_WINDOW_MM_HY_WET = (747.89, 753.56)
_WINDOW_MMM_NOMINAL = (13220.02, 14)
_WINDOW_VERDICT_TARGETS = {
    "single_elec": (False, False, False, False, False, True, True),
    "multimode": (False,) * 7,
    "hybrid": (False,) * 7,
    "mono_multimode": (False,) * 7,
}
# Maximum-electric transfer: multimode at 24 months (share, wet [kg]), hybrid
# wet masses [kg] by month, monopropellant multimode at 24 months (share,
# units) and hybrid verdicts
_WINDOW_EPTR_MM_24MO_S = 0.78967
_WINDOW_EPTR_MM_24MO_WET = 329.03
_WINDOW_EPTR_HY_WETS = {3.0: 406.85, 6.0: 189.61, 12.0: 151.72, 24.0: 132.78}
_WINDOW_EPTR_MMM_24MO = (0.72626, 24)
_WINDOW_EPTR_HY_VERDICTS = (False, False, False, False, False, False, True)
_WINDOW_EPTR_MM_ALL_EXCLUDED = True

# Hybrid with a one-month window, maximum-electric transfer: electric share,
# wet mass [kg], unit ceilings (Hall, electrodeless) without and with the
# nominal array loss, Hall units, added power system [kg], electric
# transfer propellant [kg]; and the 24-month fully electric wet mass [kg]
_SPLIT_S = 0.182129
_SPLIT_WET = 741.38
_SPLIT_N_MAX = (12, 27)
_SPLIT_N_MAX_NOMINAL = (10, 22)
_SPLIT_N_HALL = 11
_SPLIT_EPS = 27.99
_SPLIT_TR_EP = 22.07
_FULL_EP_24MO_WET = 132.7756

# Sensitivity method [kg]
_TORNADO_D_MMM_MM_GEO = 12472.14
_ETA_U_GEO_005_WET = 10514.58
_TORNADO_PF6350_SC_WET = 8.447
_TORNADO_MARGIN50_SE_DWET = 0.7646
_TORNADO_HEAD_ECR_MM_DWET = 0.6140

# Monopropellant corner grid: (head level, Isp level) -> wet [kg];
# unfavourable-corner difference and interaction term [kg]
_CORNER_GRID_TARGETS = {
    ("lo", "lo"): 5.1441,
    ("lo", "nom"): 5.1049,
    ("lo", "hi"): 5.0844,
    ("nom", "lo"): 5.4466,
    ("nom", "nom"): 5.4074,
    ("nom", "hi"): 5.3868,
    ("hi", "lo"): 5.7491,
    ("hi", "nom"): 5.7098,
    ("hi", "hi"): 5.6892,
}
_CORNER_D_UNFAV = -0.7916
_CORNER_INTERACTION = +0.000076
_CORNER_INTERACTION_ATOL = 2e-4


# Checks

_FAILURES = []


def _check(label, ok, detail=""):
    print(f"  [{'OK' if ok else 'FAIL'}] {label}" + (f" — {detail}" if detail else ""))
    if not ok:
        _FAILURES.append(label)


@task("check_spacecraft", "Spacecraft-model checks", in_all=False)
def run_spacecraft_checks(ctx):
    _FAILURES.clear()

    section("LEO mission results")
    mission, runs = nominal_designs_leo()
    for key, tgt in _LEO_TARGETS.items():
        r = runs[key][1]
        got = (r.wet, r.dry, r.prop, r.tanks, r.eps_inc)
        err = max(abs(g - t) for g, t in zip(got, tgt))
        _check(f"{key} wet/dry/prop/tanks/eps", err <= _ATOL, f"max err {err:.3f} kg")
        _check(f"{key} converged", r.converged, f"{r.n_iter} iterations")

    section("Differences between architectures")
    d_mm_se = runs["multimode"][1].wet - runs["single_elec"][1].wet
    d_hy_mm = runs["hybrid"][1].wet - runs["multimode"][1].wet
    _check(
        "mm - se delta",
        abs(d_mm_se - _LEO_DELTA_TARGETS["mm_minus_se"]) <= _ATOL,
        f"{d_mm_se:+.3f} kg (target {_LEO_DELTA_TARGETS['mm_minus_se']:+.2f})",
    )
    _check(
        "hy - mm delta",
        abs(d_hy_mm - _LEO_DELTA_TARGETS["hy_minus_mm"]) <= _ATOL,
        f"{d_hy_mm:+.3f} kg (target {_LEO_DELTA_TARGETS['hy_minus_mm']:+.2f})",
    )
    d_mmm_mm = runs["mono_multimode"][1].wet - runs["multimode"][1].wet
    _check(
        "mono_mm - mm delta",
        abs(d_mmm_mm - _LEO_DELTA_TARGETS["mmm_minus_mm"]) <= _ATOL,
        f"{d_mmm_mm:+.3f} kg (target {_LEO_DELTA_TARGETS['mmm_minus_mm']:+.2f})",
    )
    for key in ARCHITECTURE_ORDER:
        n = runs[key][1].n_units["ep"]
        expected = 1 if get_architecture(key).has_ep else 0
        _check(f"{key} n_ep == {expected}", n == expected, f"got {n}")

    section("Tank sizing, LEO")
    _, runs_copv = nominal_designs_leo(tank_tech_overrides={"n2o": "copv"})
    for key, tgt in _LEO_COPV_TARGETS.items():
        w = runs_copv[key][1].wet
        _check(f"COPV {key} wet", abs(w - tgt) <= _ATOL, f"{w:.3f} kg (target {tgt})")
    _check(
        "COPV ordering mm < hy preserved",
        runs_copv["multimode"][1].wet < runs_copv["hybrid"][1].wet,
    )

    section("Power of the two electric units")
    q = transfer_asymmetry_check()
    t = _TRANSFER_CHECK_TARGETS
    _check(
        "T_req",
        abs(q["T_req_N"] * 1e3 - t["T_req_mN"]) < 0.1,
        f"{q['T_req_N'] * 1e3:.3f} mN",
    )
    for name in ("hall", "ecr"):
        d, tt = q[name], t[name]
        _check(
            f"{name} n_units == {tt['n']}",
            d["n_units"] == tt["n"],
            f"got {d['n_units']}",
        )
        _check(
            f"{name} P_dc",
            abs(d["p_dc_w"] / 1e3 - tt["p_kw"]) < 0.01,
            f"{d['p_dc_w'] / 1e3:.4f} kW",
        )
        _check(
            f"{name} EPS increment",
            abs(d["eps_inc"] - tt["eps"]) < 0.1,
            f"{d['eps_inc']:.2f} kg",
        )
        _check(
            f"{name} electronics",
            abs(d["electronics"] - tt["elec"]) < 0.1,
            f"{d['electronics']:.2f} kg",
        )
    _check(
        "power-system limit separates the fleets",
        q["hall"]["power_ok"] and not q["ecr"]["power_ok"],
        f"hall ok={q['hall']['power_ok']}, ecr ok={q['ecr']['power_ok']}",
    )
    _check(
        "chem transfer propellant",
        abs(q["chem"]["m_prop"] - t["chem_prop"]) < 0.5,
        f"{q['chem']['m_prop']:.1f} kg",
    )

    section("Allocation limits")
    off0, off1 = expected_idle_offsets()
    r_mm0 = close("multimode", mission, sweep_allocation("multimode", mission, 0.0))
    r_sc = runs["single_chem"][1]
    gap0 = r_mm0.wet - r_sc.wet
    _check(
        "f_EP = 0 offset within 20% of expected",
        abs(gap0 - off0) <= _IDLE_OFFSET_TOL * off0,
        f"gap {gap0:.3f} kg vs expected {off0:.3f} kg",
    )
    _check(
        "expected f0 offset matches reference value",
        abs(off0 - _IDLE_OFFSET_F0) < 0.01,
        f"{off0:.3f} kg",
    )
    _check(
        "expected f1 offset matches reference value",
        abs(off1 - _IDLE_OFFSET_F1) < 0.01,
        f"{off1:.3f} kg",
    )
    r_mm1 = close("multimode", mission, sweep_allocation("multimode", mission, 1.0))
    r_se = runs["single_elec"][1]
    gap1 = r_mm1.wet - r_se.wet
    _floor = value("m_min_tank") + value("m_fixed_tank")
    _b20 = get_architecture("multimode").chem_unit.head_mass
    named = (
        off1
        + (r_mm1.prop - r_se.prop)
        + with_structure(
            r_mm1.tanks
            - r_se.tanks
            - _floor
            + r_mm1.electronics
            - r_se.electronics
            + r_mm1.eps_inc
            - r_se.eps_inc
            + r_mm1.heads
            - r_se.heads
            - _b20
            + r_mm1.feed
            - r_se.feed
        )
    )
    resid = gap1 - named
    _check(
        "f_EP = 1 gap decomposes into named terms (residual <= 20% of hardware offset)",
        abs(resid) <= _IDLE_OFFSET_TOL * off1,
        f"gap {gap1:.3f} = named {named:.3f} + residual {resid:.3f} kg",
    )

    section("Allocation sweep")
    ok_mono, ok_conv = True, True
    prev_ep, prev_ch = -1.0, float("inf")
    for i in range(11):
        f = i / 10.0
        al = sweep_allocation("multimode", mission, f)
        ep = sum(v["ep"] for v in al.values())
        ch = sum(v["chem"] for v in al.values())
        ok_mono &= (ep >= prev_ep - 1e-9) and (ch <= prev_ch + 1e-9)
        prev_ep, prev_ch = ep, ch
        ok_conv &= close("multimode", mission, al).converged
    _check(
        "EP loss-inclusive dv non-decreasing, chem non-increasing over f_EP grid",
        ok_mono,
    )
    _check("all f_EP grid closures converge", ok_conv)

    section("GEO transfer mission results")
    mission_b, runs_b = nominal_designs_geo()
    for key, (tgt, n_tgt) in _GEO_TARGETS.items():
        r = runs_b[key][1]
        got = (r.wet, r.dry, r.prop, r.tanks, r.eps_inc)
        err = max(abs(g - t) for g, t in zip(got, tgt))
        _check(
            f"GEO {key} wet/dry/prop/tanks/eps", err <= _ATOL, f"max err {err:.3f} kg"
        )
        _check(f"GEO {key} converged", r.converged, f"{r.n_iter} iterations")
        _check(
            f"GEO {key} n_ep == {n_tgt}",
            r.n_units["ep"] == n_tgt,
            f"got {r.n_units['ep']}",
        )
    r_mp = runs_b["mono_multimode"][1]
    _check(
        "GEO mono_multimode closes and is envelope-excluded",
        r_mp.converged
        and not r_mp.diverged
        and runs_b["mono_multimode"][2].excluded
        and not runs_b["mono_multimode"][2].envelope_ok,
        f"wet {r_mp.wet:.2f} kg, n_ECR {r_mp.n_units['ep']}, "
        f"{r_mp.wet / runs_b['mono_multimode'][2].envelope_cap:.1f}x "
        f"the cap",
    )
    for key, feasible in (
        ("single_elec", False),
        ("multimode", False),
        ("hybrid", False),
        ("single_chem", False),
    ):
        rep = runs_b[key][2]
        if feasible:
            _check(f"GEO {key} nominal feasible", not rep.excluded)
        else:
            _check(
                f"GEO {key} nominal excluded by the envelope only",
                rep.excluded
                and (not rep.envelope_ok)
                and rep.power_ok
                and all(v.feasible for v in rep.verdicts),
                f"wet {rep.wet:.1f} kg vs cap {rep.envelope_cap:.0f} "
                f"kg (mass still reported)",
            )
    r_cf, rep_cf = multimode_all_ep_transfer(mission_b)
    if _GEO_ALL_EP_TRANSFER_STATE == "diverged":
        _check(
            "GEO multimode all-EP transfer has no fixed point (diverged and flagged)",
            r_cf.diverged and not r_cf.converged,
            f"converged={r_cf.converged}, diverged={r_cf.diverged}, "
            f"{r_cf.n_iter} iterations",
        )
    else:
        _check(
            "GEO multimode all-EP transfer converges and "
            "is excluded by the power-system limit",
            r_cf.converged
            and not r_cf.diverged
            and (not rep_cf.power_ok)
            and r_cf.eps_inc > 2.0 * rep_cf.eps_cap,
            f"converged={r_cf.converged}, diverged={r_cf.diverged}; "
            f"EPS +{r_cf.eps_inc:.1f} kg vs {rep_cf.eps_cap:.0f} kg cap",
        )

    section("Collision-avoidance check")
    verdicts, factor_1d = [], None
    for tw in T_WARN_SET_S:
        m_tw = leo_mission(t_warn_s=tw)
        al = nominal_allocation("single_elec", m_tw)
        r = close("single_elec", m_tw, al)
        rep = assess("single_elec", m_tw, al, r)
        v = next(v for v in rep.verdicts if v.class_key == "ca")
        verdicts.append(v.feasible)
        if tw == 86400.0:
            factor_1d = r.ca_charge["factor"]
        al_mm = nominal_allocation("multimode", m_tw)
        r_mm = close("multimode", m_tw, al_mm)
        rep_mm = assess("multimode", m_tw, al_mm, r_mm)
        v_mm = next(v for v in rep_mm.verdicts if v.class_key == "ca")
        _check(
            f"chem CA feasible at t_warn = {tw / 3600:.0f} h",
            v_mm.feasible,
            f"burn {v_mm.t_burn_s:.0f} s",
        )
    _check(
        f"single_elec CA verdicts over {{2 h, 6 h, 1 d, 3 d}} == {_CA_VERDICT_TARGETS}",
        tuple(verdicts) == _CA_VERDICT_TARGETS,
        f"got {tuple(verdicts)}",
    )
    _check(
        "1 d centroid-charge factor",
        abs(factor_1d - _CA_FACTOR_1D) < 0.005,
        f"{factor_1d:.4f} (target {_CA_FACTOR_1D})",
    )

    section("Allocation sweep structure")
    ok_flags = True
    div_found = {}
    for m_i in (mission, mission_b):
        sw = sweep_f_ep(ARCHITECTURE_ORDER, m_i, n_points=21)
        for key, data in sw.items():
            for rec in data["records"]:
                if rec["diverged"]:
                    div_found.setdefault((key, m_i.key), []).append(
                        round(rec["f_ep"], 2)
                    )
                if not rec["diverged"] and not rec["converged"]:
                    ok_flags = False
                if not rec["diverged"] and math.isnan(rec["wet"]):
                    ok_flags = False
    _check("all sweep closures converge (both missions)", ok_flags)
    _check(
        "the divergences on the grid are exactly the reference set "
        "(mono_multimode and multimode, GEO only)",
        {k: tuple(v) for k, v in div_found.items()} == _F_EP_DIVERGENCE_TARGETS,
        f"diverged at {div_found}",
    )
    sw_mmm_b = sweep_f_ep(("mono_multimode",), mission_b, n_points=21)
    div_set = tuple(
        round(r["f_ep"], 2)
        for r in sw_mmm_b["mono_multimode"]["records"]
        if r["diverged"]
    )
    _check(
        "GEO mono_multimode f_EP divergence set == {0.00 .. 0.25} and {0.80 .. 1.00}",
        div_set == _F_EP_DIVERGENCE_TARGETS[("mono_multimode", "geo")],
        f"diverged at f = {div_set}",
    )
    ok_mono = True
    prev_ep, prev_ch = -1.0, float("inf")
    for i in range(21):
        al = sweep_allocation("multimode", mission_b, i / 20.0)
        ep = sum(v["ep"] for v in al.values())
        ch = sum(v["chem"] for v in al.values())
        ok_mono &= (ep >= prev_ep - 1e-9) and (ch <= prev_ch + 1e-9)
        prev_ep, prev_ch = ep, ch
    _check("GEO allocation monotone over the f_EP grid", ok_mono)
    off0, off1 = expected_idle_offsets()
    ecr = get_architecture("multimode").ep_unit
    b20_head = get_architecture("multimode").chem_unit.head_mass
    rf_rated = ecr.p_dc * value("sm_rf") / 1e3
    r_mm0 = close("multimode", mission_b, sweep_allocation("multimode", mission_b, 0.0))
    r_sc = runs_b["single_chem"][1]
    gap0 = r_mm0.wet - r_sc.wet
    named0 = (
        off0
        + (r_mm0.prop - r_sc.prop)
        + with_structure(
            r_mm0.tanks
            - r_sc.tanks
            + r_mm0.eps_inc
            - r_sc.eps_inc
            + (r_mm0.heads - r_sc.heads - ecr.head_mass)
            + (r_mm0.electronics - r_sc.electronics - rf_rated)
            + r_mm0.feed
            - r_sc.feed
        )
    )
    _check(
        "GEO f_EP = 0 gap decomposes (residual <= 20% of hardware offset)",
        abs(gap0 - named0) <= _IDLE_OFFSET_TOL * off0,
        f"gap {gap0:.3f} = named {named0:.3f} + residual {gap0 - named0:.3f} kg",
    )
    r_mm1 = close("multimode", mission_b, sweep_allocation("multimode", mission_b, 1.0))
    r_se = runs_b["single_elec"][1]
    gap1 = r_mm1.wet - r_se.wet
    _floor_b = value("m_min_tank") + value("m_fixed_tank")
    named1 = (
        off1
        + (r_mm1.prop - r_se.prop)
        + with_structure(
            r_mm1.tanks
            - r_se.tanks
            - _floor_b
            + r_mm1.electronics
            - r_se.electronics
            + r_mm1.eps_inc
            - r_se.eps_inc
            + r_mm1.heads
            - r_se.heads
            - b20_head
            + r_mm1.feed
            - r_se.feed
        )
    )
    _check(
        "GEO f_EP = 1 gap decomposes (residual <= 20% of hardware offset)",
        abs(gap1 - named1) <= _IDLE_OFFSET_TOL * off1,
        f"gap {gap1:.3f} = named {named1:.3f} + residual {gap1 - named1:.3f} kg",
    )

    section("Interface to the thruster-level model")
    from ..spacecraft.architectures import ecr_unit
    from ..thruster.ecr_unit import ecr_unit_at, eta_u_parity

    u0 = ecr_unit()
    u_par = ecr_unit(eta_u_parity())
    _check(
        "ecr_unit(0.61) == ecr_unit() on every field (exact)",
        all(
            getattr(u0, a) == getattr(u_par, a)
            for a in ("isp", "thrust", "mdot", "p_dc", "head_mass")
        ),
        "register-anchored ratios are exactly 1.0 at parity",
    )
    for eu, (isp_t, mdot_t) in _ECR_ETA_U_TARGETS.items():
        u = ecr_unit(eu)
        _check(
            f"interface Isp at eta_u = {eu}",
            abs(u.isp / isp_t - 1.0) < 1e-9,
            f"{u.isp:.4f} s",
        )
        _check(
            f"interface mdot at eta_u = {eu}",
            abs(u.mdot / mdot_t - 1.0) < 1e-9,
            f"{u.mdot * 1e9:.4f} ug/s",
        )
        _check(
            f"T and P_dc held at eta_u = {eu} (W/mN invariant)",
            u.thrust == u0.thrust and u.p_dc == u0.p_dc,
            "thrust and power held, the eta_u penalty is in Isp",
        )
    op = ecr_unit_at(eta_u_parity())
    rows = (
        (op.thrust, value("t_unit_ecr"), "thrust"),
        (op.isp, value("isp_ecr"), "Isp"),
        (op.mdot, value("mdot_unit_ecr"), "mdot"),
        (op.p_abs, value("p_abs_ecr"), "P_abs"),
    )
    _check(
        "raw pipeline parity op_point within "
        f"{_ECR_PARITY_RTOL:.0e} of the register rows",
        all(abs(a / b - 1.0) < _ECR_PARITY_RTOL for a, b, _ in rows),
        f"max rel dev {max(abs(a / b - 1.0) for a, b, _ in rows):.2e}",
    )
    eta_star = eta_u_parity() * value("isp_chem") / value("isp_ecr")
    _check(
        f"system EP Isp == chemical Isp at eta_u = {_ECR_ISP_CHEM_ETA_U}",
        abs(eta_star - _ECR_ISP_CHEM_ETA_U) < 1e-5
        and abs(ecr_unit(eta_star).isp - value("isp_chem")) < 1e-6,
        f"eta_u* = {eta_star:.6f}",
    )

    section("Utilisation sweep")
    for m_i in (mission, mission_b):
        sw = sweep_eta_u(m_i)
        recs = {r["eta_u"]: r for r in sw["records"]}
        n_nc = [r["eta_u"] for r in sw["records"] if not r["converged"]]
        if m_i.key == "leo":
            _check(
                "LEO eta_u grid: 41/41 converge, none excluded",
                not n_nc and not any(r["excluded"] for r in sw["records"]),
            )
        else:
            _check(
                "GEO eta_u grid: 41/41 converge, and eta_u = 0.05 "
                f"closes at {_ETA_U_GEO_005_WET} kg",
                not n_nc and abs(recs[0.05]["wet"] - _ETA_U_GEO_005_WET) <= 0.5,
                f"wet at eta_u = 0.05 is {recs[0.05]['wet']:.2f} kg",
            )
            _check(
                "GEO eta_u grid: every point excluded (the mass-"
                "fraction limit binds across the whole range)",
                all(r["excluded"] for r in sw["records"]),
            )
        for eu, tgt in _ETA_U_SPOT_TARGETS[m_i.key].items():
            w = recs[eu]["wet"]
            tol = 0.02 if tgt < 100 else 0.5
            _check(
                f"{m_i.key} sweep wet at eta_u = {eu}",
                abs(w - tgt) <= tol,
                f"{w:.3f} kg (target {tgt})",
            )
        _check(
            f"{m_i.key} f_EP = 0 reference",
            abs(sw["references"]["multimode_f0"] - _ETA_U_F0_REFS[m_i.key]) <= 0.02,
            f"{sw['references']['multimode_f0']:.3f} kg",
        )
        for crit, tgt in _MIN_VIABLE_TARGETS[m_i.key].items():
            x = sw["min_viable"][crit]
            _check(
                f"{m_i.key} minimum viable eta_u ({crit})",
                x is not None and abs(x - tgt) <= 1e-3,
                f"{x:.5f} (target {tgt})",
            )
        x = sw["min_viable"]["own_chem"]
        ref = sw["references"]["multimode_f0"]
        alloc_ep = nominal_allocation("multimode", m_i)
        _, _, _, w_lo = _close_arch_at_eta_u(m_i, x - 0.005, alloc_ep, "multimode")
        _, _, _, w_hi = _close_arch_at_eta_u(m_i, x + 0.005, alloc_ep, "multimode")
        _check(
            f"{m_i.key} min-viable bracket (own_chem): "
            f"wet(x - 0.005) > ref >= wet(x + 0.005)",
            w_lo > ref >= w_hi,
            f"{w_lo:.2f} > {ref:.2f} >= {w_hi:.2f} kg",
        )

    sw_mmm_a = sweep_eta_u(mission, arch_key="mono_multimode")
    _check(
        "LEO mono_multimode grid: 41/41 converge, none excluded",
        all(r["converged"] for r in sw_mmm_a["records"])
        and not any(r["excluded"] for r in sw_mmm_a["records"]),
    )
    _check(
        f"LEO mono_multimode own-chem (160 s) reference wet ({_MONO_F0_REF_LEO})",
        abs(sw_mmm_a["references"]["own_f0"] - _MONO_F0_REF_LEO) <= _ATOL,
        f"{sw_mmm_a['references']['own_f0']:.4f} kg",
    )
    rec61 = next(r for r in sw_mmm_a["records"] if r["eta_u"] == 0.61)
    _check(
        "LEO mono_multimode parity grid point == the nominal "
        f"design ({_LEO_TARGETS['mono_multimode'][0]})",
        abs(rec61["wet"] - _LEO_TARGETS["mono_multimode"][0]) <= _ATOL,
        f"{rec61['wet']:.4f} kg",
    )
    for crit, tgt in _MONO_MIN_VIABLE_LEO.items():
        x = sw_mmm_a["min_viable"][crit]
        _check(
            f"LEO mono_multimode minimum viable eta_u ({crit}) — reference {tgt}",
            x is not None and abs(x - tgt) <= 2e-4,
            f"{x:.5f}",
        )
    _check(
        "LEO crossover ordering: mono_multimode thresholds sit "
        "below the biprop multimode's (worse chemical alternative)",
        sw_mmm_a["min_viable"]["own_chem"] < _MIN_VIABLE_TARGETS["leo"]["own_chem"]
        and sw_mmm_a["min_viable"]["single_chem"]
        < _MIN_VIABLE_TARGETS["leo"]["single_chem"],
    )
    sw_mmm_b = sweep_eta_u(mission_b, arch_key="mono_multimode")
    _div_mmm_b = [r["eta_u"] for r in sw_mmm_b["records"] if r["diverged"]]
    _check(
        "GEO mono_multimode eta_u grid: divergence confined to "
        f"eta_u <= {_MONO_GEO_ETA_U_DIVERGED[-1]}",
        _div_mmm_b == _MONO_GEO_ETA_U_DIVERGED,
        f"diverged at {_div_mmm_b}",
    )
    _check(
        "GEO mono_multimode min-viable readings None (own-chem "
        "reference has no fixed point)",
        sw_mmm_b["min_viable"]["own_chem"] is None
        and sw_mmm_b["min_viable"]["single_chem"] is None,
    )
    _check(
        "default eta_u path untouched: multimode references carry "
        "the multimode_f0 alias == own_f0",
        "multimode_f0" not in sw_mmm_a["references"],
    )

    section("Bus-fraction comparison and tank type")
    ok_a = all(
        assess(k, mission, runs[k][0], runs[k][1], envelope_frac=f).envelope_ok
        for k in ("single_elec", "multimode", "hybrid", "single_chem")
        for f in (0.35, 0.50, 1.00)
    )
    _check("LEO envelope never binds over {0.35, 0.5, 1.0}", ok_a)
    ok_b = True
    for f, expected in _GEO_ENVELOPE_TARGETS.items():
        for k, excl in expected.items():
            rep = assess(k, mission_b, runs_b[k][0], runs_b[k][1], envelope_frac=f)
            ok_b &= (not rep.envelope_ok) == excl
    _check(
        "GEO envelope verdicts over {0.35, 0.5, 1.0} match the reference table",
        ok_b,
        "all four excluded at 0.35 and 0.5; single_elec passes at 1.0",
    )
    wets_copv = {}
    for key, tgt in _GEO_COPV_TARGETS.items():
        r_c = close(
            key,
            mission_b,
            nominal_allocation(key, mission_b),
            tank_tech_overrides={"n2o": "copv"},
        )
        wets_copv[key] = r_c.wet
        _check(
            f"GEO COPV {key} wet",
            abs(r_c.wet - tgt) <= _ATOL,
            f"{r_c.wet:.3f} kg (target {tgt})",
        )
        if key == "multimode":
            _check(
                "GEO COPV multimode n_ep == 2",
                r_c.n_units["ep"] == 2,
                f"got {r_c.n_units['ep']}",
            )
    gap_c = wets_copv["multimode"] - wets_copv["hybrid"]
    _check(
        f"GEO COPV multimode - hybrid gap ({_GEO_COPV_GAP:+.2f} kg)",
        abs(gap_c - _GEO_COPV_GAP) <= 0.1,
        f"{gap_c:+.2f} kg",
    )
    _check(
        "GEO COPV ordering: multimode < hybrid",
        wets_copv["multimode"] < wets_copv["hybrid"],
    )

    from ..spacecraft.parameters import override

    ok_noop = True
    ok_equiv = True
    for key in ("single_elec", "multimode", "hybrid", "single_chem"):
        alloc_k = nominal_allocation(key, mission_b)
        w_def = close(key, mission_b, alloc_k).wet
        with override(chem_tank_technology="metal"):
            w_met = close(key, mission_b, alloc_k).wet
        with override(chem_tank_technology="copv"):
            w_tog = close(key, mission_b, alloc_k).wet
        w_ovr = close(
            key,
            mission_b,
            alloc_k,
            tank_tech_overrides={"n2o": "copv", "propene": "copv"},
        ).wet
        ok_noop = ok_noop and (w_met == w_def)
        ok_equiv = ok_equiv and (w_tog == w_ovr)
    _check("toggle at metal is a bit-identical no-op (default path unchanged)", ok_noop)
    _check(
        "toggle copv == overrides {'n2o','propene'} -> 'copv' (mechanism identity)",
        ok_equiv,
    )

    section("Warning-time sweep")
    sw_tw = sweep_t_warn(ARCHITECTURE_ORDER)
    grid = sw_tw["grid"]
    idx4 = [grid.index(tw) for tw in T_WARN_SET_S]
    for key, tgt in _T_WARN_VERDICT_TARGETS.items():
        got = tuple(sw_tw["records"][key][i]["feasible"] for i in idx4)
        _check(
            f"t_warn verdicts {key} over {{2 h, 6 h, 1 d, 3 d}}",
            got == tgt,
            f"got {got}",
        )
    _check(
        f"single_elec CA flip within 0.5 h of {_T_WARN_FLIP_S / 3600:.2f} h",
        sw_tw["flip_s"] is not None
        and abs(sw_tw["flip_s"] - _T_WARN_FLIP_S) <= _T_WARN_FLIP_ATOL_S,
        f"{sw_tw['flip_s'] / 3600:.3f} h (target {_T_WARN_FLIP_S / 3600:.3f} h)",
    )
    i6 = grid.index(21600.0)
    w6 = sw_tw["records"]["multimode"][i6]["wet"]
    _check(
        "multimode 6 h spot wet",
        abs(w6 - _T_WARN_MM_6H_WET) <= _ATOL,
        f"{w6:.4f} kg (target {_T_WARN_MM_6H_WET})",
    )
    ok_mono = True
    for key, recs in sw_tw["records"].items():
        seq = [r for r in recs if key != "single_elec" or r["drift_ok"]]
        ok_mono &= all(
            a["ca_prop"] >= b["ca_prop"] - 1e-12 for a, b in zip(seq, seq[1:])
        )
    _check(
        "CA-class propellant monotone in 1/t_warn (single_elec: "
        "within the drift_ok region)",
        ok_mono,
    )
    _check(
        "every t_warn grid closure converged",
        all(r["converged"] for recs in sw_tw["records"].values() for r in recs),
    )
    _check(
        "a wet-mass-optimal feasible contender exists at every grid point",
        all(o["arch"] is not None for o in sw_tw["optimal"]),
    )
    opt4 = tuple(
        next(o for o in sw_tw["optimal"] if o["t_warn_s"] == tw)["arch"]
        for tw in T_WARN_SET_S
    )
    _check(
        f"optimal contender over {{2 h, 6 h, 1 d, 3 d}} == {_T_WARN_OPTIMAL_TARGETS}",
        opt4 == _T_WARN_OPTIMAL_TARGETS,
        f"got {opt4}",
    )

    section("Transfer-window sweep")
    for t_yr, tgt in _WINDOW_CHECK_TARGETS.items():
        q = transfer_asymmetry_check(t_max_yr=t_yr)
        for name in ("hall", "ecr"):
            n_t, p_t, eps_t, el_t, ok_t = tgt[name]
            d = q[name]
            _check(
                f"closed form {name} at {t_yr * 12:.0f} mo (n/P/EPS/elec/power_ok)",
                d["n_units"] == n_t
                and abs(d["p_dc_w"] / 1e3 - p_t) < 0.01
                and abs(d["eps_inc"] - eps_t) < 0.1
                and abs(d["electronics"] - el_t) < 0.1
                and d["power_ok"] == ok_t,
                f"n={d['n_units']} {d['p_dc_w'] / 1e3:.3f} kW "
                f"EPS {d['eps_inc']:.2f} kg ok={d['power_ok']}",
            )

    sw_t2 = sweep_t_transfer_max()
    months = tuple(round(t * 12, 2) for t in sw_t2["t_set_yr"])
    _check(
        "transfer-window set {1 wk, 2 wk, 1 mo, 3, 6, 12, 24 mo}",
        months == (0.23, 0.46, 1.0, 3.0, 6.0, 12.0, 24.0),
        f"got {months}",
    )

    recs_se = sw_t2["records"]["single_elec"]
    n_se = tuple(None if r["diverged"] else r["n_ep"] for r in recs_se)
    _check(
        "single_elec nominal n_Hall over the set == "
        f"{_WINDOW_SE_N_HALL} (None = diverged)",
        n_se == _WINDOW_SE_N_HALL,
        f"got {n_se}",
    )
    for key, tgt in _WINDOW_VERDICT_TARGETS.items():
        got = tuple(r["feasible"] for r in sw_t2["records"][key])
        _check(
            f"nominal T_max verdicts {key} over the 7-point set",
            got == tgt,
            f"got {got}",
        )
    _check(
        "multimode / hybrid nominal wets are window-invariant "
        f"(chemical transfer): {_WINDOW_MM_HY_WET[0]} / "
        f"{_WINDOW_MM_HY_WET[1]} kg at every point",
        all(
            abs(r["wet"] - _WINDOW_MM_HY_WET[0]) <= _ATOL
            for r in sw_t2["records"]["multimode"]
        )
        and all(
            abs(r["wet"] - _WINDOW_MM_HY_WET[1]) <= _ATOL
            for r in sw_t2["records"]["hybrid"]
        ),
    )

    eptr_mm = sw_t2["ep_transfer"]["multimode"]
    _check(
        "multimode EP-transfer 24 mo split (s* "
        f"{_WINDOW_EPTR_MM_24MO_S}, {_WINDOW_EPTR_MM_24MO_WET} kg)",
        abs(eptr_mm[-1]["s_ep"] - _WINDOW_EPTR_MM_24MO_S) <= 5e-4
        and abs(eptr_mm[-1]["wet"] - _WINDOW_EPTR_MM_24MO_WET) <= 0.05,
        f"s* {eptr_mm[-1]['s_ep']:.3f}, {eptr_mm[-1]['wet']:.2f} kg",
    )
    _check(
        "multimode EP-transfer excluded at every window",
        all(r["excluded"] for r in eptr_mm) == _WINDOW_EPTR_MM_ALL_EXCLUDED,
    )
    eptr_hy = sw_t2["ep_transfer"]["hybrid"]
    by_mo = {round(r["t_yr"] * 12, 2): r for r in eptr_hy}
    ok_hy = all(
        abs(by_mo[m]["wet"] - w) <= _ATOL for m, w in _WINDOW_EPTR_HY_WETS.items()
    )
    _check(
        "hybrid EP-transfer wets at 3/6/12/24 mo == "
        + "/".join(f"{w:.2f}" for w in _WINDOW_EPTR_HY_WETS.values())
        + " kg",
        ok_hy,
    )
    _check(
        "hybrid EP-transfer verdicts over the set == "
        + "("
        + ",".join("T" if v else "F" for v in _WINDOW_EPTR_HY_VERDICTS)
        + ")",
        tuple(r["feasible"] for r in eptr_hy) == _WINDOW_EPTR_HY_VERDICTS,
        f"got {tuple(r['feasible'] for r in eptr_hy)}",
    )

    recs_mmm = sw_t2["records"]["mono_multimode"]
    _check(
        "mono_multimode nominal transfer converges window-invariantly "
        "and is envelope-excluded throughout",
        all(
            not r["diverged"]
            and r["n_ep"] == _WINDOW_MMM_NOMINAL[1]
            and abs(r["wet"] - _WINDOW_MMM_NOMINAL[0]) <= _ATOL
            and not r["feasible"]
            for r in recs_mmm
        ),
        f"wets {sorted({round(r['wet'], 2) for r in recs_mmm})}, "
        f"n_ep {sorted({r['n_ep'] for r in recs_mmm})}",
    )
    eptr_mmm = sw_t2["ep_transfer"]["mono_multimode"]
    _check(
        "mono_multimode EP-transfer: a feasible split exists at "
        "every window, and every one is excluded",
        all(
            not r["no_feasible_s"] and not r["diverged"] and not r["feasible"]
            for r in eptr_mmm
        ),
        f"no_feas = {tuple(r['no_feasible_s'] for r in eptr_mmm)}",
    )
    _check(
        "mono_multimode EP-transfer 24 mo: s* "
        f"{_WINDOW_EPTR_MMM_24MO[0]}, converges, power_ok, n_EP "
        f"{_WINDOW_EPTR_MMM_24MO[1]}, envelope-excluded",
        abs(eptr_mmm[6]["s_ep"] - _WINDOW_EPTR_MMM_24MO[0]) <= 5e-4
        and not eptr_mmm[6]["diverged"]
        and eptr_mmm[6]["power_ok"]
        and eptr_mmm[6]["n_ep"] == _WINDOW_EPTR_MMM_24MO[1]
        and not eptr_mmm[6]["feasible"],
        f"s* {eptr_mmm[6]['s_ep']:.3f}, n_EP {eptr_mmm[6]['n_ep']}, "
        f"wet {eptr_mmm[6]['wet']:.1f} kg",
    )

    sw_off = sweep_t_transfer_max(ep_transfer=False)
    _check("toggle off: ep_transfer reading empty", not sw_off["ep_transfer"])
    ident = all(
        (ra["wet"] == rb["wet"] or (ra["wet"] != ra["wet"] and rb["wet"] != rb["wet"]))
        and ra["n_ep"] == rb["n_ep"]
        for k in sw_off["records"]
        for ra, rb in zip(sw_off["records"][k], sw_t2["records"][k])
    )
    _check("toggle off: nominal records bit-identical to toggle on", ident)

    from ..spacecraft.sweeps import ep_transfer_allocation, ep_max_units
    from ..spacecraft.architectures import ecr_unit as _ecr_u, hall_unit as _hall_u

    m_1mo = geo_mission(t_transfer_max_yr=1.0 / 12.0)
    _check(
        f"n_max under the power-system limit without degradation (d = 0): "
        f"Hall {_SPLIT_N_MAX[0]} | ECR {_SPLIT_N_MAX[1]}",
        ep_max_units(_hall_u(), m_1mo) == _SPLIT_N_MAX[0]
        and ep_max_units(_ecr_u(), m_1mo) == _SPLIT_N_MAX[1],
    )
    from ..spacecraft.mass_models import transfer_array_derate

    _d_full = transfer_array_derate(m_1mo, 1.0)
    _check(
        f"n_max under the power-system limit at the nominal derate (d = "
        f"{_d_full:g}): Hall {_SPLIT_N_MAX_NOMINAL[0]} | ECR "
        f"{_SPLIT_N_MAX_NOMINAL[1]}",
        ep_max_units(_hall_u(), m_1mo, _d_full) == _SPLIT_N_MAX_NOMINAL[0]
        and ep_max_units(_ecr_u(), m_1mo, _d_full) == _SPLIT_N_MAX_NOMINAL[1],
        f"Hall {ep_max_units(_hall_u(), m_1mo, _d_full)} | ECR "
        f"{ep_max_units(_ecr_u(), m_1mo, _d_full)}",
    )
    al_sp, dg_sp = ep_transfer_allocation("hybrid", m_1mo)
    _check(
        f"hybrid 1 mo split share s* (reference {_SPLIT_S})",
        abs(dg_sp["s_ep"] - _SPLIT_S) <= 5e-4,
        f"s* = {dg_sp['s_ep']:.6f}",
    )
    r_sp = close("hybrid", m_1mo, al_sp)
    rep_sp = assess("hybrid", m_1mo, al_sp, r_sp)
    _check(
        f"hybrid 1 mo split wet (reference {_SPLIT_WET})",
        abs(r_sp.wet - _SPLIT_WET) <= _ATOL,
        f"{r_sp.wet:.3f} kg",
    )
    _n_own = ep_max_units(_hall_u(), m_1mo, transfer_array_derate(m_1mo, dg_sp["s_ep"]))
    _check(
        "hybrid 1 mo split n_hall == the unit ceiling at its own derate "
        f"d x s* ({_SPLIT_N_HALL})",
        r_sp.n_units["ep"] == _n_own == _SPLIT_N_HALL,
        f"got {r_sp.n_units['ep']} (ceiling at d x s* {_n_own}; "
        f"all-electric n_max {dg_sp['n_max']})",
    )
    _check(
        f"hybrid 1 mo split EPS (reference {_SPLIT_EPS} kg, within the 30 kg limit)",
        abs(r_sp.eps_inc - _SPLIT_EPS) <= _ATOL and rep_sp.power_ok,
        f"{r_sp.eps_inc:.3f} kg",
    )
    _check(
        "hybrid 1 mo split transfer EP propellant = fleet window "
        f"capacity (reference {_SPLIT_TR_EP} kg)",
        abs(r_sp.per_class_prop["transfer"]["ep"] - _SPLIT_TR_EP) <= _ATOL,
        f"{r_sp.per_class_prop['transfer']['ep']:.3f} kg expelled",
    )
    _check(
        "hybrid 1 mo split: every class verdict feasible (the "
        "top-up rescues the window); exclusion is envelope-only",
        all(v.feasible for v in rep_sp.verdicts)
        and rep_sp.power_ok
        and not rep_sp.envelope_ok,
    )
    m_24 = geo_mission(t_transfer_max_yr=2.0)
    al_24, dg_24 = ep_transfer_allocation("hybrid", m_24)
    r_24 = close("hybrid", m_24, al_24)
    _check(
        "hybrid 24 mo: s* = 1.0 and the closure reproduces the "
        f"full-EP value ({_FULL_EP_24MO_WET:.2f} kg)",
        dg_24["s_ep"] == 1.0 and abs(r_24.wet - _FULL_EP_24MO_WET) <= _ATOL,
        f"s* = {dg_24['s_ep']}, wet {r_24.wet:.4f} kg",
    )
    al_se, dg_se = ep_transfer_allocation("single_elec", m_1mo)
    _check(
        "single_elec toggle allocation == nominal (no top-up path)",
        dg_se["forced_all_ep"] and al_se == nominal_allocation("single_elec", m_1mo),
    )

    section("Sensitivity method")
    from ..spacecraft.parameters import override, REGISTER
    from ..spacecraft.architectures import invalidate_cache

    with override(pf_metal=6350.0):
        invalidate_cache()
        m_a = leo_mission()
        r_pf = close("single_chem", m_a, nominal_allocation("single_chem", m_a))
    invalidate_cache()
    _check(
        f"pf_metal = 6350 m single_chem LEO wet (reference {_TORNADO_PF6350_SC_WET})",
        abs(r_pf.wet - _TORNADO_PF6350_SC_WET) <= _ATOL,
        f"{r_pf.wet:.4f} kg",
    )
    tor = tornado()
    se_base = tor["baseline"]["leo"]["wet_se"]
    se_lo = next(b for b in tor["bars"] if b["name"] == "power_margin_150")["edges"][
        "lo"
    ]
    _check(
        "power_margin_150 = 50 W moves single_elec LEO wet by "
        f"{_TORNADO_MARGIN50_SE_DWET} kg",
        se_lo["error"] is None
        and abs((se_lo["leo"]["wet_se"] - se_base) - _TORNADO_MARGIN50_SE_DWET) <= 0.01,
        f"d = {se_lo['leo']['wet_se'] - se_base:+.4f} kg",
    )
    mm_base = tor["baseline"]["leo"]["wet_mm"]
    mm_hi = next(b for b in tor["bars"] if b["name"] == "m_head_ecr")["edges"]["hi"]
    d_prop = mm_hi["leo"]["wet_mm"] - mm_base
    _check(
        "m_head_ecr override propagates through the architecture cache",
        mm_hi["error"] is None and abs(d_prop - _TORNADO_HEAD_ECR_MM_DWET) <= 0.01,
        f"d = {d_prop:+.4f} kg (target {_TORNADO_HEAD_ECR_MM_DWET})",
    )
    oat_rows = sorted(k for k, p in REGISTER.items() if p.sweep == "oat")
    bar_names = [b["name"] for b in tor["bars"]]
    _check(
        "every 'oat' register row appears exactly once in the tornado output",
        sorted(n for n in bar_names if n in oat_rows) == oat_rows
        and all(bar_names.count(n) == 1 for n in oat_rows),
        f"{len(oat_rows)} oat rows",
    )
    _check(
        "the three thruster band bars are present",
        [n for n in bar_names if n.startswith("thruster_")]
        == ["thruster_R_Isp_ue", "thruster_R_T", "thruster_R_Pin"],
    )
    _check("p_heater_mono is an oat bar in the tornado", "p_heater_mono" in bar_names)
    hb = next(b for b in tor["bars"] if b["name"] == "p_heater_mono")
    base_mmm = tor["baseline"]["leo"]["wet_mmm"]
    hw = max(
        abs(hb["edges"][e]["leo"]["wet_mmm"] - base_mmm)
        for e in ("lo", "hi")
        if not hb["edges"][e]["error"]
    )
    _check(
        "p_heater_mono bar zero-width on mono_mm LEO (10-50 W)",
        hw <= _ATOL,
        f"max |d wet| = {hw:.2e} kg",
    )
    _check(
        "d(mono_mm - mm) tornado baseline == the LEO difference "
        f"({_LEO_DELTA_TARGETS['mmm_minus_mm']:+.2f} kg)",
        abs(tor["baseline"]["leo"]["d_mmm_mm"] - _LEO_DELTA_TARGETS["mmm_minus_mm"])
        <= _ATOL,
        f"{tor['baseline']['leo']['d_mmm_mm']:+.4f} kg",
    )
    _d_mmm_b = tor["baseline"]["geo"]["d_mmm_mm"]
    _check(
        f"d(mono_mm - mm) GEO tornado baseline == {_TORNADO_D_MMM_MM_GEO:+.2f} kg",
        not math.isnan(_d_mmm_b) and abs(_d_mmm_b - _TORNADO_D_MMM_MM_GEO) <= 0.05,
        f"{_d_mmm_b:+.2f} kg",
    )
    m_a2 = leo_mission()
    r_post = close("multimode", m_a2, nominal_allocation("multimode", m_a2))
    _check(
        "parameter values restored after the tornado "
        "(nominal multimode reproduces "
        f"{_LEO_TARGETS['multimode'][0]})",
        abs(r_post.wet - _LEO_TARGETS["multimode"][0]) <= _ATOL,
        f"{r_post.wet:.4f} kg",
    )

    cg = mono_corner_grid()
    lv_ok = (
        all(
            (lv["lo"], lv["nom"], lv["hi"])
            == (REGISTER[key].range[0], REGISTER[key].value, REGISTER[key].range[1])
            for key, lv in cg["levels"].items()
        )
        and set(cg["levels"]) == {"m_head_mono", "isp_mono"}
        and len(cg["records"]) == 9
    )
    _check("corner grid levels read from the register brackets (3x3, both rows)", lv_ok)
    by = {(r["head_level"], r["isp_level"]): r for r in cg["records"]}
    r_nom = by[("nom", "nom")]
    d_bit = max(abs(r_nom["wets"][k] - runs[k][1].wet) for k in ARCHITECTURE_ORDER)
    _check(
        "nominal grid point bit-reproduces the nominal closures "
        "(all five architectures)",
        d_bit <= 1e-9,
        f"max |d wet| = {d_bit:.2e} kg",
    )
    err9 = max(
        abs(by[k]["wets"]["mono_multimode"] - tgt)
        for k, tgt in _CORNER_GRID_TARGETS.items()
    )
    _check(
        "mono_multimode wet matches the reference grid at all 9 points",
        err9 <= _ATOL,
        f"max err {err9:.4f} kg",
    )
    _check(
        "all 45 grid closures converged",
        all(all(r["converged"].values()) for r in cg["records"]),
    )
    spread = max(
        abs(r["wets"][k] - r_nom["wets"][k])
        for r in cg["records"]
        for k in ARCHITECTURE_ORDER
        if k != "mono_multimode"
    )
    _check(
        "the four non-mono architectures are grid-invariant "
        "(the swept rows reach only mono_unit)",
        spread <= 1e-9,
        f"max spread {spread:.2e} kg",
    )
    _check(
        f"unfavourable corner (head hi, isp lo) d(mono_mm - mm) ({_CORNER_D_UNFAV})",
        abs(cg["corner"]["computed"] - _CORNER_D_UNFAV) <= _ATOL,
        f"{cg['corner']['computed']:+.4f} kg",
    )
    _check(
        f"interaction term == affine cross term ({_CORNER_INTERACTION:+g} kg)",
        abs(cg["corner"]["interaction"] - _CORNER_INTERACTION)
        <= _CORNER_INTERACTION_ATOL,
        f"{cg['corner']['interaction']:+.6f} kg",
    )
    _check(
        "unfavourable-corner lead stays negative",
        _CORNER_D_UNFAV + _ATOL < 0.0 and cg["corner"]["computed"] < 0.0,
        f"{cg['corner']['computed']:+.4f} kg",
    )
    _check(
        "mono_multimode lightest of the five at every grid point",
        cg["lightest_everywhere"],
    )
    m_a3 = leo_mission()
    r_post_cg = close("multimode", m_a3, nominal_allocation("multimode", m_a3))
    _check(
        "parameter values restored after the corner "
        "grid (nominal multimode reproduces "
        f"{_LEO_TARGETS['multimode'][0]})",
        abs(r_post_cg.wet - _LEO_TARGETS["multimode"][0]) <= _ATOL,
        f"{r_post_cg.wet:.4f} kg",
    )

    section("Input consistency")
    _tr = mission_b.get_class("transfer")
    _dv_ep_class = _tr.dv_impulsive * _tr.ep_loss_ratio
    _dv_ep_der = value("dv_transfer_chem") / value("eta_planning_transfer")
    _check(
        "GEO transfer class EP budget == dv_transfer_chem / "
        "eta_planning_transfer (|diff| < 1e-9 m/s)",
        abs(_dv_ep_class - _dv_ep_der) < 1e-9
        and abs(dv_transfer_ep() - _dv_ep_der) < 1e-9,
        f"{_dv_ep_class:.4f} m/s vs {_dv_ep_der:.4f} m/s",
    )
    _chem_alt = REGISTER["dv_transfer_chem"].range[0]
    with override(dv_transfer_chem=_chem_alt):
        _tr_alt = geo_mission().get_class("transfer")
        _dv_ep_alt = _tr_alt.dv_impulsive * _tr_alt.ep_loss_ratio
        _q_alt = transfer_asymmetry_check()["T_req_N"]
    _q_nom = transfer_asymmetry_check()["T_req_N"]
    _check(
        "a dv_transfer_chem override carries the EP budget (transfer "
        "class and asymmetry check both scale by chem_alt / chem)",
        abs(_dv_ep_alt / _dv_ep_class - _chem_alt / value("dv_transfer_chem")) < 1e-12
        and abs(_q_alt / _q_nom - _chem_alt / value("dv_transfer_chem")) < 1e-12,
        f"chem {_chem_alt:g} m/s -> EP {_dv_ep_alt:.2f} m/s",
    )

    section("Result")
    if _FAILURES:
        print(f"  {len(_FAILURES)} check(s) FAILED:")
        for f in _FAILURES:
            print(f"    - {f}")
        raise SystemExit(1)
    print("  All spacecraft-model checks passed.")
