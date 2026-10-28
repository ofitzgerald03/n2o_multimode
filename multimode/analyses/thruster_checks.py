"""
Thruster-model checks: reproduction of the Cannat et al. (2015) operating
point, the argon ionisation cost against the fit of Petro and Sedwick
(2017), Eq. (17) against Eq. (18), an analytic limit, and reference values of
the ionisation cost and the N2O operating point (relative tolerance 1e-6).
"""

import numpy as np

from ..config import PARAMS
from ..thruster.pipeline import run_pipeline
from ..thruster.thrust_efficiency import (
    thrust_efficiency,
    thrust_efficiency_from_eq17,
    PS_Ar_curve_fit,
)
from .registry import task


# Reference values (regression targets) at the PARAMS defaults
REFERENCE = {
    # E_c [eV] at Te = 2, 5, 10, 22 eV
    "Ec_Xe": {
        2: 57.199867106963126,
        5: 26.349465539845067,
        10: 20.280051479694965,
        22: 17.20196935766925,
    },
    "Ec_N2O": {
        2: 854.0413313552535,
        5: 65.14204523056067,
        10: 37.05533353831231,
        22: 26.814797720781954,
    },
    "Ec_H2O": {
        2: 436.6346722529455,
        5: 29.948020935749085,
        10: 23.237570376761866,
        22: 22.193075537214007,
    },
    "Ec_Ar": {
        2: 66.13425669646354,
        5: 31.564057651086955,
        10: 25.17788204820193,
        22: 22.03910323486576,
    },
    # Xe at the anchor (Cannat et al. 2015, Table II)
    "Xe_T": 0.98e-3,  # N
    "Xe_mdot": 1.0e-7,  # kg/s
    "Xe_Pin": 30.0,  # W
    "Xe_Isp_gross": 999.3218887183698,  # s
    # N2O at the anchor, common Bohm velocity
    "N2O_B_T": 0.0009802319571572402,  # N
    "N2O_B_mdot": 5.5297898735885645e-08,  # kg/s
    "N2O_B_Pin": 61.30635816490098,  # W
    "N2O_B_Isp_gross": 1807.5884292935332,  # s
    "N2O_B_Isp_int": 2520.2465431342926,  # s
}

RTOL = 1e-6


def _assert_close(name, got, expected, rtol=RTOL, atol=0.0):
    ok = np.isclose(got, expected, rtol=rtol, atol=atol)
    status = "OK  " if ok else "FAIL"
    print(f"  [{status}] {name:<38} got {got:.10g}   expect {expected:.10g}")
    return bool(ok)


@task("check_thruster", "Thruster-model checks", in_all=False)
def run_thruster_checks(ctx):
    print("Running pipeline (Xe, Ar, N2O, H2O) at PARAMS defaults...")
    results = run_pipeline(
        params=PARAMS,
        propellants=("Xe", "Ar", "N2O", "H2O"),
        anchor_propellant="Xe",
    )
    Te = results["Te"]
    checks = []

    print("\nE_c at spot temperatures (reference values, rel 1e-6):")
    for prop in ("Xe", "N2O", "H2O", "Ar"):
        for t, exp in REFERENCE[f"Ec_{prop}"].items():
            i = int(np.argmin(np.abs(Te - t)))
            checks.append(
                _assert_close(f"Ec {prop} @ {t} eV", results[prop]["Ec"][i], exp)
            )

    print("\nCannat 2015 Table II reproduction (Xe at eta_u_xe = 0.61):")
    p_xe = dict(PARAMS)
    p_xe["eta_u_xe"] = p_xe["eta_u_xe_ref"]
    res_xe = run_pipeline(params=p_xe, propellants=("Xe",), anchor_propellant="Xe")
    op = res_xe["Xe"]["op_point_A"]
    checks.append(_assert_close("Xe T [N]", op["T"], REFERENCE["Xe_T"]))
    checks.append(_assert_close("Xe m_dot [kg/s]", op["m_dot"], REFERENCE["Xe_mdot"]))
    checks.append(_assert_close("Xe P_in [W]", op["P_in"], REFERENCE["Xe_Pin"]))
    checks.append(
        _assert_close("Xe Isp_gross [s]", op["Isp_gross"], REFERENCE["Xe_Isp_gross"])
    )

    print("\nXe propellant-physics ratios at the anchor (must be 1.000):")
    r = op["ratios"]
    for k in ("R_T", "R_mdot_ion", "R_ue", "R_Pin"):
        checks.append(_assert_close(f"Xe {k}", r[k], 1.0, rtol=1e-9))

    print("\nN2O at the anchor, common Bohm velocity (eta_u = 0.61):")
    opB = results["N2O"]["op_point_B"]
    checks.append(_assert_close("N2O_B T [N]", opB["T"], REFERENCE["N2O_B_T"]))
    checks.append(
        _assert_close("N2O_B m_dot [kg/s]", opB["m_dot"], REFERENCE["N2O_B_mdot"])
    )
    checks.append(_assert_close("N2O_B P_in [W]", opB["P_in"], REFERENCE["N2O_B_Pin"]))
    checks.append(
        _assert_close(
            "N2O_B Isp_gross [s]", opB["Isp_gross"], REFERENCE["N2O_B_Isp_gross"]
        )
    )
    checks.append(
        _assert_close(
            "N2O_B Isp_intensive [s]", opB["Isp_intensive"], REFERENCE["N2O_B_Isp_int"]
        )
    )

    print("\nEquation checks:")
    xe = results["Xe"]
    max_diff = 0.0
    for R_m in (1.0, 10.0, 100.0):
        eta_18 = thrust_efficiency(Te, xe["phi_sh_A"], xe["Ec"], R_m)
        eta_17 = thrust_efficiency_from_eq17(Te, xe["phi_sh_A"], xe["Ec"], R_m)
        max_diff = max(max_diff, float(np.max(np.abs(eta_18 - eta_17))))
    ok = max_diff < 1e-10
    checks.append(ok)
    print(
        f"  [{'OK  ' if ok else 'FAIL'}] Eq. 17 ≡ 18 max |diff|"
        f"{'':<21} {max_diff:.2e}  (< 1e-10)"
    )

    zero_Ec = np.zeros_like(Te)
    eta_num = thrust_efficiency(Te, xe["phi_sh_A"], zero_Ec, 1.0)
    eta_analytic = xe["phi_sh_A"] / (2.0 * (2.0 * Te + xe["phi_sh_A"]))
    lim_diff = float(np.max(np.abs(eta_num - eta_analytic)))
    ok = lim_diff < 1e-12
    checks.append(ok)
    print(
        f"  [{'OK  ' if ok else 'FAIL'}] Analytic limit (Rm=1, Ec=0)"
        f"{'':<15} {lim_diff:.2e}  (< 1e-12)"
    )

    print(
        "\nAr E_c vs P&S Eqs. (9)-(10) curve fit "
        "(expected: -1.85% @5 eV, -0.32% @10 eV):"
    )
    Ec_ar = results["Ar"]["Ec"]
    Ec_ar_PS = PS_Ar_curve_fit(Te)
    for t, exp_pct, tol in ((5, -1.85, 0.15), (10, -0.32, 0.15)):
        i = int(np.argmin(np.abs(Te - t)))
        diff = 100.0 * (Ec_ar[i] - Ec_ar_PS[i]) / Ec_ar_PS[i]
        ok = abs(diff - exp_pct) < tol
        checks.append(ok)
        print(
            f"  [{'OK  ' if ok else 'FAIL'}] Ar vs P&S @ {t} eV"
            f"{'':<24} {diff:+.2f}%  (expect {exp_pct:+.2f}%)"
        )

    n_fail = sum(1 for c in checks if not c)
    print("\n" + "=" * 78)
    if n_fail == 0:
        print(f"CHECK PASSED — all {len(checks)} checks OK.")
    else:
        print(f"CHECK FAILED — {n_fail}/{len(checks)} checks failed.")
    print("=" * 78)
    if n_fail:
        raise SystemExit(1)
