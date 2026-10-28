"""
Thruster-level results: anchor operating points, the N2O to Xe ratio table
and the comparison of the two multi-species Bohm closures.
"""

import csv

import numpy as np

from ..config import PARAMS
from ..thruster.anchoring import (
    at_anchor_operating_point,
    compose_extensive_curves,
    framework_extensive_arrays,
)
from ..thruster.cross_sections import compute_ion_fractions_species_resolved
from ..thruster.exhaust import (
    bohm_velocity_individual,
    bohm_velocity_system,
    compute_fractions_for_propellant,
    compute_per_propellant_quantities,
    effective_sheath_potential_closure_A,
)
from ..thruster.pipeline import _build_anchor_reference, run_pipeline
from ..thruster.propellants import ALL_PROPELLANTS
from ..thruster.thrust_efficiency import thrust_efficiency
from .common import csv_dir, section
from .registry import task


# Anchor operating points


def _print_anchor_reference(anchor_ref):
    print(f"  T_meas      = {anchor_ref.T_meas * 1e3:.3f} mN")
    print(f"  mdot_meas   = {anchor_ref.mdot_meas * 1e6:.3f} mg/s")
    print(f"  Pin_meas    = {anchor_ref.Pin_meas:.2f} W")
    print(f"  Isp_meas    = {anchor_ref.Isp_meas:.1f} s")
    print(f"  eta_u_meas  = {anchor_ref.eta_u_meas:.3f}")
    print(f"  Te_ref      = {anchor_ref.Te_ref:.2f} eV")
    print(f"  Rm_ref      = {anchor_ref.Rm_ref:.2f}")


def _print_op_point(label, op):
    print(f"  {label}:")
    print(f"    T             = {op['T'] * 1e3:.4f} mN")
    print(f"    m_dot_total   = {op['m_dot'] * 1e6:.4f} mg/s")
    print(f"    m_dot_ion     = {op['m_dot_ion'] * 1e6:.4f} mg/s")
    print(f"    P_in          = {op['P_in']:.3f} W")
    print(f"    Isp_gross     = {op['Isp_gross']:.1f} s")
    print(f"    Isp_intensive = {op['Isp_intensive']:.1f} s")
    r = op["ratios"]
    print(
        f"    ratios (vs Xe): "
        f"R_T={r['R_T']:.4f}  R_mdot_ion={r['R_mdot_ion']:.4f}  "
        f"R_ue={r['R_ue']:.4f}  R_Pin={r['R_Pin']:.4f}"
    )


def _comparison_table(rows, header_labels):
    cols = ["Quantity"] + header_labels
    widths = [22] + [max(len(label), 14) for label in header_labels]

    def fmt_row(cells):
        return "  " + "".join(f"{c:<{w}}" for c, w in zip(cells, widths))

    print(fmt_row(cols))
    print("  " + "-" * (sum(widths)))
    for r in rows:
        print(fmt_row(r))


@task(
    "thruster_anchor",
    "Anchor reference, Cannat (2015) reproduction, "
    "at-anchor operating points and the eta_u scenarios",
)
def run_thruster_anchor(ctx):
    section("(1) Anchor reference")
    p_default = dict(PARAMS)
    results = run_pipeline(
        params=p_default,
        propellants=("Xe", "N2O"),
        anchor_propellant="Xe",
    )
    anchor_ref = results["anchor"]["reference"]
    _print_anchor_reference(anchor_ref)
    i_anchor = results["anchor"]["i_anchor"]
    print(
        f"\n  Te grid resolution at the anchor: "
        f"Te[i_anchor={i_anchor}] = {results['Te'][i_anchor]:.3f} eV"
    )

    section("(2) Xe reproduction (eta_u_xe = eta_u_xe_ref = 0.61)")
    p_xe = dict(PARAMS)
    p_xe["eta_u_xe"] = p_xe["eta_u_xe_ref"]
    results_xe = run_pipeline(
        params=p_xe,
        propellants=("Xe",),
        anchor_propellant="Xe",
    )
    op_xe = results_xe["Xe"]["op_point_A"]
    rows = [
        (
            "Thrust (mN)",
            f"{anchor_ref.T_meas * 1e3:.3f}",
            f"{op_xe['T'] * 1e3:.4f}",
            "OK" if abs(op_xe["T"] * 1e3 - anchor_ref.T_meas * 1e3) < 1e-3 else "CHECK",
        ),
        (
            "Mass flow (mg/s)",
            f"{anchor_ref.mdot_meas * 1e6:.4f}",
            f"{op_xe['m_dot'] * 1e6:.4f}",
            "OK"
            if abs(op_xe["m_dot"] * 1e6 - anchor_ref.mdot_meas * 1e6) < 1e-3
            else "CHECK",
        ),
        (
            "P_in (W)",
            f"{anchor_ref.Pin_meas:.2f}",
            f"{op_xe['P_in']:.4f}",
            "OK" if abs(op_xe["P_in"] - anchor_ref.Pin_meas) < 1e-2 else "CHECK",
        ),
        (
            "Isp_gross (s)",
            f"{anchor_ref.Isp_meas:.1f}",
            f"{op_xe['Isp_gross']:.2f}",
            "OK" if abs(op_xe["Isp_gross"] - anchor_ref.Isp_meas) < 10 else "CHECK",
        ),
    ]
    _comparison_table(rows, ["Cannat 2015", "Composed", "Match?"])
    print("\n  Isp_gross = T / (mdot g0) from the same Table II row.")

    section(
        "(3) At-anchor operating points and propellant-physics ratios "
        "(eta_u from PARAMS)"
    )
    print(f"  eta_u_xe  = {PARAMS['eta_u_xe']}   eta_u_n2o = {PARAMS['eta_u_n2o']}")
    print()

    for label, case_key, case_name in (
        ("Xe", "op_point_A", "single species"),
        ("N2O", "op_point_A", "species-resolved Bohm velocity"),
        ("N2O", "op_point_B", "common Bohm velocity"),
    ):
        op = results[label][case_key]
        if op is None:
            continue
        _print_op_point(f"{label}, {case_name}", op)
        print()

    section("(4) N2O at the anchor over eta_u scenarios")

    def _run_and_get_n2o(eu_xe, eu_n2o):
        p = dict(PARAMS)
        p["eta_u_xe"] = eu_xe
        p["eta_u_n2o"] = eu_n2o
        r = run_pipeline(
            params=p,
            propellants=("Xe", "N2O"),
            anchor_propellant="Xe",
        )
        return r["Xe"]["op_point_A"], r["N2O"]["op_point_B"]

    scenarios = [
        ("Cannat-reproduction", 0.61, 0.61),
        ("Symmetric (eta_u=1)", 1.00, 1.00),
        ("Optimistic N2O", 0.61, 0.80),
        ("Optimistic N2O", 0.61, 0.50),
        ("Plausible N2O", 0.61, 0.45),
        ("Plausible N2O", 0.61, 0.40),
        ("Plausible N2O", 0.61, 0.35),
        ("Plausible N2O", 0.61, 0.30),
        ("Pessimistic N2O", 0.61, 0.10),
    ]

    print(
        f"  {'Scenario':<22} {'(eu_xe, eu_n2o)':<18} "
        f"{'Isp_N2O (s)':<14} {'mdot_N2O (mg/s)':<18} "
        f"{'Pin_N2O (W)':<14} {'Pin/T_N2O':<12}"
    )
    print("  " + "-" * 96)
    for name, eu_xe, eu_n2o in scenarios:
        op_xe_i, op_n2o_i = _run_and_get_n2o(eu_xe, eu_n2o)
        print(
            f"  {name:<22} ({eu_xe:.2f}, {eu_n2o:.2f})        "
            f"{op_n2o_i['Isp_gross']:<14.0f} "
            f"{op_n2o_i['m_dot'] * 1e6:<18.4f} "
            f"{op_n2o_i['P_in']:<14.2f} "
            f"{op_n2o_i['P_in'] / (op_n2o_i['T'] * 1e3):<12.2f}"
        )


# N2O to Xe ratio table at the spot electron temperatures


def build_ratio_table(results: dict) -> list[dict]:
    """N2O and Xe values and their ratios at each params['spot_Te'].

    The N2O columns use params['n2o_case']; the *_band_% columns give the
    relative difference between the two Bohm closures.
    """
    Te = results["Te"]
    params = results["params"]
    case = params["n2o_case"]

    case_q_key = "q_A" if case == "A" else "q_B"
    case_ext_key = "ext_A" if case == "A" else "ext_B"
    case_eta_key = "eta_T_A" if case == "A" else "eta_T_B"

    xe = results["Xe"]
    n2o = results["N2O"]

    xe_qA = xe["q_A"]
    n2o_q_pri = n2o[case_q_key]
    n2o_qA = n2o["q_A"]
    n2o_qB = n2o["q_B"]

    xe_ext = xe["ext_A"]
    n2o_ext_pri = n2o[case_ext_key]
    n2o_ext_A = n2o["ext_A"]
    n2o_ext_B = n2o["ext_B"]

    eta_xe = xe["eta_T_A"]
    eta_n2o_pri = n2o[case_eta_key]
    eta_n2o_A = n2o["eta_T_A"]
    eta_n2o_B = n2o["eta_T_B"]

    def _bandpct(a, b):
        mid = 0.5 * (a + b)
        return 100.0 * abs(a - b) / mid if mid > 0 else 0.0

    rows = []
    for t in params["spot_Te"]:
        i = int(np.argmin(np.abs(Te - t)))

        Isp_xe = xe_qA["I_sp"][i]
        Isp_n2o = n2o_q_pri["I_sp"][i]
        Ec_xe = xe["Ec"][i]
        Ec_n2o = n2o["Ec"][i]
        eta_xe_i = eta_xe[i]
        eta_n2o_i = eta_n2o_pri[i]

        T_xe = xe_ext["T"][i]
        T_n2o = n2o_ext_pri["T"][i]
        mdot_xe = xe_ext["m_dot"][i]
        mdot_n2o = n2o_ext_pri["m_dot"][i]
        Pin_xe = xe_ext["P_in"][i]
        Pin_n2o = n2o_ext_pri["P_in"][i]

        rows.append(
            {
                "Te_eV": float(Te[i]),
                "n2o_case": case,
                "Isp_Xe_s": float(Isp_xe),
                "Isp_N2O_s": float(Isp_n2o),
                "Isp_ratio": float(Isp_n2o / Isp_xe),
                "Ec_Xe_eV": float(Ec_xe),
                "Ec_N2O_eV": float(Ec_n2o),
                "Ec_ratio": float(Ec_n2o / Ec_xe),
                "etaT_Xe": float(eta_xe_i),
                "etaT_N2O": float(eta_n2o_i),
                "etaT_ratio": (
                    float(eta_n2o_i / eta_xe_i) if eta_xe_i > 0 else float("nan")
                ),
                "T_Xe_mN": float(T_xe * 1e3),
                "T_N2O_mN": float(T_n2o * 1e3),
                "T_ratio": (float(T_n2o / T_xe) if T_xe > 0 else float("nan")),
                "mdot_Xe_mgs": float(mdot_xe * 1e6),
                "mdot_N2O_mgs": float(mdot_n2o * 1e6),
                "mdot_ratio": (
                    float(mdot_n2o / mdot_xe) if mdot_xe > 0 else float("nan")
                ),
                "Pin_Xe_W": float(Pin_xe),
                "Pin_N2O_W": float(Pin_n2o),
                "Pin_ratio": (float(Pin_n2o / Pin_xe) if Pin_xe > 0 else float("nan")),
                "Isp_N2O_band_%": float(_bandpct(n2o_qA["I_sp"][i], n2o_qB["I_sp"][i])),
                "etaT_N2O_band_%": float(_bandpct(eta_n2o_A[i], eta_n2o_B[i])),
                "T_N2O_band_%": float(_bandpct(n2o_ext_A["T"][i], n2o_ext_B["T"][i])),
                "Pin_N2O_band_%": float(
                    _bandpct(n2o_ext_A["P_in"][i], n2o_ext_B["P_in"][i])
                ),
            }
        )
    return rows


def print_ratio_table(rows: list[dict], params: dict) -> None:
    print("=" * 112)
    print(
        f"Ratio table — N2O vs Xe "
        f"(R_m = {params['Rm_nominal']:g}, "
        f"T_ref = {params['T_ref'] * 1e3:g} mN @ "
        f"Te_ref = {params['Te_ref']:g} eV)"
    )
    print("=" * 112)
    hdr = (
        f"{'Te':>4} | "
        f"{'Isp_Xe':>7} {'Isp_N2O':>8} {'ratio':>6} | "
        f"{'Ec_Xe':>7} {'Ec_N2O':>8} {'ratio':>6} | "
        f"{'etaT_Xe':>7} {'etaT_N2O':>8} {'ratio':>6} | "
        f"{'T_Xe':>6} {'T_N2O':>7} {'P_Xe':>7} {'P_N2O':>7}"
    )
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        print(
            f"{r['Te_eV']:>4.1f} | "
            f"{r['Isp_Xe_s']:>7.0f} {r['Isp_N2O_s']:>8.0f} "
            f"{r['Isp_ratio']:>6.2f} | "
            f"{r['Ec_Xe_eV']:>7.2f} {r['Ec_N2O_eV']:>8.2f} "
            f"{r['Ec_ratio']:>6.2f} | "
            f"{r['etaT_Xe']:>7.3f} {r['etaT_N2O']:>8.3f} "
            f"{r['etaT_ratio']:>6.2f} | "
            f"{r['T_Xe_mN']:>6.2f} {r['T_N2O_mN']:>7.2f} "
            f"{r['Pin_Xe_W']:>7.1f} {r['Pin_N2O_W']:>7.1f}"
        )
    print()


def write_ratio_table_csv(rows: list[dict], path) -> None:
    fieldnames = list(rows[0].keys())
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"Saved: {path}")


@task(
    "ratio_table",
    "N2O to Xe ratio table at the spot electron temperatures (ratio_table.csv)",
)
def run_ratio_table(ctx):
    results = run_pipeline(
        params=PARAMS,
        propellants=("Xe", "Ar", "N2O", "H2O"),
        anchor_propellant="Xe",
    )
    rows = build_ratio_table(results)
    print_ratio_table(rows, results["params"])
    if ctx.save_csvs:
        write_ratio_table_csv(rows, csv_dir(ctx) / "ratio_table.csv")


# Multi-species Bohm closures for N2O
# Species-resolved closure: each ion species leaves at its own Bohm
# velocity, with the sqrt(m)-weighted ion fractions this implies.
# Common closure: one system Bohm velocity, kappa-ratio fractions.


def build_closure_A(results, params, fractions):
    """N2O operating point with species-resolved Bohm velocities and the
    given ion fractions, composed against the anchor of `results`."""
    Te = results["Te"]
    n2o = ALL_PROPELLANTS["N2O"]
    eta_RF, eta_D = params["eta_RF"], params["eta_D"]
    Rm_ref = params["Rm_ref"]
    eta_u_n2o = params.get("eta_u_n2o", 1.0)
    frac_A = fractions

    q_A = compute_per_propellant_quantities(
        Te, "multi", fractions=frac_A, species_masses_kg=n2o.species_masses_kg, case="A"
    )
    phi_A = effective_sheath_potential_closure_A(Te, frac_A, n2o.species_masses_kg)
    Ec = results["N2O"]["Ec"]
    eta_T_A = thrust_efficiency(Te, phi_A, Ec, Rm_ref)

    scale_factors = results["anchor"]["scale_factors"]
    i_anchor = results["anchor"]["i_anchor"]
    anchor_ref = _build_anchor_reference(params)
    fw_xe = framework_extensive_arrays(
        results["Xe"]["q_A"], results["Xe"]["eta_T_A_at_Rm_ref"], eta_RF, eta_D
    )
    fw_A = framework_extensive_arrays(q_A, eta_T_A, eta_RF, eta_D)
    composed_A = compose_extensive_curves(fw_A, scale_factors, eta_u_n2o)
    op_A = at_anchor_operating_point(
        composed_A, fw_xe, fw_A, anchor_ref, eta_u_n2o, i_anchor
    )

    return {
        "q_A": q_A,
        "eta_T_A": eta_T_A,
        "phi_A": phi_A,
        "ext_A": composed_A,
        "op_point_A": op_A,
        "frac_A": frac_A,
        "masses_kg": n2o.species_masses_kg,
    }


def vb_scalar(results, consistent, te0):
    """Mean species-resolved Bohm velocity over the system Bohm velocity
    at Te = te0."""
    Te = results["Te"]
    i = int(np.argmin(np.abs(Te - te0)))
    frac_A = consistent["frac_A"]
    frac_B = results["N2O"]["fractions"]
    masses_kg = consistent["masses_kg"]
    vb_eff_A = sum(
        frac_A[sp][i] * bohm_velocity_individual(Te, masses_kg[sp])[i] for sp in frac_A
    )
    v_sys = bohm_velocity_system(Te, frac_B, masses_kg)[i]
    return vb_eff_A / v_sys


def flux_distribution_gap(results, consistent, te0):
    """Largest difference between the per-species mass-flux shares of the
    two closures at Te = te0."""
    Te = results["Te"]
    i = int(np.argmin(np.abs(Te - te0)))
    frac_A, frac_B = consistent["frac_A"], results["N2O"]["fractions"]
    masses_kg = consistent["masses_kg"]
    sp_list = list(frac_A.keys())
    vbi = {sp: bohm_velocity_individual(Te, masses_kg[sp])[i] for sp in sp_list}
    v_sys = bohm_velocity_system(Te, frac_B, masses_kg)[i]
    pfA = np.array([frac_A[sp][i] * vbi[sp] for sp in sp_list])
    pfA /= pfA.sum()
    pfB = np.array([frac_B[sp][i] * v_sys for sp in sp_list])
    pfB /= pfB.sum()
    m = np.array([masses_kg[sp] for sp in sp_list])
    mfA = (pfA * m) / (pfA * m).sum()
    mfB = (pfB * m) / (pfB * m).sum()
    return float(np.abs(mfA - mfB).max())


def fmt_spread(a, b):
    """Symmetric relative spread |b - a| / mean [%]."""
    return 100.0 * abs(b - a) / (0.5 * (a + b))


@task(
    "closure_spread",
    "N2O species-resolved vs common Bohm closure: "
    "intensive and extensive differences (closure_spread.csv)",
)
def run_closure_spread(ctx):
    params = dict(PARAMS)
    results = run_pipeline(
        params, propellants=("Xe", "Ar", "N2O"), anchor_propellant="Xe"
    )

    n2o = ALL_PROPELLANTS["N2O"]
    _, kappas = compute_fractions_for_propellant(
        n2o.ion_channels, n2o.ion_energies, results["Te"]
    )
    masses_amu = {sp: n2o.ion_channels[sp]["mass_amu"] for sp in n2o.ion_channels}
    frac_consistent = compute_ion_fractions_species_resolved(kappas, masses_amu)

    consistent = build_closure_A(results, params, frac_consistent)

    opB = results["N2O"]["op_point_B"]
    opAc = consistent["op_point_A"]

    check = {
        "Isp_intensive": (opAc["Isp_intensive"], 2520.2, 1.0),
        "R_ue": (opAc["ratios"]["R_ue"], 1.809, 0.01),
        "m_dot": (opAc["m_dot"] * 1e6, 0.0520, 0.001),  # mg/s
        "R_T": (opAc["ratios"]["R_T"], 0.940, 0.005),
        "R_Pin": (opAc["ratios"]["R_Pin"], 1.832, 0.01),
    }
    print("Species-resolved closure at the anchor:")
    ok = True
    for k, (got, exp, tol) in check.items():
        good = abs(got - exp) <= tol
        ok = ok and good
        print(
            f"  {k:14s} got {got:10.4f}  expect ~{exp:<8} {'OK' if good else 'MISMATCH'}"
        )
    print(f"  -> {'PASS' if ok else 'FAIL'}\n")

    print("Intensive quantities, species-resolved vs common closure")
    Te = results["Te"]
    mask = (Te >= 5.0) & (Te <= 25.0)
    Ec = results["N2O"]["Ec"]
    intensive = {
        "E_c [eV]": (Ec, Ec),
        "eta_T (Rm_ref)": (consistent["eta_T_A"], results["N2O"]["eta_T_B_at_Rm_ref"]),
        "u_e [m/s]": (consistent["q_A"]["u_e"], results["N2O"]["q_B"]["u_e"]),
        "I_sp,int [s]": (consistent["q_A"]["I_sp"], results["N2O"]["q_B"]["I_sp"]),
    }
    for name, (a, b) in intensive.items():
        maxgap = float(np.nanmax(np.abs(np.asarray(a)[mask] - np.asarray(b)[mask])))
        rel = maxgap / (float(np.nanmean(np.abs(np.asarray(b)[mask]))) + 1e-30)
        print(
            f"  {name:16s} max|resolved - common| over [5,25] eV = {maxgap:.3e}  (rel {rel:.1e})"
        )
    print()

    print(
        "Mean resolved Bohm velocity over the system Bohm velocity (scales T, mdot and P_in)"
    )
    for te0 in (5, 10, 22, 25):
        s = vb_scalar(results, consistent, te0)
        print(
            f"  Te={te0:2d} eV: v_b ratio = {s:.4f}  "
            f"-> extensive spread = {100 * (1 - s):+.1f}%"
        )
    print()

    print("Per-species mass-flux distribution, species-resolved vs common")
    for te0 in (5, 22, 25):
        g = flux_distribution_gap(results, consistent, te0)
        print(f"  Te={te0:2d} eV: max|resolved - common| = {g:.2e}")
    print()

    print("At-anchor operating point (N2O)")
    hdr = f"  {'quantity':14s} {'resolved':>14s} {'common':>12s}"
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    rows = [
        ("T [mN]", opAc["T"] * 1e3, opB["T"] * 1e3),
        ("mdot_ion [mg/s]", opAc["m_dot_ion"] * 1e6, opB["m_dot_ion"] * 1e6),
        ("mdot_tot [mg/s]", opAc["m_dot"] * 1e6, opB["m_dot"] * 1e6),
        ("P_in [W]", opAc["P_in"], opB["P_in"]),
        ("Isp_gross [s]", opAc["Isp_gross"], opB["Isp_gross"]),
        ("Isp_int [s]", opAc["Isp_intensive"], opB["Isp_intensive"]),
        ("R_T", opAc["ratios"]["R_T"], opB["ratios"]["R_T"]),
        ("R_mdot_ion", opAc["ratios"]["R_mdot_ion"], opB["ratios"]["R_mdot_ion"]),
        ("R_ue", opAc["ratios"]["R_ue"], opB["ratios"]["R_ue"]),
        ("R_Pin", opAc["ratios"]["R_Pin"], opB["ratios"]["R_Pin"]),
    ]
    for name, a_c, b in rows:
        print(f"  {name:14s} {a_c:14.4f} {b:12.4f}")
    print()
    print("  Spreads |common - resolved| / mean [%]:")
    print(
        f"    intensive (Isp_int, R_ue): {fmt_spread(opAc['Isp_intensive'], opB['Isp_intensive']):.2f} %, "
        f"{fmt_spread(opAc['ratios']['R_ue'], opB['ratios']['R_ue']):.2f} %"
    )
    print(
        f"    extensive (mdot_tot, P_in): {fmt_spread(opAc['m_dot'], opB['m_dot']):.2f} %, "
        f"{fmt_spread(opAc['P_in'], opB['P_in']):.2f} %"
    )

    if ctx.save_csvs:
        csv_path = csv_dir(ctx) / "closure_spread.csv"
        with open(csv_path, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(
                [
                    "Te_eV",
                    "vb_ratio_resolved_over_common",
                    "extensive_spread_pct",
                    "massflux_dist_gap",
                ]
            )
            for te0 in np.arange(5.0, 25.0 + 0.5, 0.5):
                s = vb_scalar(results, consistent, te0)
                g = flux_distribution_gap(results, consistent, te0)
                w.writerow(
                    [f"{te0:.1f}", f"{s:.6f}", f"{100 * (1 - s):.4f}", f"{g:.3e}"]
                )
        print(f"\nWrote {csv_path}")
