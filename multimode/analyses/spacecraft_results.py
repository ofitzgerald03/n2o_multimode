"""
Spacecraft-level results: the parameter register, the nominal designs of
the LEO and GEO missions with their sensitivities, and the sweeps over
warning time, transfer window, eta_u, planning efficiency and the
one-at-a-time inputs. Tables are printed and written as CSV.
"""

import pandas as pd

from ..spacecraft.allocation import nominal_allocation, derived_f_ep
from ..spacecraft.architectures import ARCHITECTURE_ORDER, get_architecture
from ..spacecraft.closure import close, expected_idle_offsets
from ..spacecraft.feasibility import assess, transfer_asymmetry_check
from ..spacecraft.mass_models import copv_pf, floor_crossover
from ..spacecraft.mission import leo_mission, geo_mission, ca_dv_per_event
from ..spacecraft.parameters import (
    REGISTER,
    export_register_csv,
    override,
    value,
    dv_transfer_ep,
    dv_transfer_ep_range,
)
from ..spacecraft.sweeps import (
    sweep_eta_u,
    sweep_t_warn,
    sweep_t_transfer_max,
    tornado,
    ca_flip_band,
    ep_transfer_allocation,
    sweep_eta_planning,
    eta_planning_crossover,
    ETA_PLANNING_SET,
    STRATEGY_LABEL,
    T_WARN_SET_S,
)
from .common import csv_dir, section
from .registry import task

# Measured water utilisation range, Sheppard and Little (2022)
MEASURED_H2O_ETA_U = (0.05, 0.20)


def nominal_designs_leo(tank_tech_overrides=None):
    """LEO mission and {arch: (allocation, closure, assessment)} at the
    nominal allocation."""
    mission = leo_mission()
    out = {}
    for key in ARCHITECTURE_ORDER:
        alloc = nominal_allocation(key, mission)
        res = close(key, mission, alloc, tank_tech_overrides=tank_tech_overrides)
        rep = assess(key, mission, alloc, res)
        out[key] = (alloc, res, rep)
    return mission, out


def nominal_designs_geo():
    """GEO mission and {arch: (allocation, closure, assessment or None)}
    at the nominal allocation."""
    mission = geo_mission()
    out = {}
    for key in ARCHITECTURE_ORDER:
        alloc = nominal_allocation(key, mission)
        res = close(key, mission, alloc)
        rep = assess(key, mission, alloc, res) if res.converged else None
        out[key] = (alloc, res, rep)
    return mission, out


def multimode_all_ep_transfer(mission):
    """Multimode on the GEO mission with the whole transfer electric."""
    alloc = nominal_allocation("multimode", mission)
    tr = mission.get_class("transfer")
    alloc["transfer"] = {"chem": 0.0, "ep": tr.dv_impulsive * tr.ep_loss_ratio}
    res = close("multimode", mission, alloc)
    rep = assess("multimode", mission, alloc, res)
    return res, rep


def _sk_duty_rows(mission, runs):
    """Station-keeping burn time against the available sunlit window."""
    rows = []
    for key in ARCHITECTURE_ORDER:
        _, r, _rep = runs[key]
        if r.diverged or not r.converged:
            print(f"  {key:12s}: closure diverged (no fixed point) — no duty row")
            continue
        for mode in ("chem", "ep"):
            mp = r.per_class_prop["sk"][mode]
            if mp <= 0.0:
                continue
            unit = (
                get_architecture(key).chem_unit
                if mode == "chem"
                else get_architecture(key).ep_unit
            )
            n = r.n_units[mode]
            t_burn = mp / (n * unit.mdot)
            cls_sk = next(c for c in mission.classes if c.key == "sk")
            avail = mission.duration_s * cls_sk.w_duty
            duty = t_burn / avail
            rows.append(
                {
                    "architecture": key,
                    "mode": mode,
                    "sk_prop_kg": round(mp, 3),
                    "t_burn_d": round(t_burn / 86400.0, 2),
                    "available_window_d": round(avail / 86400.0, 1),
                    "duty": round(duty, 4),
                    "frac_of_limit": round(duty / value("sk_duty_max"), 4),
                }
            )
    return rows


def _design_row(mission, key, alloc, r, rep):
    if r.diverged:
        return {
            "architecture": get_architecture(key).display,
            "wet_kg": float("nan"),
            "dry_kg": float("nan"),
            "prop_kg": float("nan"),
            "tanks_kg": float("nan"),
            "eps_inc_kg": float("nan"),
            "n_ep": r.n_units["ep"],
            "n_chem": r.n_units["chem"],
            "f_ep_derived": round(derived_f_ep(mission, alloc), 3),
            "feasible": "DOES NOT CLOSE",
        }
    return {
        "architecture": get_architecture(key).display,
        "wet_kg": round(r.wet, 2),
        "dry_kg": round(r.dry, 2),
        "prop_kg": round(r.prop, 2),
        "tanks_kg": round(r.tanks, 2),
        "eps_inc_kg": round(r.eps_inc, 2),
        "n_ep": r.n_units["ep"],
        "n_chem": r.n_units["chem"],
        "f_ep_derived": round(derived_f_ep(mission, alloc), 3),
        "feasible": "yes" if not rep.excluded else "EXCLUDED",
    }


@task(
    "spacecraft_results",
    "Parameter register; LEO and GEO designs, "
    "tank type, feasibility, duty, transfer check, eta_u sweeps",
)
def run_spacecraft_results(ctx):
    path = export_register_csv()
    print(f"parameter register -> {path}")
    print(
        f"COPV PF (from the Xe tank fraction "
        f"{value('alpha_xe_copv'):.3f}): {copv_pf():.0f} m"
    )

    # LEO mission
    section("LEO mission: 500 km, 150 kg bus, 5 yr (nominal, N2O metal tank)")
    mission, runs = nominal_designs_leo()
    df = pd.DataFrame([_design_row(mission, k, *runs[k]) for k in ARCHITECTURE_ORDER])
    print(df.to_string(index=False))
    out = csv_dir(ctx) / "leo_designs.csv"
    df.to_csv(out, index=False)
    print(f"-> {out}")

    d_mm_se = runs["multimode"][1].wet - runs["single_elec"][1].wet
    d_hy_mm = runs["hybrid"][1].wet - runs["multimode"][1].wet
    print(
        f"\nDifferences: multimode - single_elec = {d_mm_se:+.2f} kg | "
        f"hybrid - multimode = {d_hy_mm:+.2f} kg"
    )

    section("LEO, composite N2O tank")
    _, runs_copv = nominal_designs_leo(tank_tech_overrides={"n2o": "copv"})
    for key in ("multimode", "hybrid", "single_chem"):
        print(
            f"  {key:12s}: {runs_copv[key][1].wet:6.2f} kg wet "
            f"(metal: {runs[key][1].wet:.2f})"
        )

    section("LEO feasibility (t_warn = 1 d)")
    for key in ARCHITECTURE_ORDER:
        rep = runs[key][2]
        if rep.excluded:
            print(f"  {key:12s}: EXCLUDED — " + "; ".join(rep.reasons))
        else:
            print(f"  {key:12s}: ok")

    section("LEO station-keeping duty")
    df_duty = pd.DataFrame(_sk_duty_rows(mission, runs))
    print(df_duty.to_string(index=False))
    out_duty = csv_dir(ctx) / "leo_sk_duty.csv"
    df_duty.to_csv(out_duty, index=False)
    print(f"  -> {out_duty}")

    section("LEO collision avoidance against warning time, single_elec")
    print(
        f"  {'t_warn':>8s} {'dv/event':>9s} {'charge':>8s} "
        f"{'t_burn':>8s} {'limit':>8s}  verdict"
    )
    _ca_ok = {}
    for tw in T_WARN_SET_S:
        m_tw = leo_mission(t_warn_s=tw)
        al = nominal_allocation("single_elec", m_tw)
        r = close("single_elec", m_tw, al)
        rep = assess("single_elec", m_tw, al, r)
        v = next(v for v in rep.verdicts if v.class_key == "ca")
        _ca_ok[tw] = v.feasible
        fac = f"{r.ca_charge['factor']:.3f}" if r.ca_charge["drift_ok"] else "no root"
        print(
            f"  {tw / 3600:6.0f} h {ca_dv_per_event(tw):7.2f} m/s "
            f"{fac:>8s} {v.t_burn_s / 3600:6.2f} h "
            f"{v.t_limit_s / 3600:6.2f} h  "
            f"{'feasible' if v.feasible else 'INFEASIBLE'}"
        )
    _fail_hi = max((t for t, ok in _ca_ok.items() if not ok), default=None)
    _pass_lo = min((t for t, ok in _ca_ok.items() if ok), default=None)
    _brk = (
        "no change on this set"
        if _fail_hi is None or _pass_lo is None
        else f"between {_fail_hi / 3600:.0f} h and {_pass_lo / 3600:.0f} h"
    )
    print(f"  verdict changes {_brk}")

    # GEO mission
    section("GEO mission: transfer, SK and graveyard (300 kg bus, T_max = 6 mo)")
    mission_b, runs_b = nominal_designs_geo()
    df_b = pd.DataFrame(
        [_design_row(mission_b, k, *runs_b[k]) for k in ARCHITECTURE_ORDER]
    )
    print(df_b.to_string(index=False))
    out_b = csv_dir(ctx) / "geo_designs.csv"
    df_b.to_csv(out_b, index=False)
    print(f"-> {out_b}")

    d_mm_hy = runs_b["multimode"][1].wet - runs_b["hybrid"][1].wet
    print(
        f"\nDifferences: multimode - hybrid = {d_mm_hy:+.2f} kg | "
        f"hybrid - single_elec = "
        f"{runs_b['hybrid'][1].wet - runs_b['single_elec'][1].wet:+.2f}"
        f" kg"
    )
    print(
        f"  mass-fraction limit: {value('envelope_frac') * mission_b.bus_mass:.0f} kg"
    )

    r_cf, rep_cf = multimode_all_ep_transfer(mission_b)
    print(
        f"  multimode with an all-electric transfer: n_ECR = "
        f"{r_cf.n_units['ep']}, EPS +{r_cf.eps_inc:.1f} kg vs "
        f"{rep_cf.eps_cap:.0f} kg limit "
        f"(x{r_cf.eps_inc / rep_cf.eps_cap:.1f})"
    )

    section("Closed-form transfer check, GEO (T_transfer_max = 6 months)")
    q = transfer_asymmetry_check()
    print(
        f"  T_req = {q['T_req_N'] * 1e3:.1f} mN on "
        f"{value('geo_m_avg'):.0f} kg average stack"
    )
    for name in ("hall", "ecr"):
        d = q[name]
        print(
            f"  {name:4s}: n = {d['n_units']:3d} units | "
            f"{d['p_dc_w'] / 1e3:5.2f} kW DC | EPS +{d['eps_inc']:5.1f} kg "
            f"| electronics {d['electronics']:5.1f} kg | power "
            f"{'ok' if d['power_ok'] else 'EXCEEDS the limit'}"
        )
    print(
        f"  chem: {q['chem']['m_prop']:.1f} kg transfer propellant on "
        f"the {q['chem']['stack']:.0f} kg stack"
    )

    section("GEO, composite N2O tank")
    runs_b_copv = {}
    for key in ("single_elec", "multimode", "hybrid", "single_chem"):
        alloc = nominal_allocation(key, mission_b)
        r_c = close(key, mission_b, alloc, tank_tech_overrides={"n2o": "copv"})
        runs_b_copv[key] = r_c
        r_m = runs_b[key][1]
        print(
            f"  {key:12s}: metal {r_m.wet:8.2f} kg -> COPV "
            f"{r_c.wet:8.2f} kg (tanks {r_m.tanks:7.2f} -> "
            f"{r_c.tanks:6.2f} kg)"
        )
    gap_m = runs_b["multimode"][1].wet - runs_b["hybrid"][1].wet
    gap_c = runs_b_copv["multimode"].wet - runs_b_copv["hybrid"].wet
    print(
        f"  multimode - hybrid: {gap_m:+.2f} kg (metal) -> {gap_c:+.2f} kg (N2O COPV)"
    )

    section("GEO, all chemical tanks composite")
    m_cross_metal = floor_crossover("propene", "metal")
    print(
        f"  propene load at which a metal tank leaves its minimum mass: "
        f"{m_cross_metal:.2f} kg"
    )
    runs_b_tog = {}
    with override(chem_tank_technology="copv"):
        for key in ("single_elec", "multimode", "hybrid", "single_chem"):
            alloc = nominal_allocation(key, mission_b)
            runs_b_tog[key] = close(key, mission_b, alloc)
    rows_pair = []
    for key in ("single_elec", "multimode", "hybrid", "single_chem"):
        r_m = runs_b[key][1]
        r_t = runs_b_tog[key]
        r_c = runs_b_copv[key]
        rows_pair.append(
            {
                "architecture": get_architecture(key).display,
                "wet_metal_kg": round(r_m.wet, 2),
                "wet_copv_kg": round(r_t.wet, 2),
                "tanks_metal_kg": round(r_m.tanks, 2),
                "tanks_copv_kg": round(r_t.tanks, 2),
                "delta_vs_n2o_only_kg": round(r_t.wet - r_c.wet, 2),
            }
        )
        print(
            f"  {key:12s}: metal {r_m.wet:8.2f} kg -> copv "
            f"{r_t.wet:8.2f} kg (vs N2O-only override "
            f"{r_t.wet - r_c.wet:+7.2f} kg)"
        )
    df_pair = pd.DataFrame(rows_pair)
    out_pair = csv_dir(ctx) / "geo_tank_type.csv"
    df_pair.to_csv(out_pair, index=False)
    print(f"-> {out_pair}")
    gap_t = runs_b_tog["multimode"].wet - runs_b_tog["hybrid"].wet
    print(
        f"  multimode - hybrid: {gap_m:+.2f} kg (metal) -> "
        f"{gap_t:+.2f} kg (all composite; N2O composite only "
        f"{gap_c:+.2f} kg)"
    )

    section("GEO station-keeping duty")
    df_duty_b = pd.DataFrame(_sk_duty_rows(mission_b, runs_b))
    print(df_duty_b.to_string(index=False))
    out_duty_b = csv_dir(ctx) / "geo_sk_duty.csv"
    df_duty_b.to_csv(out_duty_b, index=False)
    print(f"  -> {out_duty_b}")

    section("Idle hardware of the multimode")
    off0, off1 = expected_idle_offsets()
    print(
        f"  f_EP = 0: multimode = single_chem + {off0:.2f} kg "
        f"(idle ECR head and rated RF chain, with structure)"
    )
    print(
        f"  f_EP = 1: multimode = single_elec + {off1:.2f} kg "
        f"(B20 head and empty propene tank, with structure) "
        f"+ Isp and hardware differences"
    )

    # eta_u sweeps
    section("eta_u of N2O, bipropellant multimode")
    for mission_i in (mission, mission_b):
        sw = sweep_eta_u(mission_i)
        mv = sw["min_viable"]
        print(f"\n  {mission_i.display}")
        print(
            f"    references [kg]: multimode f_EP=0 "
            f"{sw['references']['multimode_f0']:.2f} | single_chem "
            f"{sw['references']['single_chem']:.2f} | single_elec "
            f"{sw['references']['single_elec']:.2f} | hybrid "
            f"{sw['references']['hybrid']:.2f}"
        )
        for crit in ("own_chem", "single_chem"):
            x = mv[crit]
            print(
                f"    MINIMUM VIABLE eta_u ({crit:11s}): "
                + (f"{x:.5f}" if x is not None else "none in range")
            )
        n_nc = sum(1 for r in sw["records"] if not r["converged"])
        if n_nc:
            print(f"    {n_nc} grid point(s) did not converge")
        df_eu = pd.DataFrame(
            [
                {
                    k: r[k]
                    for k in (
                        "eta_u",
                        "in_range",
                        "wet",
                        "isp",
                        "n_ep",
                        "converged",
                        "excluded",
                    )
                }
                for r in sw["records"]
            ]
        )
        df_eu["n_ep"] = df_eu["n_ep"].astype("Int64")
        out_eu = csv_dir(ctx) / f"sweep_eta_u_{mission_i.key}.csv"
        df_eu.to_csv(out_eu, index=False)
        print(f"    -> {out_eu}")

    section("eta_u of N2O, monopropellant multimode")
    for mission_i in (mission, mission_b):
        sw = sweep_eta_u(mission_i, arch_key="mono_multimode")
        mv = sw["min_viable"]
        print(f"\n  {mission_i.display}")
        n_div = sum(1 for r in sw["records"] if r["diverged"])
        if mission_i.key == "geo":
            n_tot = len(sw["records"])
            div_eu = [r["eta_u"] for r in sw["records"] if r["diverged"]]
            if n_div == n_tot:
                head = f"    every grid point diverged ({n_div}/{n_tot})"
            elif n_div:
                head = (
                    f"    {n_div}/{n_tot} grid points diverged "
                    f"(eta_u <= {max(div_eu):.2f}); the rest close"
                )
            else:
                head = f"    all {n_tot} grid points close"
            print(head)
            for crit in ("own_chem", "single_chem"):
                x = mv[crit]
                print(
                    f"    MINIMUM VIABLE eta_u ({crit:11s}): "
                    + (
                        f"{x:.5f}"
                        if x is not None
                        else "none in range (reference diverges)"
                    )
                )
        else:
            print(
                f"    references [kg]: own all-chem (160 s, f_EP=0) "
                f"{sw['references']['own_f0']:.2f} | single_chem "
                f"{sw['references']['single_chem']:.2f}"
            )
            for crit in ("own_chem", "single_chem"):
                x = mv[crit]
                print(
                    f"    MINIMUM VIABLE eta_u ({crit:11s}): "
                    + (f"{x:.5f}" if x is not None else "none in range")
                )
        df_eu = pd.DataFrame(
            [
                {
                    k: r[k]
                    for k in (
                        "eta_u",
                        "in_range",
                        "wet",
                        "isp",
                        "n_ep",
                        "converged",
                        "diverged",
                        "excluded",
                    )
                }
                for r in sw["records"]
            ]
        )
        df_eu["n_ep"] = df_eu["n_ep"].astype("Int64")
        out_eu = csv_dir(ctx) / (f"sweep_eta_u_mono_multimode_{mission_i.key}.csv")
        df_eu.to_csv(out_eu, index=False)
        print(f"    -> {out_eu}")


@task(
    "spacecraft_sweeps",
    "Warning-time, transfer-window, one-at-a-time "
    "and planning-efficiency sweeps; minimum viable eta_u",
)
def run_spacecraft_sweeps(ctx):
    section("Warning-time sweep, LEO")
    sw = sweep_t_warn(ARCHITECTURE_ORDER)
    rows = []
    for key, recs in sw["records"].items():
        for r in recs:
            rows.append(
                {
                    "architecture": key,
                    "t_warn_h": r["t_warn_s"] / 3600.0,
                    "wet_kg": r["wet"],
                    "feasible": r["feasible"],
                    "ca_prop_kg": r["ca_prop"],
                    "ca_factor": r["ca_factor"],
                    "drift_ok": r["drift_ok"],
                    "n_ep": r["n_ep"],
                    "reasons": "; ".join(r["reasons"]),
                }
            )
    df = pd.DataFrame(rows)
    out = csv_dir(ctx) / "sweep_t_warn_leo.csv"
    df.to_csv(out, index=False)
    print(
        f"  {len(sw['grid'])} grid points x {len(sw['records'])} architectures -> {out}"
    )
    print(f"  single_elec CA limit (bisection): {sw['flip_s'] / 3600:.2f} h")
    print("  Lightest feasible design:")
    optimal_at = {
        tw: o for o in sw["optimal"] for tw in T_WARN_SET_S if o["t_warn_s"] == tw
    }
    for tw in T_WARN_SET_S:
        o = optimal_at[tw]
        print(f"    t_warn = {tw / 3600:5.0f} h: {o['arch']} ({o['wet']:.2f} kg)")

    section("Single-electric warning-time limit, sensitivity")
    fb = ca_flip_band()
    rows = [
        {
            "param": "nominal",
            "value": float("nan"),
            "flip_h": fb["nominal"] / 3600.0,
            "hand_sqrt_h": float("nan"),
        }
    ]
    for e in fb["band"] + fb["w_sensitivity"]:
        rows.append(
            {
                "param": e["param"],
                "value": e["value"],
                "flip_h": e["flip_s"] / 3600.0,
                "sqrt_scaling_h": e["sqrt_scaling_s"] / 3600.0,
            }
        )
    df_fb = pd.DataFrame(rows)
    out_fb = csv_dir(ctx) / "ca_flip_band_leo.csv"
    df_fb.to_csv(out_fb, index=False)
    print(f"  -> {out_fb}")
    print(
        f"  nominal {fb['nominal'] / 3600:.2f} h; range over the "
        f"dv_ca_budget_1d bracket "
        f"{REGISTER['dv_ca_budget_1d'].range}: "
        f"{fb['band_lo_s'] / 3600:.1f} - {fb['band_hi_s'] / 3600:.1f} h"
    )
    for e in fb["w_sensitivity"]:
        print(
            f"  w_duty_leo = {e['value']}: {e['flip_s'] / 3600:.2f} h "
            f"(sqrt scaling {e['sqrt_scaling_s'] / 3600:.2f} h)"
        )

    section("Transfer-window sweep, GEO")
    sw_t = sweep_t_transfer_max()
    rows = []
    for key, recs in sw_t["records"].items():
        for r in recs:
            rows.append(
                {
                    "architecture": key,
                    "mode": "nominal",
                    "t_max_months": r["t_yr"] * 12,
                    "wet_kg": r["wet"],
                    "s_ep": float("nan"),
                    "n_ep": r["n_ep"],
                    "n_max": float("nan"),
                    "p_dc_w": r["p_dc_w"],
                    "eps_inc_kg": r["eps_inc"],
                    "electronics_kg": r["electronics"],
                    "tr_prop_ep_kg": float("nan"),
                    "tr_prop_chem_kg": float("nan"),
                    "t_commission_d": r["t_comm_s"] / 86400.0,
                    "feasible": r["feasible"],
                    "reasons": "; ".join(r["reasons"]),
                }
            )
    for key, recs in sw_t["ep_transfer"].items():
        for r in recs:
            rows.append(
                {
                    "architecture": key,
                    "mode": "ep_transfer",
                    "t_max_months": r["t_yr"] * 12,
                    "wet_kg": r["wet"],
                    "s_ep": r["s_ep"],
                    "n_ep": r["n_ep"],
                    "n_max": r["n_max"],
                    "p_dc_w": r["p_dc_w"],
                    "eps_inc_kg": r["eps_inc"],
                    "electronics_kg": float("nan"),
                    "tr_prop_ep_kg": r["tr_prop_ep"],
                    "tr_prop_chem_kg": r["tr_prop_chem"],
                    "t_commission_d": r["t_comm_s"] / 86400.0,
                    "feasible": r["feasible"],
                    "reasons": "; ".join(r["reasons"]),
                }
            )
    df_t = pd.DataFrame(rows)
    out_t = csv_dir(ctx) / "sweep_t_transfer_max_geo.csv"
    df_t.to_csv(out_t, index=False)
    print(f"  -> {out_t}")
    print(
        "  Time-optimal single_elec and the maximum-electric split "
        "(s* = largest electric share within the power-system limit):"
    )
    print(
        f"  {'T_max':>8s} | {'se nominal n / verdict':>22s} | "
        f"{'mm split s* / wet kg':>21s} | "
        f"{'hy split s* / wet kg / verdict':>30s}"
    )
    by_mm = {r["t_yr"]: r for r in sw_t["ep_transfer"]["multimode"]}
    by_hy = {r["t_yr"]: r for r in sw_t["ep_transfer"]["hybrid"]}
    for t, r_se in zip(sw_t["t_set_yr"], sw_t["records"]["single_elec"]):
        se = (
            "DIVERGES"
            if r_se["diverged"]
            else f"{r_se['n_ep']:d} / "
            + ("feasible" if r_se["feasible"] else "excluded")
        )
        r_m, r_h = by_mm[t], by_hy[t]
        mm = f"{r_m['s_ep']:.3f} / {r_m['wet']:7.1f}"
        hy = f"{r_h['s_ep']:.3f} / {r_h['wet']:7.1f} / " + (
            "feasible" if r_h["feasible"] else "excluded"
        )
        print(f"  {t * 12:5.1f} mo | {se:>22s} | {mm:>21s} | {hy:>30s}")

    section("One-at-a-time sensitivity (tornado)")
    tor = tornado()
    base = tor["baseline"]
    print(
        f"  baseline: leo d(mm - se) = "
        f"{base['leo']['d_mm_se']:+.3f} kg | geo d(mm - hy) = "
        f"{base['geo']['d_mm_hy']:+.1f} kg"
    )
    stats = {}
    for mkey, metric in (("leo", "d_mm_se"), ("geo", "d_mm_hy")):
        rows, flips, errors, vals = [], [], [], []
        for bar in tor["bars"]:
            for e in ("lo", "hi"):
                d = bar["edges"][e]
                row = {
                    "bar": bar["name"],
                    "group": bar["group"],
                    "flag": bar["flag"],
                    "edge": e,
                    "value": d["value_label"],
                    "error": d["error"] or "",
                }
                if d["error"]:
                    errors.append(f"{bar['name']} {e} ({d['value_label']})")
                else:
                    row.update(
                        {
                            k: d[mkey][k]
                            for k in (
                                "wet_mm",
                                "wet_se",
                                "wet_hy",
                                "wet_mmm",
                                "d_mm_se",
                                "d_mm_hy",
                                "d_mmm_mm",
                            )
                        }
                    )
                    v = d[mkey][metric]
                    if v == v:
                        vals.append((v, bar["name"], e))
                        if v * base[mkey][metric] < 0:
                            flips.append(f"{bar['name']} {e} ({d['value_label']})")
                rows.append(row)
        vals.sort()
        stats[mkey] = {
            "baseline": base[mkey][metric],
            "n_eval": len(vals),
            "lo": vals[0] if vals else None,
            "hi": vals[-1] if vals else None,
            "flips": flips,
        }
        df_tor = pd.DataFrame(rows)
        out_tor = csv_dir(ctx) / f"tornado_{mkey}.csv"
        df_tor.to_csv(out_tor, index=False)
        print(f"\n  {mkey} -> {out_tor}")
        if flips:
            print(f"  sign changes of {metric} ({len(flips)}): " + "; ".join(flips))
        else:
            print(f"  no sign change of {metric}")
        if mkey == "leo" and errors:
            print("  not evaluable: " + "; ".join(sorted(set(errors))))
    for mkey, metric in (("leo", "d(mm - se)"), ("geo", "d(mm - hy)")):
        s = stats[mkey]
        print(
            f"\n  {mkey} {metric}: baseline {s['baseline']:+.3f} kg, "
            f"range {s['lo'][0]:+.3f} ({s['lo'][1]} {s['lo'][2]}) .. "
            f"{s['hi'][0]:+.3f} ({s['hi'][1]} {s['hi'][2]}) kg, "
            f"{len(s['flips'])} of {s['n_eval']} cases change sign"
        )
    base_a = base["leo"]
    hb = next((b for b in tor["bars"] if b["name"] == "p_heater_mono"), None)
    if hb is not None:
        widths = []
        for e in ("lo", "hi"):
            d = hb["edges"][e]
            if not d["error"] and d["leo"]["wet_mmm"] == d["leo"]["wet_mmm"]:
                widths.append(abs(d["leo"]["wet_mmm"] - base_a["wet_mmm"]))
        w_max = max(widths) if widths else float("nan")
        print(
            f"\n  mono_multimode, LEO: baseline wet "
            f"{base_a['wet_mmm']:.4f} kg, d(mono_mm - mm) = "
            f"{base_a['d_mmm_mm']:+.4f} kg; p_heater_mono 10-50 W "
            f"moves the wet mass by at most {w_max:.2e} kg"
        )

    section("Minimum viable eta_u")
    sw_eu = {}
    for mission_i in (leo_mission(), geo_mission()):
        sw_eu[mission_i.key] = (mission_i, sweep_eta_u(mission_i))
    for mkey, (mission_i, sw) in sw_eu.items():
        mv = sw["min_viable"]
        print(
            f"  {mkey} minimum viable eta_u: own_chem "
            + (f"{mv['own_chem']:.4f}" if mv["own_chem"] else "none")
            + ", single_chem "
            + (f"{mv['single_chem']:.4f}" if mv["single_chem"] else "none")
        )
    print(f"  measured H2O utilisation range: {MEASURED_H2O_ETA_U}")
    _m_b = sw_eu["geo"][0]
    _alloc_dm, _diag_dm = ep_transfer_allocation("multimode", _m_b)
    _r_dm = close("multimode", _m_b, _alloc_dm)
    print(
        f"  GEO multimode with the maximum-electric transfer at parity: "
        f"{_r_dm.wet:.2f} kg, s* = {_diag_dm['s_ep']:.4f}"
    )

    section("Planning efficiency and array power loss, GEO")
    dv_chem = value("dv_transfer_chem")
    lo_b, hi_b = dv_transfer_ep_range()
    print(
        f"  eta_planning = dv_transfer_chem / dv_transfer_ep = "
        f"{dv_chem:g} / {dv_transfer_ep():g} m/s = "
        f"{dv_chem / dv_transfer_ep():.4f}; "
        f"dv_transfer_ep bracket ({lo_b:g}, {hi_b:g}) m/s; "
        f"register bracket {REGISTER['eta_planning_transfer'].range}"
    )

    sep = sweep_eta_planning(
        eta_set=tuple(
            sorted(
                set(ETA_PLANNING_SET) | {round(value("eta_planning_transfer"), 6)},
                reverse=True,
            )
        )
    )
    combos = sep["combos"]
    m_ep = sep["mission"]
    cap_env = value("envelope_frac") * m_ep.bus_mass

    print(
        f"\n  Designs within the mass-fraction limit ({cap_env:g} kg on "
        f"a {m_ep.bus_mass:g} kg bus):"
    )
    for strategy in sep["strategies"]:
        print(f"    {strategy} -- {STRATEGY_LABEL[strategy]}")
        for eta, d in combos:
            surv = sep["survivors"][strategy][(eta, d)]
            base = sep["survivors"][strategy][combos[0]]
            print(
                f"      eta {eta:.2f}  d {d:.2f}  -> "
                f"{list(surv) if surv else 'EMPTY'}"
                + ("" if surv == base else "   <- CHANGED")
            )
        fc = sep["first_change"].get(strategy)
        print(
            f"      {
                'FIRST CHANGE at eta %.2f, d %.2f' % (fc[0], fc[1])
                if fc
                else 'no combination changes the set'
            }"
        )

    print("\n  single_elec margin to the mass-fraction limit:")
    print(
        f"    {'eta':>5s} {'d':>5s} {'dv_ep':>7s} {'wet[kg]':>9s} "
        f"{'env frac':>9s} {'margin[kg]':>11s} {'verdict':>9s}"
    )
    for eta, d in combos:
        r = next(
            x
            for x in sep["grid"][(eta, d)]["time-optimal"]
            if x["arch"] == "single_elec"
        )
        print(
            f"    {eta:5.2f} {d:5.2f} {r['dv_ep']:7.1f} {r['wet']:9.3f} "
            f"{r['env_frac']:9.4f} {cap_env - r['wet']:11.3f} "
            f"{('PASS' if r['envelope_ok'] else 'EXCLUDED'):>9s}"
        )

    print("\n  Lowest planning efficiency within the limit (bisection, tol 1e-4):")
    for strategy in sep["strategies"]:
        for key in ("single_elec", "hybrid"):
            for d in sep["degr_set"]:
                e = eta_planning_crossover(key, strategy, d)
                if e is None:
                    continue
                print(
                    f"    {strategy} {key:<12s} d = {d:.2f}: down to "
                    f"eta = {e:.4f}  (dv_transfer_ep "
                    f"{dv_chem / e:.0f} m/s)"
                )

    rows = []
    for strategy in sep["strategies"]:
        for eta, d in combos:
            for r in sep["grid"][(eta, d)][strategy]:
                rows.append(
                    {
                        "strategy": strategy,
                        "eta_planning": eta,
                        "arr_degr_transfer": d,
                        "architecture": r["arch"],
                        "dv_transfer_ep_ms": r["dv_ep"],
                        "loss_ratio": r["loss_ratio"],
                        "s_ep": r["s_ep"],
                        "arr_degr_effective": r["d_eff"],
                        "p_dc_peak_w": r["p_dc_w"],
                        "n_ep": r["n_ep"],
                        "eps_inc_kg": r["eps_inc"],
                        "wet_kg": r["wet"],
                        "envelope_frac": r["env_frac"],
                        "envelope_cap_kg": r["envelope_cap"],
                        "power_ok": r["power_ok"],
                        "envelope_ok": r["envelope_ok"],
                        "converged": r["converged"],
                        "diverged": r["diverged"],
                    }
                )
    out_ep = csv_dir(ctx) / "sweep_eta_planning_geo.csv"
    pd.DataFrame(rows).to_csv(out_ep, index=False)
    print(f"\n  -> {out_ep}")
