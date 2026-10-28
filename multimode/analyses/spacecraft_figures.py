"""
Spacecraft-level thesis figures (chapter 6) and the LEO response-variant
sweep. Figures are written as <label>.pdf to PATHS['figures'], or to the
directory in the environment variable THESIS_FIG_DIR; THESIS6_ONLY (comma
separated: wet, ca, eta, fep, twarn, response, tornado) limits the set.
"""

import csv
import json
import os
import sys

import numpy as np

from ..config import PARAMS, PATHS
from ..plotting import thesis_figures as tf
from ..plotting.thesis_print import save_pdf
from .common import section
from .registry import task
from .thruster_figures import figure_dir

# eta_u grid to 1.0: 0.05 to 0.99 in steps of 0.02, and 1.00
ETA_U_GRID_TO_1 = tuple(i / 100.0 for i in range(5, 100, 2)) + (1.0,)


@task(
    "response_variant",
    "LEO response-variant sweep "
    "(scripts/response_variant_sweep.py), input to fig_res_response",
)
def run_response_variant(ctx):
    sys.path.insert(0, str(PATHS["project_root"] / "scripts"))
    import response_variant_sweep

    response_variant_sweep.main()


# Inputs left out of the LEO tornado: they do not act on the multimode and
# single-electric LEO designs, or (the LEO bus mass) have no power margin
# at their bracket ends
_TORNADO_A_DROP = (
    "bus_mass_geo",
    "power_margin_300",
    "dv_transfer_chem",
    "eta_planning_transfer",
    "w_duty_transfer",  # GEO mission only
    "m_head_mono",
    "isp_mono",
    "p_heater_mono",  # monopropellant unit
    "eps_cap_frac",
    "envelope_frac",  # verdicts only
    "bus_mass_leo",  # not evaluable
)


def _reasons(records):
    out = {}
    for r in records:
        for s in r.get("reasons", ()):
            out[s] = out.get(s, 0) + 1
    return out


def _ca_feasible(key, t_warn_h):
    """Collision-avoidance verdict of `key` flown all-electric (f_EP = 1)
    at this warning time [h]."""
    from ..spacecraft.allocation import sweep_allocation
    from ..spacecraft.closure import close
    from ..spacecraft.feasibility import assess
    from ..spacecraft.mission import leo_mission

    m = leo_mission(t_warn_h * 3600.0)
    a = sweep_allocation(key, m, 1.0)
    r = close(key, m, a)
    if not r.converged or r.diverged:
        return False
    rep = assess(key, m, a, r)
    v = [x for x in rep.verdicts if x.class_key == "ca"][0]
    return bool(v.feasible)


def _ca_min_twarn(key, lo=2.0, hi=72.0, tol=0.01):
    """Shortest warning time [h] at which `key` passes the check."""
    if _ca_feasible(key, lo):
        return float("-inf")
    if not _ca_feasible(key, hi):
        return float("nan")
    a, b = lo, hi
    while b - a > tol:
        mid = 0.5 * (a + b)
        if _ca_feasible(key, mid):
            b = mid
        else:
            a = mid
    return b


def _ca_feasible_event(key, dv_event, t_warn_h):
    """Collision-avoidance verdict at a given per-event delta-v at the
    warning time: the budget at the reference warning time is set so that
    ca_dv_per_event(t_warn) = dv_event."""
    from ..spacecraft.parameters import value

    b1d = dv_event * t_warn_h * 3600.0 / value("t_warn_ref_ca")
    return _with_budget(b1d, _ca_feasible, key, t_warn_h)


def _ca_min_twarn_event(key, dv_event, lo=2.0, hi=72.0, tol=0.01):
    """Shortest warning time [h] at which the electric burn of this
    per-event delta-v fits; -inf if it fits at lo, nan if not at hi."""
    if _ca_feasible_event(key, dv_event, lo):
        return float("-inf")
    if not _ca_feasible_event(key, dv_event, hi):
        return float("nan")
    a, b = lo, hi
    while b - a > tol:
        mid = 0.5 * (a + b)
        if _ca_feasible_event(key, dv_event, mid):
            b = mid
        else:
            a = mid
    return b


def _with_budget(dv, fn, *args):
    from ..spacecraft.architectures import invalidate_cache
    from ..spacecraft.parameters import override

    try:
        with override(dv_ca_budget_1d=dv):
            invalidate_cache()
            return fn(*args)
    finally:
        invalidate_cache()


@task(
    "spacecraft_figures",
    "Chapter 6 figures: wet-mass bars, "
    "collision-avoidance map, eta_u thresholds, allocation sweep, "
    "warning-time pivot, response variant, tornado",
)
def run_spacecraft_figures(ctx):
    from ..spacecraft.allocation import nominal_allocation
    from ..spacecraft.architectures import ARCHITECTURE_ORDER
    from ..spacecraft.closure import close
    from ..spacecraft.feasibility import assess
    from ..spacecraft.mission import leo_mission, geo_mission
    from ..spacecraft.parameters import value
    from ..spacecraft.sweeps import (
        ETA_U_RANGE,
        ca_flip_band,
        sweep_eta_u,
        sweep_f_ep,
        sweep_t_warn,
        tornado,
    )

    out = figure_dir()
    only = os.environ.get("THESIS6_ONLY", "")

    def want(k):
        return not only or k in only.split(",")

    if want("wet"):
        section("Wet-mass bars, nominal designs, LEO and GEO")
        for mk, mission_fn, name in (
            ("leo", leo_mission, "fig_res_leo_wet"),
            ("geo", geo_mission, "fig_res_geo_wet"),
        ):
            m = mission_fn()
            env = value("envelope_frac") * m.bus_mass
            rows = []
            for k in ARCHITECTURE_ORDER:
                a = nominal_allocation(k, m)
                r = close(k, m, a)
                ok = r.converged and not r.diverged
                rep = assess(k, m, a, r) if ok else None
                rows.append(
                    dict(key=k, dry=r.dry, prop=r.prop, wet=r.wet, converged=ok)
                )
                print(
                    f"  {mk} {k:15s} dry {r.dry:9.3f} + prop {r.prop:9.3f}"
                    f" = wet {r.wet:9.3f} kg  converged {ok}  "
                    + (
                        ""
                        if rep is None
                        else f"envelope_ok {rep.envelope_ok}  excluded "
                        f"{rep.excluded}  {'; '.join(rep.reasons)}"
                    )
                )
            print(
                f"  {mk} mass-fraction limit {env:g} kg "
                f"(envelope_frac {value('envelope_frac'):g} x bus "
                f"{m.bus_mass:g} kg)"
            )
            wets = sorted(r["wet"] for r in rows if r["converged"])
            cut = 1.25 * wets[-2] if wets[-1] > 2.0 * wets[-2] else None
            save_pdf(tf.fig_res_wet_bars(rows, truncate_above=cut), out, f"{name}.pdf")

    if want("ca"):
        section(
            "LEO collision-avoidance map (multimode at f_EP = 1; "
            "per-event delta-v at the warning time)"
        )
        t_ref_h = value("t_warn_ref_ca") / 3600.0
        dv_ev = sorted(set(np.round(np.geomspace(0.01, 1.0, 25), 6)))
        boundary = []
        for d in dv_ev:
            t = _ca_min_twarn_event("multimode", d)
            boundary.append((d, t))
            print(
                f"  dv_event {d:7.4f} m/s: shortest warning for the ECR burn {t:8.3f} h"
            )
        # Boundary point at the left edge (t_warn = 2 h), by bisection
        t_edge = 2.0
        for (d0, t0), (d1, t1) in zip(boundary, boundary[1:]):
            if np.isneginf(t0) and np.isfinite(t1):
                a, b = d0, d1
                while b / a - 1.0 > 1e-5:
                    mid = float(np.sqrt(a * b))
                    if _ca_feasible_event("multimode", mid, t_edge):
                        a = mid
                    else:
                        b = mid
                boundary.insert(boundary.index((d1, t1)), (b, t_edge))
                print(
                    f"  left edge: dv_event {b:.5f} m/s needs {t_edge:g} h "
                    f"(fits at {t_edge:g} h up to {a:.5f} m/s)"
                )
                break
        grid = []
        for d in (0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0):
            row = []
            for th in (2, 4, 6, 12, 24, 48, 72):
                ok = _ca_feasible_event("multimode", d, th)
                grid.append((th, d, ok))
                row.append("ECR" if ok else "chem")
            print(
                f"  spot grid dv_event {d:4.2f}: "
                + " ".join(
                    f"{th}h:{v}" for th, v in zip((2, 4, 6, 12, 24, 48, 72), row)
                )
            )
        b1d = value("dv_ca_budget_1d")
        t_nom = value("t_warn") / 3600.0
        # Minimum warning time at the nominal budget
        t_cross = _with_budget(b1d, _ca_min_twarn, "multimode")
        nominal = (t_nom, b1d * t_ref_h / t_nom)
        print(
            f"  nominal point: {nominal[1]:g} m/s per event at "
            f"{nominal[0]:g} h; budget line {b1d:g} m/s x "
            f"({t_ref_h:g} h / t_warn); crosses the boundary at "
            f"{t_cross:.2f} h ({b1d * t_ref_h / t_cross:.4f} m/s)"
        )
        save_pdf(
            tf.fig_res_ca_map(
                boundary,
                nominal=nominal,
                t_ref_h=t_ref_h,
                crossing_h=t_cross,
                grid_pts=grid,
            ),
            out,
            "fig_res_ca_map.pdf",
        )

    if want("eta"):
        section("Minimum viable eta_u, grid to 1.0")
        sw = {}
        for mk, mission_fn in (("leo", leo_mission), ("geo", geo_mission)):
            s = sweep_eta_u(mission_fn(), eta_u_grid=ETA_U_GRID_TO_1)
            sw[mk] = s
            refs = s["references"]
            print(
                f"  {mk}: min viable own_chem {s['min_viable']['own_chem']}"
                f", single_chem {s['min_viable']['single_chem']}; refs "
                + ", ".join(f"{k} {v:.3f}" for k, v in refs.items())
            )
            ext = [(r["eta_u"], r["wet"]) for r in s["records"] if r["eta_u"] >= 0.61]
            print("    extension: " + ", ".join(f"{e:.2f}:{w:.2f}" for e, w in ext))
        save_pdf(
            tf.fig_res_min_eta_u(
                sw["leo"],
                sw["geo"],
                eta_range=ETA_U_RANGE,
                eta_measured=PARAMS["eta_u_xe_ref"],
            ),
            out,
            "fig_res_min_eta_u.pdf",
        )

    if want("fep"):
        section("Allocation sweep, GEO mission")
        swf = sweep_f_ep(ARCHITECTURE_ORDER, geo_mission())
        allrecs = [r for k in swf for r in swf[k]["records"]]
        kinds = {}
        for k, v in _reasons([r for r in allrecs if r["excluded"]]).items():
            kinds[k.split(":")[0]] = kinds.get(k.split(":")[0], 0) + v
        print(
            "  exclusion reasons by kind (count): "
            + "; ".join(f"{k} [{v}]" for k, v in kinds.items())
            + f"; excluded points {sum(r['excluded'] for r in allrecs)}"
            f" of {len(allrecs)}"
        )
        print(
            "  derived f_EP: "
            + ", ".join(f"{k} {swf[k].get('f_ep_derived')}" for k in swf)
        )
        both = sum(
            1
            for r in allrecs
            if r["excluded"]
            and any(x.startswith("power") for x in r["reasons"])
            and any(x.startswith("envelope") for x in r["reasons"])
        )
        print(
            f"  points over the power-system limit and the mass-fraction limit: {both}"
        )
        print(f"  points not excluded: {sum(not r['excluded'] for r in allrecs)}")
        save_pdf(
            tf.fig_res_fep_geo(swf, power_label="over power-system limit"),
            out,
            "fig_res_fep_geo.pdf",
        )

    if want("twarn"):
        section("LEO warning-time pivot")
        swt = sweep_t_warn(ARCHITECTURE_ORDER)
        fb = ca_flip_band()
        bad = [
            r for k in swt["records"] for r in swt["records"][k] if not r["feasible"]
        ]
        print(
            "  infeasibility reasons (count): "
            + "; ".join(f"{k} [{v}]" for k, v in _reasons(bad).items())
        )
        print(
            f"  flip {swt['flip_s'] / 3600:.2f} h, band "
            f"{fb['band_lo_s'] / 3600:.2f}-{fb['band_hi_s'] / 3600:.2f} h"
        )
        save_pdf(
            tf.fig_res_twarn_pivot(
                swt,
                fb,
                nominal_h=value("t_warn") / 3600.0,
                open_label="burn exceeds warning time",
            ),
            out,
            "fig_res_twarn_pivot.pdf",
        )

    if want("response"):
        section("LEO response variant")
        rdir = PATHS["output"] / "response_variant"
        csv_path = rdir / "response_variant.csv"
        if not csv_path.exists():
            raise FileNotFoundError(
                f"{csv_path} not found: run `python run.py response_variant` first"
            )
        with open(csv_path, encoding="utf-8") as fh:
            # an empty wet_kg is an excluded design without a mass
            rrows = [
                dict(
                    key=r["architecture"],
                    total=float(r["total_response_dv_mps"]),
                    wet=(float(r["wet_kg"]) if r["wet_kg"] else float("nan")),
                    feasible=(r["feasible"] == "True" and r["wet_kg"] != ""),
                )
                for r in csv.DictReader(fh)
            ]
        summary = json.loads(
            (rdir / "response_variant.json").read_text(encoding="utf-8")
        )
        se_lim = summary["single_elec_first_excluded_total"]
        changes = [
            (0.5 * (c["total_lo"] + c["total_hi"]), c["before"], c["after"])
            for c in summary["lightest_changes"]
        ]
        print(f"  single-electric first excluded at {se_lim!r} m/s")
        for x, b, a in changes:
            print(f"  lightest feasible: {b} -> {a} at {x!r} m/s")
        save_pdf(
            tf.fig_res_response(rrows, se_limit_total=se_lim),
            out,
            "fig_res_response.pdf",
        )

    if want("tornado"):
        section("Tornado, LEO mission")
        tor = tornado()
        labels = tf.TORNADO_ROW_LABELS
        names = [b["name"] for b in tor["bars"]]
        print(
            f"  {len(names)} bars; baseline d(mm - se) "
            f"{tor['baseline']['leo']['d_mm_se']:+.4f} kg"
        )
        drop = _TORNADO_A_DROP
        base, rows = tf.tornado_rows(tor, "leo", drop=drop)
        n = len(rows)
        print(
            f"  dropped: {', '.join(drop)}; {n} rows drawn; "
            f"{sum(rr['flips'] for rr in rows)} cross zero"
        )
        print("  row  span [kg]  flips  label (1 = top)")
        for k, rr in enumerate(reversed(rows), start=1):
            print(
                f"  {k:3d}  {rr['span']:9.4f}  {str(rr['flips']):5s}  "
                f"{labels[rr['name']]}"
            )
        small = [
            (k, labels[rr["name"]])
            for k, rr in enumerate(reversed(rows), start=1)
            if rr["span"] < 0.01 and rr["devs"]
        ]
        print(
            "  rows with an effect below 0.01 kg (evaluated): "
            + "; ".join(f"{k} {l}" for k, l in small)
        )
        save_pdf(
            tf.fig_res_tornado(tor, "leo", row_labels=labels, drop=drop),
            out,
            "fig_res_tornado_a.pdf",
        )
