"""
Thruster-level thesis figures: chapter 4 results and the chapter 3 check
against Sheppard and Little (2022). Figures are written as <label>.pdf to
PATHS['figures'], or to the directory in the environment variable
THESIS_FIG_DIR.
"""

import json
import os
from pathlib import Path

import numpy as np

from ..config import PARAMS, PATHS
from ..thruster.pipeline import run_pipeline, baseline_at_anchor
from ..thruster.ion_population import analyse_ion_populations, compute_neg_ion_ratios
from ..thruster.propellants import N2O
from ..thruster import sl2022_validation as sl
from ..plotting import thesis_figures as tf
from ..plotting.thesis_figures import ELECTROPOSITIVE_RATIO_MAX
from ..plotting.thesis_print import save_pdf, REPRESENTATIVE_TE
from .common import csv_dir, section
from .registry import task


def figure_dir():
    env = os.environ.get("THESIS_FIG_DIR")
    return Path(env) if env else PATHS["figures"]


def _te_grid():
    return np.arange(
        PARAMS["Te_min"], PARAMS["Te_max"] + PARAMS["Te_step"] / 2, PARAMS["Te_step"]
    )


def electropositive_boundary():
    """T_e above which the N2O negative-ion ratio stays below
    ELECTROPOSITIVE_RATIO_MAX (last downward crossing, interpolated).
    Returns (T_e, ion analysis, total ratio)."""
    Te = _te_grid()
    an = analyse_ion_populations(
        N2O.ion_channels, N2O.ion_energies, Te, neg_ion_channels=N2O.neg_ion_channels
    )
    _ind, total = compute_neg_ion_ratios(an["kappa_neg"], Te)
    xs = tf.crossings(Te, total - ELECTROPOSITIVE_RATIO_MAX)
    return (xs[-1] if xs else None), an, total


@task(
    "thruster_figures",
    "Chapter 4 figures: plasma composition, "
    "ionisation cost, thrust efficiency, specific impulse, power per thrust",
)
def run_thruster_figures(ctx):
    out = figure_dir()
    te_ref = PARAMS["Te_ref"]
    te_b, an, total = electropositive_boundary()
    section("Electron-temperature markers")
    print(
        f"  electropositive boundary (N2O ratio = "
        f"{ELECTROPOSITIVE_RATIO_MAX:g}, last crossing): {te_b:.3f} eV"
    )
    print(
        f"  reference T_e: {te_ref:g} eV; representative ranges: "
        + ", ".join(f"{n} {a:g}-{b:g} eV" for n, a, b in REPRESENTATIVE_TE)
    )
    mk = dict(te_boundary=te_b, te_ref=te_ref)

    save_pdf(tf.fig_ep_res_n2o_ions(an, **mk), out, "fig_ep_res_n2o_ions.pdf")
    only = next(iter(an["kappa_neg"]), "total")
    save_pdf(
        tf.fig_ep_res_neg_ion(
            an, total, criterion=ELECTROPOSITIVE_RATIO_MAX, species_label=only, **mk
        ),
        out,
        "fig_ep_res_neg_ion.pdf",
    )

    r = run_pipeline(
        PARAMS, propellants=("Xe", "Ar", "N2O", "H2O"), anchor_propellant="Xe"
    )
    save_pdf(tf.fig_ep_res_ec(r, **mk), out, "fig_ep_res_ec.pdf")
    save_pdf(tf.fig_ep_res_ec_linear(r, **mk), out, "fig_ep_res_ec_linear.pdf")
    save_pdf(tf.fig_ep_res_h1_categories(r, **mk), out, "fig_ep_res_h1_categories.pdf")

    Te, N, X, inc = tf.h1_groups(r)
    section("Loss groups: crossings drawn")
    for name, y in (
        ("diss = internal (increments)", inc["dissociation"] - inc["internal"]),
        ("ionisation = internal (increments)", inc["ionisation"] - inc["internal"]),
        ("internal increment = 0", inc["internal"]),
    ):
        print(f"  {name}: " + ", ".join(f"{x:.2f} eV" for x in tf.crossings(Te, y)))

    rm_set = (1, 2, 4, 5, 6, 8, 10)
    save_pdf(tf.fig_ep_res_etaT_ratio(r, **mk), out, "fig_ep_res_etaT_ratio.pdf")
    section("eta_T ratio to Xe at the sensitivity set, T_e_ref")
    i = int(np.argmin(np.abs(r["Te"] - te_ref)))
    for num in ("N2O", "Ar"):
        band = tf.eta_T_ratio_band(r, num, "Xe", rm_set)
        print(
            f"  {num}/Xe: " + ", ".join(f"R_m {R:g} {band[R][i]:.4f}" for R in rm_set)
        )
        dense = tf.eta_T_ratio_band(r, num, "Xe", np.linspace(1.0, 10.0, 37))
        arr = np.vstack(list(dense.values()))
        width = np.nanmax(arr, 0) - np.nanmin(arr, 0)
        for lo_t, hi_t in ((1.0, 30.0), (5.0, 30.0)):
            m = (r["Te"] >= lo_t) & (r["Te"] <= hi_t)
            j = int(np.nanargmax(np.where(m, width, -1.0)))
            print(
                f"    band width (R_m 1-10), {lo_t:g}-{hi_t:g} eV: "
                f"max {width[j]:.4f} at {r['Te'][j]:.1f} eV"
            )
    save_pdf(tf.fig_ep_res_etaT_abs(r, **mk), out, "fig_ep_res_etaT_abs.pdf")
    section("eta_T at R_m 5 (absolute), at T_e_ref")
    case = r["params"].get("n2o_case", "B")
    for p in ("Xe", "Ar", "N2O", "H2O"):
        e = tf._eta_T_at_Rm(r[p], r["Te"], 5.0, case)
        print(f"  {p}: {e[i]:.4f}")
    save_pdf(tf.fig_ep_res_isp(r, **mk), out, "fig_ep_res_isp.pdf")
    save_pdf(tf.fig_ep_res_decomposition(r), out, "fig_ep_res_decomposition.pdf")
    tp = tf._thrust_to_power_vs_xe(r)
    print(
        f"  decomposition, power per thrust vs Xe: Xe 1, "
        f"Ar {1 / tp[1]:.4f}, N2O {1 / tp[2]:.4f}; steps Xe->Ar "
        f"{tp[0] / tp[1]:.4f}, Ar->N2O {tp[1] / tp[2]:.4f}, net "
        f"{tp[0] / tp[2]:.4f}; W/mN "
        + ", ".join(f"{w:.2f}" for w in tf._watts_per_mn(r))
    )

    base = baseline_at_anchor(PARAMS)
    save_pdf(tf.fig_ep_res_pt(base, **mk), out, "fig_ep_res_pt.pdf")


@task(
    "sl2022_validation",
    "Chapter 3 exit-velocity check against Sheppard "
    "and Little (2022): per-point factor k, pair offsets, fig_sl_validation",
)
def run_sl2022_validation(ctx):
    R = sl.per_point_k()

    section("Per-point k = u_e,meas / u_e,model(T_e)")
    for p, s in R["per_propellant"].items():
        print(
            f"  {p:4s} n {s['n']:2d}  mean k {s['mean']:.4f}  sd {s['sd']:.4f}  "
            f"se {s['se']:.4f}  T_e {s['Te_range'][0]:.2f}-{s['Te_range'][1]:.2f} eV  "
            f"dropped {s['dropped']}  groups {s['groups']}"
        )
    print(f"  model checks: {R['model_checks']}")

    section("Offsets, mean(k_num) / mean(k_den) - 1")
    for name, e in R["pairs"].items():
        print(
            f"  {name:7s} {100 * e['offset']:+7.2f} %  SE {100 * e['se']:.2f} %  "
            f"bootstrap 95 % {100 * e['boot_interval95'][0]:+.2f} to "
            f"{100 * e['boot_interval95'][1]:+.2f} %  cluster {e['n_groups']}: "
            f"{100 * e['cluster_boot_interval95'][0]:+.2f} to "
            f"{100 * e['cluster_boot_interval95'][1]:+.2f} %"
        )

    section("H2O / Ar offset, method variants")
    res, _tables, _check = sl.sensitivity_variants()
    for name, e in res.items():
        print(f"  {name}: {100 * e['offset']:+.2f} % (SE {100 * e['se']:.2f} %)")

    out = csv_dir(ctx)
    rows = [
        "propellant,mdot_mg_s,mdot_over_P_ug_per_J,ue_meas_m_s,Te_interp_eV,ue_model_m_s,k"
    ]
    for p, d in R["points"].items():
        for _, r in d.iterrows():
            rows.append(
                f"{p},{r.mdot_mg_s},{r.x!r},{r.ue_m_s!r},{r.Te_eV!r},{r.ue_model!r},{r.k!r}"
            )
    (out / "sl2022_perpoint_k.csv").write_text("\n".join(rows) + "\n", encoding="utf-8")
    (out / "sl2022_perpoint.json").write_text(
        json.dumps(
            dict(
                per_propellant=R["per_propellant"],
                pairs=R["pairs"],
                model_checks=R["model_checks"],
                h2o_ar_variants={
                    k: dict(offset=v["offset"], se=v["se"]) for k, v in res.items()
                },
            ),
            indent=1,
        ),
        encoding="utf-8",
    )

    fig_dir = figure_dir()
    pts = {
        p: (d["Te_eV"].to_numpy(), d["k"].to_numpy()) for p, d in R["points"].items()
    }
    save_pdf(
        tf.fig_sl_validation(pts, R["per_propellant"]), fig_dir, "fig_sl_validation.pdf"
    )
    print(f"  written {out} and {fig_dir / 'fig_sl_validation.pdf'}")
