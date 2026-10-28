"""
Reproduction of Petro and Sedwick (2017) Figs. 1, 2, 3, 5, 6 and 7 with this
work's water and argon model, against the digitised curves in data/ps2017/.
"""

# Writes the panels of the two chapter 3 figures to --out (default
# output/ps2017):
#
#     ps_repro_composition_a.pdf   rate constants (Fig. 1)
#     ps_repro_composition_b.pdf   negative-ion-to-neutral ratio (Fig. 2)
#     ps_repro_composition_c.pdf   positive-ion fractions (Fig. 3)
#     ps_repro_efficiency_a.pdf    ionisation cost against their fits (Fig. 6)
#     ps_repro_efficiency_b.pdf    fractional energy loss (Fig. 5)
#     ps_repro_efficiency_c.pdf    eta_T against T_e (Fig. 7a)
#     ps_repro_efficiency_d.pdf    eta_T against R_m (Fig. 7b)
#
# with the deviation tables ps_repro_deviation_<fig>.csv and
# ps_repro_deviation.md.
#
# Inputs that differ from the model, applied here on copies:
# 1. OH+ ionisation cross section x0.01, as read in Petro and Sedwick
#    (2017) Figs. 1 and 3; only the OH+ rate constant and the ion fractions
#    use it.
# 2. Positive-ion fractions from their Eq. (4), n_i/n = kappa_i sqrt(m_i) /
#    sum_j kappa_j sqrt(m_j); the model uses n_i/n = kappa_i / sum_j kappa_j.
# 3. Single-species H2O+ sheath drop, Eq. (6) with m_i = m_H2O.
# 4. Their water ionisation-cost fit, Eqs. (9) and (10) with Tables 2 and 3.
#
# Deviation metric, at each digitised x: relative for Figs. 1 and 2 (Fig. 1
# points below 1e-3 of the curve maximum excluded), absolute for Figs. 3, 5,
# 7a and 7b, and relative against the fits on a 0.01 eV grid for Fig. 6. A
# point is beyond digitising error when |dev| > 3 sigma, sigma being the
# curve's digitising error from its CSV header.
#
# Usage: python scripts/ps2017_reproduction.py [--out DIR] [--diagnostics]

import argparse
import csv
import re
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from multimode.plotting.style import (  # noqa: E402
    PROPELLANT_STYLE,
    PROPELLANT_SHADE,
    CATEGORY_COLOURS,
)
from multimode.thruster.ionisation_cost import fractional_loss_by_category  # noqa: E402
from multimode.thruster.cross_sections import (  # noqa: E402
    compute_rate_constants_vs_Te,
    compute_ion_fractions,
    compute_ion_fractions_species_resolved,
    compute_detach_rate_constants_vs_Te,
    sigma_detach_Hminus,
)
from multimode.thruster.ion_population import compute_neg_ion_ratios  # noqa: E402
from multimode.thruster.ionisation_cost import solve_Ec, build_h2o_channels_full  # noqa: E402
from multimode.thruster.propellants import AR  # noqa: E402
from multimode.thruster.propellant_data import h2o as H2O_DATA  # noqa: E402
from multimode.thruster.exhaust import sheath_potential_PS  # noqa: E402
from multimode.thruster.thrust_efficiency import thrust_efficiency, PS_Ar_curve_fit  # noqa: E402
from multimode.constants import M_AR  # noqa: E402

PS_DIR = REPO / "data" / "ps2017"
OH_SCALE_PS = 0.01  # OH+ cross-section factor of their Figs. 1 and 3
SIGMA_K = 3.0  # |dev| > SIGMA_K * digitising sigma = beyond error
FIG1_FLOOR = 1e-3  # Fig. 1 threshold-region exclusion (x curve max)

TE = np.linspace(1.0, 10.0, 181)  # model line grid [eV]
RM = np.geomspace(1.0, 100.0, 121)  # model line grid, Fig. 7b
TE_FIG7B = (2.0, 5.0, 10.0)
RM_FIG7A = (1.0, 10.0, 100.0)


# Water E_c fit, Petro and Sedwick Eqs. (9)-(10), Tables 2 and 3

PS_WATER_LOW = dict(c0=0.4202, c1=13.72, c2=-0.9948, c3=-9.343e-3, c4=7.115, c5=52.23)
PS_WATER_HIGH = dict(k0=38.54, k1=1.989, k2=-18.19, k3=2597.0, k4=-6.865)


def ps_water_curve_fit(Te):
    """Water E_c fit [eV], Eqs. (9)-(10); low branch for Te <= 3 eV [eV]."""
    Te = np.asarray(Te, dtype=float)
    out = np.zeros_like(Te)
    low, high = Te <= 3.0, Te > 3.0
    c, k = PS_WATER_LOW, PS_WATER_HIGH
    out[low] = (
        c["c0"] * np.exp(c["c1"] * Te[low] ** c["c2"])
        + c["c3"] * Te[low] ** c["c4"]
        + c["c5"]
    )
    t = Te[high]
    out[high] = (
        k["k0"] * t / (t - k["k1"])
        + k["k2"] * t / (t + k["k1"])
        + (t - k["k3"]) / (t - k["k4"]) ** 2
    )
    return out


# Species bookkeeping

_ASCII = str.maketrans("₂⁺⁻", "2+-")  # subscript 2, sup +, sup -


def token(sp):
    """'H₂O⁺' -> 'H2Oplus' (the canonical-file species token)."""
    return sp.translate(_ASCII).replace("+", "plus").replace("-", "minus")


SPECIES_TEX = {
    "H2Oplus": r"H$_2$O$^+$",
    "OHplus": r"OH$^+$",
    "Hplus": r"H$^+$",
    "Oplus": r"O$^+$",
    "H2plus": r"H$_2^+$",
    "Oplusplus": r"O$^{++}$",
    "Hminus": r"H$^-$",
    "Ominus": r"O$^-$",
    "OHminus": r"OH$^-$",
}
# One colour per species across the composition panels
SPECIES_COLOURS = {
    "H2Oplus": "#1f77b4",
    "Hplus": "#2ca02c",
    "Oplus": "#d62728",
    "OHplus": "#9467bd",
    "H2plus": "#8c564b",
    "Oplusplus": "#7a7a00",
    "Hminus": "#e377c2",
    "Ominus": "#0e7f8a",
    "OHminus": "#ff7f0e",
}
TOTAL_COLOUR = "#222222"
SINGLE_DETACH_COLOUR = "0.6"  # single-detachment-rate line

# Propellant colours from style.py, a second shade for the fits; argon
# solid, water dashed
PROP_COLOURS = {p: PROPELLANT_STYLE[p]["color"] for p in ("Ar", "H2O")}
PROP_SHADES = {p: PROPELLANT_SHADE[p] for p in ("Ar", "H2O")}
PROP_LS = {"Ar": "-", "H2O": "--"}
PROP_MARKERS = {"Ar": "o", "H2O": "s"}
PROP_TEX = {"Ar": "Ar", "H2O": r"H$_2$O"}

# Fig. 5 categories: file token, model category, label
FIG5_CATS = (
    ("ionisation", "ionisation", "ionisation"),
    ("electronic_excitation", "electronic excitation", "electronic\nexcitation"),
    ("vibration", "vibration", "vibration"),
    ("rotation", "rotation", "rotation"),
    ("elastic", "elastic", "elastic"),
    ("attachment", "attachment", "attachment"),
)
CAT_COLOURS = dict(CATEGORY_COLOURS)


# Model quantities, with the inputs above applied on copies


def ion_channels(oh_scale):
    """Copy of the water positive-ion channels with OH+ scaled by oh_scale."""
    out = {}
    for sp, d in H2O_DATA.H2O_ION_CHANNELS.items():
        s = np.array(d["sigma"], dtype=float, copy=True)
        if token(sp) == "OHplus":
            s = s * oh_scale
        out[sp] = {"sigma": s, "mass_amu": d["mass_amu"]}
    return out


def kappa_pos(Te, oh_scale):
    """{token: kappa_i(Te)} for the six positive ions [m^3/s]."""
    Te = np.atleast_1d(np.asarray(Te, dtype=float))
    return {
        token(sp): compute_rate_constants_vs_Te(
            H2O_DATA.H2O_ION_ENERGIES, ch["sigma"], Te
        )
        for sp, ch in ion_channels(oh_scale).items()
    }


def _kappa_neg_raw(Te):
    """{unicode species: kappa_attach(Te)}, as compute_neg_ion_ratios expects."""
    Te = np.atleast_1d(np.asarray(Te, dtype=float))
    return {
        sp: compute_rate_constants_vs_Te(d["energies"], d["sigma"], Te)
        for sp, d in H2O_DATA.H2O_NEG_ION_CHANNELS.items()
    }


def kappa_neg(Te):
    return {token(sp): k for sp, k in _kappa_neg_raw(Te).items()}


def ion_fractions(Te, oh_scale, closure):
    """closure 'eq4': Eq. (4); 'common': the model's kappa ratio."""
    Te = np.atleast_1d(np.asarray(Te, dtype=float))
    chans = ion_channels(oh_scale)
    kap = {
        sp: compute_rate_constants_vs_Te(H2O_DATA.H2O_ION_ENERGIES, ch["sigma"], Te)
        for sp, ch in chans.items()
    }
    if closure == "eq4":
        f = compute_ion_fractions_species_resolved(
            kap, {sp: ch["mass_amu"] for sp, ch in chans.items()}
        )
    elif closure == "common":
        f = compute_ion_fractions(kap)
    else:
        raise ValueError(closure)
    return {token(sp): v for sp, v in f.items()}


def negion_ratio(Te):
    """({token: ratio}, total) via Eq. (3) (H- and O- only)."""
    ind, tot = compute_neg_ion_ratios(_kappa_neg_raw(Te), np.atleast_1d(Te))
    return {token(sp): v for sp, v in ind.items()}, tot


def single_detach_ratio(Te, species=None):
    """Eq. (3) with one detachment rate: total attachment over the H-
    detachment rate. `species` limits the attachment sum (default H-, O-
    and OH-)."""
    Te = np.atleast_1d(np.asarray(Te, dtype=float))
    katt = _kappa_neg_raw(Te)
    att = sum(k for sp, k in katt.items() if species is None or token(sp) in species)
    return att / compute_detach_rate_constants_vs_Te(sigma_detach_Hminus, Te)


# Eq. (3) with the Fridman (2021) Fig. 2.9 detachment cross sections
# Detachment rate through Eq. (1) from the digitised H- and O- curves:
# sigma = 0 below the lowest digitised energy, log-log interpolation inside,
# held at the last value above; cm^2 x 1e16 -> 1e-16 cm^2 on load.

FRIDMAN_DIR = (
    REPO / "data" / "cross_sections" / "n2o" / "digitised_raw" / "Fridman2021_fig2_9"
)


def fridman_curve(species):
    """(E [eV], sigma [cm^2]) of the digitised Fridman Fig. 2.9 curve,
    species 'Hminus' or 'Ominus'."""
    d = np.loadtxt(FRIDMAN_DIR / f"{species}_detach.csv", delimiter=",")
    return d[:, 0], d[:, 1]


def fridman_sigma(species):
    """sigma(E) in 1e-16 cm^2 from the Fridman curve (rules above)."""
    E, s = fridman_curve(species)
    lE, ls = np.log(E), np.log(s * 1e16)

    def f(Ee):
        Ee = np.asarray(Ee, dtype=float)
        out = np.zeros_like(Ee)
        inside = (Ee >= E[0]) & (Ee <= E[-1])
        out[inside] = np.exp(np.interp(np.log(Ee[inside]), lE, ls))
        out[Ee > E[-1]] = s[-1] * 1e16
        return out

    return f


def kappa_detach_fridman(Te, species):
    Te = np.atleast_1d(np.asarray(Te, dtype=float))
    return compute_detach_rate_constants_vs_Te(fridman_sigma(species), Te)


def fridman_ratio(Te, variant):
    """n-/n_n in the single-ratio form of Eq. (3) with Fridman detachment.
    variant 'H': kappa_attach(H- + O- + OH-) / kappa_det,H-(Fridman).
    variant 'HO': the same numerator over the H- and O- Fridman rates
    weighted by each species' share of the H- + O- attachment (OH- has no
    detachment curve; its attachment stays in the numerator)."""
    Te = np.atleast_1d(np.asarray(Te, dtype=float))
    katt = kappa_neg(Te)
    att = sum(katt.values())
    kd_h = kappa_detach_fridman(Te, "Hminus")
    if variant == "H":
        return att / kd_h
    kd_o = kappa_detach_fridman(Te, "Ominus")
    wh = katt["Hminus"] / (katt["Hminus"] + katt["Ominus"])
    return att / (wh * kd_h + (1.0 - wh) * kd_o)


def tail_share(Te, species):
    """Share of the Eq. (1) integral from beyond the highest digitised energy."""
    E, _ = fridman_curve(species)
    g = np.linspace(0.0, 1000.0, 10000)  # compute_detach_rate_constant grid
    y = g * fridman_sigma(species)(g) * np.exp(-g / Te)
    tot = np.trapezoid(y, g)
    m = g >= E[-1]
    return float(np.trapezoid(y[m], g[m]) / tot)


def fig2_metrics(curve):
    """Comparison with the digitised Fig. 2 curve: model / digitised at the
    digitised T_e (median, rms log10 deviation, max |ratio - 1| below
    2 eV), the model peak T_e (0.01 eV grid) and whether the model stays
    below 1e-3 over 1-10 eV."""
    ps = load_ps("PS2017_fig2_negion_total")
    q = curve(ps["x"]) / ps["y"]
    g = np.linspace(1.0, 10.0, 901)
    y = curve(g)
    lo = ps["x"] < 2.0
    return dict(
        median=float(np.median(q)),
        rms_log=float(np.sqrt(np.mean(np.log10(q) ** 2))),
        peak_Te=float(g[int(np.argmax(y))]),
        peak=float(y.max()),
        max_dev_below_2=float(np.max(np.abs(q[lo] - 1.0))),
        below_1e3=bool(np.all(y < 1e-3)),
        ps_peak_Te=float(ps["x"][int(np.argmax(ps["y"]))]),
    )


def fridman_detachment_report():
    """Fig. 2 metrics of the negative-ion ratio variants; printed."""
    curves = {
        "Fridman H- detachment only": lambda T: fridman_ratio(T, "H"),
        "Fridman H- and O- detachment, attachment-weighted": lambda T: fridman_ratio(
            T, "HO"
        ),
        "Janev single rate (grey line)": lambda T: single_detach_ratio(T),
        "this work, per-species sum (H- Janev, O- Vejby-Christensen)": lambda T: (
            negion_ratio(T)[1]
        ),
    }
    R = {name: fig2_metrics(c) for name, c in curves.items()}
    R["_tail_share_10eV"] = {sp: tail_share(10.0, sp) for sp in ("Hminus", "Ominus")}
    R["_weights_H_at"] = {
        f"{t:g}": float(
            (lambda k: k["Hminus"] / (k["Hminus"] + k["Ominus"]))(kappa_neg([t]))[0]
        )
        for t in (1.0, 2.0, 5.0, 10.0)
    }
    R["_OH_share_of_attachment_at"] = {
        f"{t:g}": float((lambda k: k["OHminus"] / sum(k.values()))(kappa_neg([t]))[0])
        for t in (1.0, 2.0, 5.0, 10.0)
    }
    for name, m in R.items():
        print(f"  {name}: {m}")
    return R


def Ec_water(Te):
    ch, ion = build_h2o_channels_full()
    Ec, K, _ = solve_Ec(ch, ion, np.atleast_1d(np.asarray(Te, dtype=float)))
    return Ec, K, ch


def Ec_argon(Te):
    ch, ion = AR.channel_builder()
    Ec, K, _ = solve_Ec(ch, ion, np.atleast_1d(np.asarray(Te, dtype=float)))
    return Ec


def chi_water(Te):
    Te = np.atleast_1d(np.asarray(Te, dtype=float))
    _, K, ch = Ec_water(Te)
    return fractional_loss_by_category(K, ch, Te)


MASS = {"Ar": M_AR, "H2O": H2O_DATA.M_H2O}


def eta_T(prop, Te, R_m):
    """Eq. (18), single-species sheath drop (Eq. 6) for both propellants."""
    Te = np.atleast_1d(np.asarray(Te, dtype=float))
    Ec = Ec_argon(Te) if prop == "Ar" else Ec_water(Te)[0]
    phi = sheath_potential_PS(Te, MASS[prop])
    return thrust_efficiency(Te, phi, Ec, R_m)


def eta_T_vs_Rm(prop, Te, R_m_array):
    Ec = (Ec_argon(Te) if prop == "Ar" else Ec_water(Te)[0])[0]
    phi = sheath_potential_PS(Te, MASS[prop])
    R = np.asarray(R_m_array, dtype=float)
    return np.array([thrust_efficiency(Te, phi, Ec, r) for r in R])


# Digitised curves

_ERR_RE = re.compile(r"interior points\): ([0-9.eE+-]+) (decades|\(absolute)")


def load_ps(name):
    """Digitised CSV -> dict(x, y, err, err_log), or None if absent."""
    path = PS_DIR / f"{name}.csv"
    if not path.exists():
        return None
    text = path.read_text(encoding="ascii")
    hdr = " ".join(ln[1:].strip() for ln in text.splitlines() if ln.startswith("#"))
    m = _ERR_RE.search(hdr)
    rows = [ln for ln in text.splitlines() if ln and not ln.startswith("#")]
    data = np.array([[float(v) for v in ln.split(",")] for ln in rows[1:]])
    return dict(
        x=data[:, 0],
        y=data[:, 1],
        name=name,
        err=float(m.group(1)) if m else float("nan"),
        err_log=bool(m and m.group(2) == "decades"),
    )


def marker_subset(x, n=16, logx=False):
    """Indices of about n digitised points spread evenly in axis space, for
    display only (the deviation tables use every point)."""
    u = np.log10(x) if logx else x
    targets = np.linspace(u.min(), u.max(), n)
    return np.unique([int(np.argmin(np.abs(u - t))) for t in targets])


# Deviations


def deviation(ps, model_y, metric, floor_frac=None):
    """Per-point deviation and the summary row. metric 'rel' or 'abs'."""
    x, y = ps["x"], ps["y"]
    use = np.ones(len(x), dtype=bool)
    if floor_frac is not None:
        use &= y >= floor_frac * y.max()
    if metric == "rel":
        use &= y != 0
        dev = np.full(len(x), np.nan)
        dev[use] = (model_y[use] - y[use]) / y[use]
        if ps["err_log"]:
            sig = np.full(len(x), 10.0 ** ps["err"] - 1.0)
        else:
            sig = np.where(y != 0, ps["err"] / np.abs(np.where(y != 0, y, 1.0)), np.inf)
    else:
        dev = np.where(use, model_y - y, np.nan)
        sig = np.full(len(x), ps["err"])
    a = np.abs(dev)
    ok = np.isfinite(a)
    k = int(np.nanargmax(a)) if ok.any() else 0
    beyond = ok & (a > SIGMA_K * sig)
    return dict(
        curve=ps["name"],
        n=int(ok.sum()),
        n_excl=int((~ok).sum()),
        max=float(a[k]) if ok.any() else float("nan"),
        median=float(np.median(a[ok])) if ok.any() else float("nan"),
        signed_median=float(np.median(dev[ok])) if ok.any() else float("nan"),
        x_at_max=float(x[k]) if ok.any() else float("nan"),
        err=ps["err"],
        err_log=ps["err_log"],
        metric=metric,
        n_beyond=int(beyond.sum()),
        beyond_ranges=_ranges(x, beyond),
        dev=dev,
        beyond=beyond,
    )


def _ranges(x, mask):
    """Contiguous x runs where mask is True -> 'a-b' strings."""
    out, i = [], 0
    while i < len(x):
        if mask[i]:
            j = i
            while j + 1 < len(x) and mask[j + 1]:
                j += 1
            out.append(f"{x[i]:.3g}" if i == j else f"{x[i]:.3g}-{x[j]:.3g}")
            i = j + 1
        else:
            i += 1
    return out


def _fmt(v, metric):
    if not np.isfinite(v):
        return "n/a"
    return f"{100 * v:.2f} %" if metric == "rel" else f"{v:.3g}"


def _err_txt(r):
    if not np.isfinite(r["err"]):
        return "n/a"
    if r["err_log"]:
        return f"{r['err']:.3g} dec ({100 * (10 ** r['err'] - 1):.2g} %)"
    return f"{r['err']:.3g}"


def write_tables(out, fig, rows, xname, note):
    """One CSV and one markdown block per figure."""
    with open(
        out / f"ps_repro_deviation_{fig}.csv", "w", newline="", encoding="ascii"
    ) as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(
            [
                "curve",
                "metric",
                "n_points",
                "n_excluded",
                "max_abs_dev",
                "median_abs_dev",
                "signed_median_dev",
                f"{xname}_at_max",
                "digitising_error",
                "digitising_error_units",
                f"n_beyond_{SIGMA_K:g}sigma",
                f"{xname}_ranges_beyond",
            ]
        )
        for r in rows:
            w.writerow(
                [
                    r["curve"],
                    "relative" if r["metric"] == "rel" else "absolute",
                    r["n"],
                    r["n_excl"],
                    f"{r['max']:.6g}",
                    f"{r['median']:.6g}",
                    f"{r['signed_median']:.6g}",
                    f"{r['x_at_max']:.6g}",
                    f"{r['err']:.6g}",
                    "decades" if r["err_log"] else "axis units",
                    r["n_beyond"],
                    " ".join(r["beyond_ranges"]),
                ]
            )
    md = [
        f"### {fig}",
        "",
        note,
        "",
        f"| curve | n (excl.) | max abs dev | median abs dev | signed median | "
        f"{xname} at max | digitising error | n > {SIGMA_K:g} sigma | {xname} beyond |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        md.append(
            f"| {r['curve']} | {r['n']} ({r['n_excl']}) | {_fmt(r['max'], r['metric'])} | "
            f"{_fmt(r['median'], r['metric'])} | {_fmt(r['signed_median'], r['metric'])} | "
            f"{r['x_at_max']:.3g} | {_err_txt(r)} | {r['n_beyond']} | "
            f"{', '.join(r['beyond_ranges']) or '-'} |"
        )
    return md + [""]


# Figure helpers (print size, no tight bbox so the page size is exact)

CM = 1 / 2.54
FULL_W, HALF_W = 16.0, 7.8
PRINT_RC = {
    "font.size": 9,
    "axes.labelsize": 9,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "axes.titlesize": 9,
    "lines.linewidth": 1.3,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.minor.width": 0.4,
    "ytick.minor.width": 0.4,
    "grid.linewidth": 0.4,
    "grid.alpha": 0.35,
    "savefig.bbox": None,
    "pdf.fonttype": 42,
    "legend.framealpha": 0.9,
    "legend.borderpad": 0.3,
    "legend.handlelength": 1.6,
}
INK = "0.15"
MK = dict(linestyle="none", mfc="none", ms=3.6, mew=0.8)


def new_fig(
    w_cm, h_cm, l=1.45, r=0.3, b=1.0, t=0.25, nrows=1, ratios=None, hspace=0.08
):
    fig = plt.figure(figsize=(w_cm * CM, h_cm * CM))
    gs = fig.add_gridspec(
        nrows,
        1,
        left=l / w_cm,
        right=1 - r / w_cm,
        bottom=b / h_cm,
        top=1 - t / h_cm,
        height_ratios=ratios,
        hspace=hspace,
    )
    axes = [fig.add_subplot(gs[0])]
    for i in range(1, nrows):
        axes.append(fig.add_subplot(gs[i], sharex=axes[0]))
    return fig, axes


def end_labels(ax, items, min_sep=0.075, x_text=1.03, top=0.99, bottom=0.01):
    """Labels in the right margin at each curve's right-hand end, spread to
    min_sep (axes fraction) with a thin leader in the series colour.

    items: list of (y_end_data, text, colour[, leader linestyle])."""
    to_ax = ax.transData + ax.transAxes.inverted()
    x_end = ax.get_xlim()[1]
    pts = []
    for it in items:
        y, text, colour = it[:3]
        if not np.isfinite(y):
            continue
        yf = to_ax.transform((x_end, y))[1]
        pts.append([yf, yf, text, colour, it[3] if len(it) > 3 else "-"])
    pts.sort(key=lambda p: p[0])
    lab = [min(max(p[0], bottom), top) for p in pts]
    for i in range(1, len(lab)):
        lab[i] = max(lab[i], lab[i - 1] + min_sep)
    if lab and lab[-1] > top:
        lab[-1] = top
        for i in range(len(lab) - 2, -1, -1):
            lab[i] = min(lab[i], lab[i + 1] - min_sep)
    for p, yl in zip(pts, lab):
        y0 = min(max(p[1], 0.0), 1.0)
        ax.plot(
            [1.0, x_text - 0.008],
            [y0, yl],
            transform=ax.transAxes,
            color=p[3],
            lw=0.7,
            ls=p[4],
            clip_on=False,
            solid_capstyle="butt",
        )
        ax.text(
            x_text,
            yl,
            p[2],
            transform=ax.transAxes,
            ha="left",
            va="center",
            color=INK,
            fontsize=9,
            linespacing=0.95,
        )


def save(fig, out, name):
    fig.savefig(out / name, bbox_inches=None)
    plt.close(fig)
    print(f"    wrote {name}")


def te_axis(ax):
    ax.set_xlim(1, 10)
    ax.set_xticks([1, 2, 4, 6, 8, 10])
    ax.set_xlabel(r"$T_e$ [eV]")


# Panels, deviations

PS_LABEL = "P & S (2017)"
PS_MK = "o"  # digitised points


def _ps_handle(label=PS_LABEL, marker=PS_MK, colour="0.3"):
    from matplotlib.lines import Line2D

    return Line2D([], [], color=colour, marker=marker, label=label, **MK)


def _line_handle(label, colour="0.3", ls="-", lw=None):
    from matplotlib.lines import Line2D

    kw = {} if lw is None else {"lw": lw}
    return Line2D([], [], color=colour, ls=ls, label=label, **kw)


# Negative ions dashed, positive ions solid
def _species_ls(sp):
    return "--" if sp.endswith("minus") else "-"


def fig1(out, missing):
    """Composition panel (a): rate constants, Fig. 1."""
    order = (
        "H2Oplus",
        "Hplus",
        "Oplus",
        "OHplus",
        "H2plus",
        "Oplusplus",
        "Hminus",
        "Ominus",
        "OHminus",
    )
    kp, kn = kappa_pos(TE, OH_SCALE_PS), kappa_neg(TE)
    kall = {**kp, **kn}
    fig, (ax,) = new_fig(FULL_W, 9.4, l=1.75, r=1.55, b=1.95)
    ax.set_yscale("log")
    rows, labels = [], []
    for sp in order:
        c = SPECIES_COLOURS[sp]
        ax.plot(TE, kall[sp], color=c, ls=_species_ls(sp))
        ps = load_ps(f"PS2017_fig1_k_{sp}")
        if ps is None:
            missing.append(f"PS2017_fig1_k_{sp}")
        else:
            i = marker_subset(ps["x"], 18)
            ax.plot(ps["x"][i], ps["y"][i], marker=PS_MK, mec=c, **MK)
            km = (
                kappa_pos(ps["x"], OH_SCALE_PS)
                if sp.endswith("plus")
                else kappa_neg(ps["x"])
            )[sp]
            rows.append(deviation(ps, km, "rel", floor_frac=FIG1_FLOOR))
        labels.append((kall[sp][-1], SPECIES_TEX[sp], c, _species_ls(sp)))
    te_axis(ax)
    ax.set_ylim(1e-22, 1e-14)
    ax.set_ylabel(r"$\kappa$ [m$^3$/s]")
    end_labels(ax, labels, min_sep=0.082, x_text=1.03)
    # Legend for panels (a) to (c)
    fig.legend(
        handles=[
            _line_handle("this work, positive ion"),
            _line_handle("this work, negative ion", ls="--"),
            _line_handle(
                "single detachment rate (P & S Eq. 3), panel (b)", SINGLE_DETACH_COLOUR
            ),
            _ps_handle(),
        ],
        loc="lower left",
        ncol=2,
        frameon=False,
        bbox_to_anchor=(1.0 / FULL_W, 0.0),
        borderaxespad=0.15,
    )
    save(fig, out, "ps_repro_composition_a.pdf")
    return rows


# Negative-ion ratio contributions: dashed H-, dash-dot O-
NEG_LS = {"Hminus": "--", "Ominus": (0, (6.0, 1.6, 1.4, 1.6))}


def fig2(out, missing):
    """Composition panel (b): negative-ion-to-neutral ratio, Fig. 2."""
    ind, tot = negion_ratio(TE)
    fig, (ax,) = new_fig(HALF_W, 6.2, l=1.25, r=0.95)
    s = 1e4
    ax.plot(TE, tot * s, color=TOTAL_COLOUR, ls="--")
    labels = [(tot[-1] * s, "total", TOTAL_COLOUR, "--")]
    handles = [_line_handle("this work", TOTAL_COLOUR, "--")]
    for sp in ("Hminus", "Ominus"):
        if sp in ind:
            ax.plot(TE, ind[sp] * s, color=SPECIES_COLOURS[sp], lw=1.0, ls=NEG_LS[sp])
            labels.append(
                (ind[sp][-1] * s, SPECIES_TEX[sp], SPECIES_COLOURS[sp], NEG_LS[sp])
            )
    # Ratio with a single detachment rate
    single = single_detach_ratio(TE)
    ax.plot(TE, single * s, color=SINGLE_DETACH_COLOUR, lw=1.0, zorder=1)
    handles.append(
        _line_handle("single detachment rate\n(P & S Eq. 3)", SINGLE_DETACH_COLOUR)
    )
    rows = []
    ps = load_ps("PS2017_fig2_negion_total")
    if ps is None:
        missing.append("PS2017_fig2_negion_total")
    else:
        i = marker_subset(ps["x"], 16)
        ax.plot(ps["x"][i], ps["y"][i] * s, marker=PS_MK, mec=TOTAL_COLOUR, **MK)
        rows.append(deviation(ps, negion_ratio(ps["x"])[1], "rel"))
        handles.append(_ps_handle(colour=TOTAL_COLOUR))
    te_axis(ax)
    if ps is not None:
        for name, sp in (("H-, O-, OH-", None), ("H-, O-", ("Hminus", "Ominus"))):
            q = single_detach_ratio(ps["x"], sp) / ps["y"]
            print(
                f"  single detachment rate, attachment {name}: model / P+S "
                f"median {np.median(q):.4f}, min {np.min(q):.4f}, "
                f"max {np.max(q):.4f}"
            )
    top = max(np.nanmax(tot), np.nanmax(single), np.nanmax(ps["y"]) if ps else 0) * s
    ax.set_ylim(0, np.ceil(top * 1.08) + 1.0)
    ax.set_yticks(np.arange(0, ax.get_ylim()[1] + 0.1, 2))
    # plotted values are 1e4 x the ratio
    ax.set_ylabel(r"$n_{i-}/n_n$ [$\times 10^{-4}$]")
    end_labels(ax, labels, min_sep=0.09)
    save(fig, out, "ps_repro_composition_b.pdf")
    return rows


def fig3(out, missing):
    """Composition panel (c): positive-ion fractions, Fig. 3, with x0.01
    OH+ and Eq. (4)."""
    f_ps = ion_fractions(TE, OH_SCALE_PS, "eq4")
    fig, (ax,) = new_fig(HALF_W, 6.2, l=1.6, r=1.35)
    ax.set_yscale("log")
    labels, rows = [], []
    for sp in ("H2Oplus", "Hplus", "Oplus", "OHplus"):
        c = SPECIES_COLOURS[sp]
        ax.plot(TE, f_ps[sp], color=c)
        labels.append((f_ps[sp][-1], SPECIES_TEX[sp], c))
        ps = load_ps(f"PS2017_fig3_frac_{sp}")
        if ps is None:
            missing.append(f"PS2017_fig3_frac_{sp}")
            continue
        i = marker_subset(ps["x"], 14)
        i = i[(ps["y"][i] >= 1e-4) & (ps["y"][i] <= 1.0)]  # drawn unclipped at the top
        ax.plot(ps["x"][i], ps["y"][i], marker=PS_MK, mec=c, clip_on=False, **MK)
        rows.append(
            deviation(ps, ion_fractions(ps["x"], OH_SCALE_PS, "eq4")[sp], "abs")
        )
    te_axis(ax)
    ax.set_ylim(1e-4, 1.0)
    ax.set_ylabel(r"$f_i$ [-]")
    end_labels(ax, labels, min_sep=0.075, top=1.0)
    save(fig, out, "ps_repro_composition_c.pdf")
    return rows


def fig6(out):
    """Efficiency panel (a): E_c against the fits of Fig. 6 (dotted)."""
    ec = {"Ar": Ec_argon(TE), "H2O": Ec_water(TE)[0]}
    fit = {"Ar": PS_Ar_curve_fit(TE), "H2O": ps_water_curve_fit(TE)}
    fig, (ax,) = new_fig(HALF_W, 6.6, l=1.35, r=0.3)
    ax.set_yscale("log")
    handles = []
    for p in ("Ar", "H2O"):
        ax.plot(TE, ec[p], color=PROP_COLOURS[p], ls=PROP_LS[p])
        ax.plot(TE, fit[p], color=PROP_SHADES[p], ls=":", lw=1.6)
        handles += [
            _line_handle(f"{PROP_TEX[p]}, this work", PROP_COLOURS[p], PROP_LS[p]),
            _line_handle(f"{PROP_TEX[p]}, P & S fit", PROP_SHADES[p], ":", lw=1.6),
        ]
    ax.set_ylim(10, 1e4)
    ax.set_ylabel(r"$E_c$ [eV]")
    te_axis(ax)
    ax.legend(handles=[handles[i] for i in (0, 2, 1, 3)], loc="upper right")
    save(fig, out, "ps_repro_efficiency_a.pdf")


# Fig. 5 marker per channel
FIG5_MARKERS = {
    "ionisation": "o",
    "electronic excitation": "s",
    "vibration": "^",
    "rotation": "v",
    "elastic": "D",
    "attachment": "P",
}
# Channels labelled at the line end; the rest are named in the legend only.
FIG5_INLINE = ("ionisation", "electronic excitation")


def fig5(out, missing):
    """Efficiency panel (b): fractional energy loss chi_j, Fig. 5."""
    chi = chi_water(TE)
    fig, (ax,) = new_fig(HALF_W, 6.6, l=1.25, r=1.95)
    labels, rows, handles = [], [], []
    for tok, cat, lab in FIG5_CATS:
        c = CAT_COLOURS[cat]
        mk = FIG5_MARKERS[cat]
        if cat in chi:
            ax.plot(TE, chi[cat], color=c)
            if cat in FIG5_INLINE:
                labels.append((chi[cat][-1], lab, c))
        handles.append(_line_handle(cat, c))
        handles[-1].set_marker(mk)
        handles[-1].set_markerfacecolor("none")
        handles[-1].set_markersize(MK["ms"])
        handles[-1].set_markeredgewidth(MK["mew"])
        ps = load_ps(f"PS2017_fig5_chi_{tok}")
        if ps is None:
            missing.append(f"PS2017_fig5_chi_{tok}")
            continue
        i = marker_subset(ps["x"], 14)
        ax.plot(ps["x"][i], ps["y"][i], marker=mk, mec=c, **MK)
        if cat in chi:
            rows.append(deviation(ps, chi_water(ps["x"])[cat], "abs"))
        else:
            print(f"    Fig. 5 {cat}: no model category; digitised markers only")
    te_axis(ax)
    ax.set_ylim(0, 0.7)
    ax.set_ylabel(r"$\chi_j$ [-]")
    end_labels(ax, labels, min_sep=0.085)
    minor = [h for h in handles if h.get_label() not in FIG5_INLINE]
    minor.append(_ps_handle())
    # legend corner at T_e = 5.4 eV, between the major and minor channels
    ax.legend(
        handles=minor,
        loc="upper left",
        bbox_to_anchor=(5.4, 0.37),
        bbox_transform=ax.transData,
        fontsize=8,
        labelspacing=0.25,
        handlelength=1.0,
        handletextpad=0.4,
        borderpad=0.25,
        borderaxespad=0.0,
    )
    save(fig, out, "ps_repro_efficiency_b.pdf")
    return rows


def _prop_legend(ax, loc):
    """Lines for this work, one entry for both digitised marker shapes."""
    from matplotlib.lines import Line2D
    from matplotlib.legend_handler import HandlerTuple

    h = [
        Line2D(
            [],
            [],
            color=PROP_COLOURS[p],
            ls=PROP_LS[p],
            label=f"{PROP_TEX[p]}, this work",
        )
        for p in ("Ar", "H2O")
    ]
    mk = tuple(
        Line2D([], [], color=PROP_COLOURS[p], marker=PROP_MARKERS[p], **MK)
        for p in ("Ar", "H2O")
    )
    ax.legend(
        handles=h + [mk],
        labels=[x.get_label() for x in h] + [PS_LABEL],
        handler_map={tuple: HandlerTuple(ndivide=None, pad=0.6)},
        loc=loc,
        frameon=True,
        fontsize=8,
        labelspacing=0.25,
    )


def fig7a(out, missing):
    """Efficiency panel (c): eta_T against T_e at R_m = 1, 10, 100."""
    fig, (ax,) = new_fig(HALF_W, 6.6, l=1.25, r=1.75)
    rows, labels = [], []
    for R in RM_FIG7A:
        ends = []
        for p in ("Ar", "H2O"):
            c = PROP_COLOURS[p]
            e = eta_T(p, TE, R)
            ax.plot(TE, e, color=c, ls=PROP_LS[p])
            ends.append(e[-1])
            name = f"PS2017_fig7a_etaT_{p}_Rm{R:g}"
            ps = load_ps(name)
            if ps is None:
                missing.append(name)
                continue
            i = marker_subset(ps["x"], 12)
            ax.plot(ps["x"][i], ps["y"][i], marker=PROP_MARKERS[p], mec=c, **MK)
            rows.append(deviation(ps, eta_T(p, ps["x"], R), "abs"))
        labels.append((float(np.mean(ends)), rf"$R_m = {R:g}$", "0.5"))
    te_axis(ax)
    ax.set_ylim(0, 0.6)
    ax.set_ylabel(r"$\eta_T$ [-]")
    _prop_legend(ax, "upper left")
    end_labels(ax, labels, min_sep=0.09)
    save(fig, out, "ps_repro_efficiency_c.pdf")
    return rows


def fig7b(out, missing):
    """Efficiency panel (d): eta_T against R_m at T_e = 2, 5, 10 eV."""
    fig, (ax,) = new_fig(HALF_W, 6.6, l=1.25, r=1.75, b=1.2)
    ax.set_xscale("log")
    rows, labels = [], []
    for T in TE_FIG7B:
        ends = []
        for p in ("Ar", "H2O"):
            c = PROP_COLOURS[p]
            e = eta_T_vs_Rm(p, T, RM)
            ax.plot(RM, e, color=c, ls=PROP_LS[p])
            ends.append(e[-1])
            name = f"PS2017_fig7b_etaT_{p}_Te{T:g}"
            ps = load_ps(name)
            if ps is None:
                missing.append(name)
                continue
            i = marker_subset(ps["x"], 12, logx=True)
            ax.plot(ps["x"][i], ps["y"][i], marker=PROP_MARKERS[p], mec=c, **MK)
            rows.append(deviation(ps, eta_T_vs_Rm(p, T, ps["x"]), "abs"))
        labels.append((float(np.mean(ends)), rf"$T_e = {T:g}$ eV", "0.5"))
    ax.set_xlim(1, 100)
    ax.set_xlabel(r"$R_m$ [-]")
    ax.set_ylim(0, 0.6)
    ax.set_ylabel(r"$\eta_T$ [-]")
    end_labels(ax, labels, min_sep=0.09)
    save(fig, out, "ps_repro_efficiency_d.pdf")
    return rows


def ec_table(out):
    """Model E_c against the fits, relative, over 1-3 eV and 3-10 eV."""
    g = np.round(np.arange(1.0, 10.0 + 1e-9, 0.01), 2)
    ec = {"Ar": Ec_argon(g), "H2O": Ec_water(g)[0]}
    fit = {"Ar": PS_Ar_curve_fit(g), "H2O": ps_water_curve_fit(g)}
    md = [
        "### Fig. 6 (E_c): model against the fit, relative, 0.01 eV grid",
        "",
        "| propellant | T_e range | n | max abs dev | median abs dev | signed median | "
        "T_e at max |",
        "|---|---|---|---|---|---|---|",
    ]
    rows = []
    for p in ("Ar", "H2O"):
        dev = (ec[p] - fit[p]) / fit[p]
        for lo, hi, m in ((1.0, 3.0, g <= 3.0), (3.0, 10.0, g > 3.0)):
            a = np.abs(dev[m])
            k = int(np.argmax(a))
            r = (
                p,
                f"{lo:g}-{hi:g} eV",
                int(m.sum()),
                a[k],
                np.median(a),
                np.median(dev[m]),
                g[m][k],
            )
            rows.append(r)
            md.append(
                f"| {p} | {r[1]} | {r[2]} | {100 * r[3]:.2f} % | {100 * r[4]:.2f} % | "
                f"{100 * r[5]:+.2f} % | {r[6]:.2f} |"
            )
    with open(
        out / "ps_repro_deviation_fig6_Ec.csv", "w", newline="", encoding="ascii"
    ) as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(
            [
                "propellant",
                "Te_range",
                "n_points",
                "max_abs_rel_dev",
                "median_abs_rel_dev",
                "signed_median_rel_dev",
                "Te_at_max",
            ]
        )
        for r in rows:
            w.writerow(
                [
                    r[0],
                    r[1],
                    r[2],
                    f"{r[3]:.6g}",
                    f"{r[4]:.6g}",
                    f"{r[5]:.6g}",
                    f"{r[6]:.6g}",
                ]
            )
    md += [
        "",
        "Spot values [eV] (quoted: Ar 43 at 3.2 eV, 25.2 at 10 eV; "
        "H2O 65 at 3.2 eV, 23.9 at 10 eV):",
        "",
        "| T_e | Ar model | Ar fit | H2O model | H2O fit |",
        "|---|---|---|---|---|",
    ]
    for t in (2.0, 3.0, 3.2, 5.0, 10.0):
        i = int(np.argmin(np.abs(g - t)))
        md.append(
            f"| {t:g} | {ec['Ar'][i]:.2f} | {fit['Ar'][i]:.2f} | "
            f"{ec['H2O'][i]:.2f} | {fit['H2O'][i]:.2f} |"
        )
    return md + [""]


# Printed values


def h2o_fraction_and_fit_numbers():
    """H2O+ fraction for the four input combinations, and fit checks."""
    print("\n== H2O+ fraction, four combinations")
    t = np.array([2.0, 5.0, 10.0])
    g = np.linspace(1.0, 10.0, 91)
    print(
        "    OH+ input      closure          Te=2      Te=5      Te=10     "
        "min/max over 1-10 eV"
    )
    for scale, sname in ((1.0, "tabulated"), (OH_SCALE_PS, "x0.01")):
        for clo, cname in (("common", "common-velocity"), ("eq4", "Eq. (4) sqrt(m)")):
            f = ion_fractions(t, scale, clo)["H2Oplus"]
            fg = ion_fractions(g, scale, clo)["H2Oplus"]
            print(
                f"    {sname:<14} {cname:<16} "
                + "".join(f"{v:<10.4f}" for v in f)
                + f"{fg.min():.4f} / {fg.max():.4f}"
            )

    print("\n== Water E_c fit check (Eqs. 9-10, Tables 2-3)")
    for T, quoted in ((3.2, 65.0), (10.0, 23.9), (100.0, 20.8)):
        v = ps_water_curve_fit(np.array([T]))[0]
        print(
            f"    Te = {T:g} eV: fit {v:.3f} eV; quoted {quoted:g} eV; "
            f"{100 * (v - quoted) / quoted:+.2f} %"
        )
    k = PS_WATER_HIGH
    print(f"    high branch limit Te -> inf: k0 + k2 = {k['k0'] + k['k2']:.3f} eV")
    for name, f in (("H2O", ps_water_curve_fit), ("Ar", PS_Ar_curve_fit)):
        lo = f(np.array([3.0]))[0]
        hi = f(np.array([3.0 + 1e-9]))[0]
        print(
            f"    {name} fit at 3 eV: low branch {lo:.3f}, high branch {hi:.3f} eV "
            f"(step {100 * (hi - lo) / lo:+.2f} %)"
        )
    print(
        f"    H2O fit low branch at 1 eV: {ps_water_curve_fit(np.array([1.0]))[0]:.4g} eV; "
        f"model {Ec_water(np.array([1.0]))[0][0]:.4g} eV"
    )


def diagnostics():
    """Pipeline multi-species water sheath drop against the single-species
    drop used in the figures."""
    from multimode.config import PARAMS
    from multimode.thruster.pipeline import run_pipeline

    print("\n== Multi-species water sheath vs single-species H2O+ (run_pipeline)")
    res = run_pipeline(
        params=PARAMS, propellants=("Xe", "Ar", "H2O"), anchor_propellant="Xe"
    )
    Te = res["Te"]
    w = res["H2O"]
    phi_ss = sheath_potential_PS(Te, H2O_DATA.M_H2O)
    print(
        "    Te   phi_A(multi-species) phi_single  eta_T R_m=1/10/100: "
        "multi-species | single-species"
    )
    for t in (2.0, 5.0, 10.0):
        i = int(np.argmin(np.abs(Te - t)))
        e_a = [
            thrust_efficiency(
                Te[i : i + 1], w["phi_sh_A"][i : i + 1], w["Ec"][i : i + 1], R
            )[0]
            for R in RM_FIG7A
        ]
        e_s = [
            thrust_efficiency(Te[i : i + 1], phi_ss[i : i + 1], w["Ec"][i : i + 1], R)[
                0
            ]
            for R in RM_FIG7A
        ]
        print(
            f"    {Te[i]:<4.1f} {w['phi_sh_A'][i]:<20.4f} {phi_ss[i]:<11.4f} "
            + " ".join(f"{v:.4f}" for v in e_a)
            + " | "
            + " ".join(f"{v:.4f}" for v in e_s)
        )
    m = (Te >= 1.0) & (Te <= 10.0)
    print(
        f"    max |phi_A - phi_single| / phi_single over 1-10 eV: "
        f"{np.max(np.abs(w['phi_sh_A'][m] - phi_ss[m]) / phi_ss[m]) * 100:.3f} %"
    )


# Main


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(REPO / "output" / "ps2017"))
    ap.add_argument("--diagnostics", action="store_true")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(PRINT_RC)

    missing = []
    print("== Panels")
    r1 = fig1(out, missing)
    r2 = fig2(out, missing)
    r3 = fig3(out, missing)
    fig6(out)
    r5 = fig5(out, missing)
    r7a = fig7a(out, missing)
    r7b = fig7b(out, missing)

    md = [
        "# Petro and Sedwick (2017) reproduction: deviation tables",
        "",
        "Written by scripts/ps2017_reproduction.py. Model evaluated at each "
        "digitised x. Digitising error from each curve's CSV header. "
        f'"Beyond" = |dev| > {SIGMA_K:g} x digitising sigma in the metric\'s units.',
        "",
    ]
    md += write_tables(
        out,
        "fig1",
        r1,
        "Te",
        "Relative deviation; points below 1e-3 of the digitised curve's "
        "maximum excluded (threshold region). OH+ with the x0.01 factor.",
    )
    md += write_tables(
        out,
        "fig2",
        r2,
        "Te",
        "Relative deviation; model total = H- + O- (no OH- detachment model).",
    )
    md += write_tables(
        out,
        "fig3",
        r3,
        "Te",
        "Absolute deviation (fraction); model with x0.01 OH+ and Eq. (4) closure.",
    )
    md += write_tables(out, "fig5", r5, "Te", "Absolute deviation (chi_j).")
    md += write_tables(
        out,
        "fig7a",
        r7a,
        "Te",
        "Absolute deviation (eta_T); single-species sheath drop.",
    )
    md += write_tables(out, "fig7b", r7b, "R_m", "Absolute deviation (eta_T).")
    md += ec_table(out)
    md += [
        "## Missing digitised curves",
        "",
        (", ".join(missing) if missing else "None."),
        "",
    ]
    (out / "ps_repro_deviation.md").write_text("\n".join(md), encoding="ascii")
    print("    wrote ps_repro_deviation.md and ps_repro_deviation_*.csv")
    print(f"\n== Missing digitised curves: {', '.join(missing) if missing else 'none'}")

    h2o_fraction_and_fit_numbers()
    print(
        "\n== Eq. (3) with Fridman detachment, against the digitised Fig. 2 (not drawn)"
    )
    fridman_detachment_report()
    if args.diagnostics:
        diagnostics()


if __name__ == "__main__":
    main()
