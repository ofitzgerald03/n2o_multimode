"""
Thesis figures at print size, one builder per figure, named after its
LaTeX label. The tasks in analyses/ compute the inputs and save each figure
as <label>.pdf.
"""

# Chapter 3: fig_sl_validation.
# Chapter 4: fig_ep_res_n2o_ions, fig_ep_res_neg_ion, fig_ep_res_ec,
# fig_ep_res_ec_linear, fig_ep_res_h1_categories, fig_ep_res_etaT_ratio,
# fig_ep_res_etaT_abs, fig_ep_res_isp, fig_ep_res_pt,
# fig_ep_res_decomposition.
# Chapter 6: fig_res_tornado, fig_res_twarn_pivot, fig_res_fep_geo,
# fig_res_min_eta_u, fig_res_ca_map, fig_res_wet_bars, fig_res_response.

from __future__ import annotations

import math

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FixedLocator

from matplotlib.legend_handler import HandlerTuple

from .style import (
    PROPELLANT_STYLE,
    ARCHITECTURE_COLOURS,
    ARCHITECTURE_LINESTYLES,
    ARCHITECTURE_NAMES,
)
from .thesis_print import (
    FULL_W,
    HALF_W,
    INK,
    log_ticks,
    mark_te,
    mark_te_stacked,
    new_figure,
    print_style,
    shade_beyond_measured,
    stacked_figure,
    te_line_handles,
    legend_below,
    below_margin,
)
from ..thruster.thrust_efficiency import thrust_efficiency


TE_LABEL = r"$T_e$ [eV]"
TE_XLIM = (1.0, 30.0)
TE_TICKS = (1, 5, 10, 15, 20, 25, 30)

# Negative-ion-to-neutral ratio below which the plasma is treated as
# electropositive
ELECTROPOSITIVE_RATIO_MAX = 1e-3

# Linestyles of the positive-ion species
SPECIES_LS = (
    "-",
    "--",
    (0, (6.0, 1.6, 1.4, 1.6)),
    (0, (1.2, 1.3)),
    (0, (3.0, 1.2, 1.0, 1.2, 1.0, 1.2)),
)

# Margins of an electron-temperature figure; the top margin carries the
# representative-range ruler
_TE_MARGINS = dict(l=1.75, r=0.3, b=1.05, t=0.62)

# Legends inside the frame
LEGEND_IN = dict(
    fontsize=7.5, labelspacing=0.25, handlelength=1.6, borderpad=0.3, borderaxespad=0.35
)
CRITERION_COLOUR = "#B2182B"

# Loss groups of hypothesis H1
H1_ORDER = ("dissociation", "ionisation", "internal")
H1_COLOURS = {"dissociation": "#1f77b4", "ionisation": "#d62728", "internal": "#2ca02c"}
H1_LABELS = {
    "dissociation": "dissociation",
    "ionisation": "ionisation",
    "internal": "internal-energy losses",
}
H1_LS = {"dissociation": "-", "ionisation": "--", "internal": (0, (6.0, 1.6, 1.4, 1.6))}

# The swarm electronic channel counted as dissociation
H1_DISSOCIATION_CHANNEL = "Electronic: swarm 9.6 eV"


def _te_axis(ax, xlim=TE_XLIM, ticks=TE_TICKS, label=True):
    ax.set_xlim(*xlim)
    ax.set_xticks([t for t in ticks if xlim[0] <= t <= xlim[1]])
    if label:
        ax.set_xlabel(TE_LABEL)


def _style(label):
    return PROPELLANT_STYLE[label]


def _n2o_case(results):
    return results["params"].get("n2o_case", "B")


def _q(results, label):
    """Per-particle bundle in the closure of params['n2o_case']."""
    b = results[label]
    if _n2o_case(results) == "B" and b.get("q_B") is not None:
        return b["q_B"]
    return b["q_A"]


def crossings(x, y):
    """Linearly interpolated zero crossings of y(x)."""
    out = []
    for i in range(len(x) - 1):
        if y[i] == 0.0:
            out.append(float(x[i]))
        elif y[i] * y[i + 1] < 0.0:
            out.append(float(x[i] + y[i] / (y[i] - y[i + 1]) * (x[i + 1] - x[i])))
    return out


def _primary_phi(bundle, n2o_case="A"):
    """Sheath drop of a propellant in the closure n2o_case."""
    if bundle["phi_sh_B"] is not None and n2o_case == "B":
        return bundle["phi_sh_B"]
    return bundle["phi_sh_A"]


def _eta_T_at_Rm(bundle, Te, R_m, n2o_case="A"):
    """Thrust efficiency against Te at one mirror ratio."""
    return thrust_efficiency(
        Te,
        _primary_phi(bundle, n2o_case),
        bundle["Ec"],
        R_m,
    )


def _thrust_to_power_vs_xe(results):
    """Thrust per absorbed power at the anchor relative to Xe, for
    (Xe, Ar, N2O): R_T / R_Pin."""
    case = results["params"].get("n2o_case", "B")
    n2o_op = results["N2O"]["op_point_B" if case == "B" else "op_point_A"]
    ar_op = results["Ar"]["op_point_A"]
    tp_xe = 1.0
    tp_ar = ar_op["ratios"]["R_T"] / ar_op["ratios"]["R_Pin"]
    tp_n2o = n2o_op["ratios"]["R_T"] / n2o_op["ratios"]["R_Pin"]
    return tp_xe, tp_ar, tp_n2o


def _watts_per_mn(results):
    """Absorbed power per thrust at the anchor [W/mN] for (Xe, Ar, N2O)."""
    case = results["params"].get("n2o_case", "B")
    ops = (
        results["Xe"]["op_point_A"],
        results["Ar"]["op_point_A"],
        results["N2O"]["op_point_B" if case == "B" else "op_point_A"],
    )
    return tuple(op["P_in"] / (op["T"] * 1e3) for op in ops)


# Loss groups of hypothesis H1


def h1_groups(results):
    """Energy lost per ion by group, N2O and xenon.

    N2O: dissociation = the 9.6 eV swarm channel; ionisation = the five
    ionisation channels; internal = the 5.6 and 8.5 eV swarm channels,
    vibration, elastic and attachment. Xenon: ionisation = Xe+; internal =
    excitation and elastic; dissociation = 0.

    Returns (Te, n2o, xe, inc), each a {group: array}, inc = n2o - xe.
    """
    Te = results["Te"]

    def per_channel(p):
        b = results[p]
        return {
            c["name"]: (c["category"], b["K_by_name"][c["name"]] / b["kappa_ion"])
            for c in b["channels"]
        }

    n2o, xe = per_channel("N2O"), per_channel("Xe")
    zero = np.zeros_like(Te)
    N = {
        "dissociation": n2o[H1_DISSOCIATION_CHANNEL][1],
        "ionisation": sum((v for c, v in n2o.values() if c == "ionisation"), zero),
        "internal": sum(
            (
                v
                for k, (c, v) in n2o.items()
                if c in ("electronic", "vibration", "elastic", "attachment")
                and k != H1_DISSOCIATION_CHANNEL
            ),
            zero,
        ),
    }
    tot = N["dissociation"] + N["ionisation"] + N["internal"]
    assert np.allclose(tot, results["N2O"]["Ec"], rtol=1e-12)
    X = {
        "dissociation": zero,
        "ionisation": sum((v for c, v in xe.values() if c == "ionisation"), zero),
        "internal": sum(
            (v for c, v in xe.values() if c in ("electronic", "elastic")), zero
        ),
    }
    inc = {g: N[g] - X[g] for g in H1_ORDER}
    return Te, N, X, inc


# Chapter 4: plasma composition


def fig_ep_res_n2o_ions(analysis, *, te_boundary, te_ref, kappa_ylim=(1e-18, 1e-13)):
    """(a) Rate constants; (b) positive-ion fractions."""
    Te = analysis["Te"]
    with print_style():
        fig, (a1, a2) = stacked_figure(12.0, **_TE_MARGINS)
        cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
        sty = {
            sp: dict(color=cycle[i % len(cycle)], ls=SPECIES_LS[i % len(SPECIES_LS)])
            for i, sp in enumerate(analysis["kappa_pos"])
        }
        for sp, k in analysis["kappa_pos"].items():
            m = k > 0
            a1.semilogy(Te[m], k[m], label=sp, **sty[sp])
        for sp, k in analysis["kappa_neg"].items():
            m = k > 0
            a1.semilogy(
                Te[m],
                k[m],
                color="0.25",
                ls=(0, (4.0, 1.5)),
                lw=1.1,
                label=f"{sp} (attachment)",
            )
        a1.set_ylim(*kappa_ylim)
        log_ticks(a1)
        a1.set_ylabel(r"$\kappa$ [m$^3$/s]")

        for sp, f in analysis["fractions"].items():
            a2.plot(Te, f, **sty.get(sp, {}))
        _te_axis(a2)
        a2.set_ylim(0.0, 1.0)
        a2.set_ylabel(r"$f_i$ [-]")
        h, _ = a1.get_legend_handles_labels()
        h += te_line_handles(te_boundary, te_ref)
        a2.legend(handles=h, loc="upper right", ncol=2, columnspacing=1.0)
        mark_te_stacked((a1, a2), te_boundary=te_boundary, te_ref=te_ref)
    return fig


def fig_ep_res_neg_ion(
    analysis, total_ratio, *, criterion, te_boundary, te_ref, species_label
):
    """Negative-ion-to-neutral ratio with the electropositive criterion."""
    Te = analysis["Te"]
    with print_style():
        fig, (ax,) = new_figure(HALF_W, 6.4, **_TE_MARGINS)
        m = total_ratio > 0
        ax.semilogy(Te[m], total_ratio[m], color="black", lw=1.5, label=species_label)
        _te_axis(ax)
        ax.set_ylim(1e-5, 1e0)
        log_ticks(ax)
        ax.set_ylabel(r"$n_{i-}/n_n$ [-]")
        ax.axhline(
            criterion,
            color=CRITERION_COLOUR,
            lw=1.0,
            ls=(0, (6.0, 1.6, 1.4, 1.6)),
            zorder=1,
            label=f"criterion, $10^{{{math.log10(criterion):.0f}}}$",
        )
        h, _ = ax.get_legend_handles_labels()
        ax.legend(
            handles=h + te_line_handles(te_boundary, te_ref),
            loc="upper right",
            **LEGEND_IN,
        )
        mark_te(ax, te_boundary=te_boundary, te_ref=te_ref)
    return fig


# Chapter 4: ionisation cost

_EC_PROPS = ("Xe", "Ar", "N2O", "H2O")


def fig_ep_res_ec(results, *, te_boundary, te_ref, ylim=(1e1, 1e4)):
    """Ionisation cost of the four propellants against T_e, log y."""
    Te = results["Te"]
    with print_style():
        fig, (ax,) = new_figure(FULL_W, 7.0, **_TE_MARGINS)
        for p in _EC_PROPS:
            s = _style(p)
            ax.semilogy(
                Te,
                results[p]["Ec"],
                color=s["color"],
                ls=s["linestyle"],
                label=s["label"],
            )
        _te_axis(ax)
        ax.set_ylim(*ylim)
        log_ticks(ax)
        ax.set_ylabel(r"$E_c$ [eV]")
        h, _ = ax.get_legend_handles_labels()
        ax.legend(
            handles=h + te_line_handles(te_boundary, te_ref),
            loc="upper right",
            ncol=2,
            columnspacing=1.2,
        )
        mark_te(ax, te_boundary=te_boundary, te_ref=te_ref)
    return fig


def fig_ep_res_ec_linear(results, *, te_boundary, te_ref, ylim=(10.0, 100.0)):
    """Ionisation cost against T_e, linear y from 10 to 100 eV."""
    Te = results["Te"]
    with print_style():
        fig, (ax,) = new_figure(FULL_W, 7.0, **_TE_MARGINS)
        for p in _EC_PROPS:
            s = _style(p)
            ax.plot(
                Te,
                results[p]["Ec"],
                color=s["color"],
                ls=s["linestyle"],
                label=s["label"],
            )
        _te_axis(ax)
        ax.set_ylim(*ylim)
        ax.set_ylabel(r"$E_c$ [eV]")
        mark_te(ax, te_boundary=te_boundary, te_ref=te_ref, ruler=False)
    return fig


def fig_ep_res_h1_categories(
    results,
    *,
    te_boundary,
    te_ref,
    ylim_abs=(0.0, 120.0),
    ylim_inc=(-5.0, 60.0),
    mark_crossings=True,
):
    """(a) N2O energy lost per ion by group, stacked (top = E_c);
    (b) each group as an increment over xenon."""
    Te, N, X, inc = h1_groups(results)
    with print_style():
        fig, (a1, a2) = stacked_figure(11.5, **dict(_TE_MARGINS, l=1.95))
        base = np.zeros_like(Te)
        for g in H1_ORDER:
            a1.fill_between(
                Te, base, base + N[g], color=H1_COLOURS[g], alpha=0.55, lw=0
            )
            base = base + N[g]
            a1.plot(Te, base, color=H1_COLOURS[g], lw=0.8)
        a1.plot(Te, results["N2O"]["Ec"], color="black", lw=1.0)
        a1.set_ylim(*ylim_abs)
        a1.set_ylabel("Energy lost per ion [eV]")

        for g in H1_ORDER:
            a2.plot(Te, inc[g], color=H1_COLOURS[g], ls=H1_LS[g])
        a2.axhline(0.0, color="0.2", lw=0.7, zorder=1.6)
        _te_axis(a2)
        a2.set_ylim(*ylim_inc)
        a2.annotate(
            "xenon baseline",
            xy=(TE_XLIM[0], 0.0),
            xytext=(3.0, -1.5),
            textcoords="offset points",
            ha="left",
            va="top",
            color="0.2",
            fontsize=8,
        )
        a2.set_ylabel("N$_2$O increment over xenon\n[eV per ion]")
        handles = [
            (
                Patch(color=H1_COLOURS[g], alpha=0.55, lw=0),
                Line2D([], [], color=H1_COLOURS[g], ls=H1_LS[g]),
            )
            for g in H1_ORDER
        ]
        handles.append(Line2D([], [], color="black", lw=1.0))
        lines = te_line_handles(te_boundary, te_ref)
        a1.legend(
            handles + lines,
            [H1_LABELS[g] for g in H1_ORDER]
            + [r"total, $E_c$"]
            + [x.get_label() for x in lines],
            loc="upper right",
            ncol=2,
            handler_map={tuple: HandlerTuple(ndivide=None, pad=0.4)},
            handlelength=3.2,
        )
        if mark_crossings:
            # Crossings of the dissociation and ionisation increments with
            # the internal increment
            for ya, yb in (
                (inc["dissociation"], inc["internal"]),
                (inc["ionisation"], inc["internal"]),
            ):
                for xc in crossings(Te, ya - yb):
                    if TE_XLIM[0] <= xc <= TE_XLIM[1]:
                        yc = float(np.interp(xc, Te, ya))
                        if ylim_inc[0] <= yc <= ylim_inc[1]:
                            a2.plot(
                                xc,
                                yc,
                                "o",
                                mfc="white",
                                mec="0.35",
                                mew=0.8,
                                ms=4.0,
                                zorder=5,
                            )
        mark_te_stacked((a1, a2), te_boundary=te_boundary, te_ref=te_ref)
    return fig


# Chapter 4: thrust efficiency and specific impulse


def eta_T_ratio_band(results, num, den, rm_values):
    """eta_T(num)/eta_T(den) at each R_m in rm_values: {R_m: array}."""
    Te = results["Te"]
    case = _n2o_case(results)
    out = {}
    for R in rm_values:
        en = _eta_T_at_Rm(results[num], Te, R, case)
        ed = _eta_T_at_Rm(results[den], Te, R, case)
        with np.errstate(divide="ignore", invalid="ignore"):
            out[R] = np.where(ed > 0, en / ed, np.nan)
    return out


def fig_ep_res_etaT_ratio(results, *, te_boundary, te_ref, rm_main=5.0, rm_band=None):
    """N2O/Xe and Ar/Xe thrust-efficiency ratios at R_m = rm_main, with the
    range over the mirror ratios in rm_band shaded."""
    Te = results["Te"]
    if rm_band is None:
        rm_band = tuple(
            np.unique(
                np.concatenate([[1, 2, 4, 5, 6, 8, 10], np.linspace(1.0, 10.0, 37)])
            )
        )
    with print_style():
        fig, (ax,) = new_figure(HALF_W, 6.6, **_TE_MARGINS)
        handles = []
        for num in ("N2O", "Ar"):
            s = _style(num)
            band = eta_T_ratio_band(results, num, "Xe", rm_band)
            arr = np.vstack(list(band.values()))
            ax.fill_between(
                Te,
                np.nanmin(arr, 0),
                np.nanmax(arr, 0),
                color=s["color"],
                alpha=0.18,
                lw=0,
            )
            main = eta_T_ratio_band(results, num, "Xe", (rm_main,))[rm_main]
            ax.plot(Te, main, color=s["color"], ls=s["linestyle"])
            handles.append(
                Line2D(
                    [],
                    [],
                    color=s["color"],
                    ls=s["linestyle"],
                    label=f"{s['label']} / Xe",
                )
            )
        ax.axhline(1.0, color="0.2", lw=0.7, ls=":", zorder=1.6)
        _te_axis(ax)
        ax.set_ylim(0.0, 1.05)
        ax.set_ylabel(r"$\eta_T/\eta_{T,\mathrm{Xe}}$ [-]")
        ax.legend(
            handles=handles + te_line_handles(te_boundary, te_ref),
            loc="lower right",
            **LEGEND_IN,
        )
        mark_te(ax, te_boundary=te_boundary, te_ref=te_ref)
    return fig


def fig_ep_res_etaT_abs(results, *, te_boundary, te_ref, rm=5.0):
    """Thrust efficiency of the four propellants against T_e at R_m = rm."""
    Te = results["Te"]
    case = _n2o_case(results)
    with print_style():
        n_leg = 6
        fig, (ax,) = new_figure(
            HALF_W,
            6.6 + below_margin(n_leg) - 1.05,
            **dict(_TE_MARGINS, b=below_margin(n_leg)),
        )
        for p in _EC_PROPS:
            s = _style(p)
            ax.plot(
                Te,
                _eta_T_at_Rm(results[p], Te, rm, case),
                color=s["color"],
                ls=s["linestyle"],
                label=s["label"],
            )
        _te_axis(ax)
        ax.set_ylim(0.0, 0.6)
        ax.set_ylabel(r"$\eta_T$ [-]")
        h, _ = ax.get_legend_handles_labels()
        legend_below(fig, ax, h + te_line_handles(te_boundary, te_ref))
        mark_te(ax, te_boundary=te_boundary, te_ref=te_ref)
    return fig


def fig_ep_res_isp(results, *, te_boundary, te_ref):
    """(a) Ion specific impulse of the four propellants; (b) N2O and Ar
    over Xe."""
    Te = results["Te"]
    with print_style():
        fig, (a1, a2) = stacked_figure(11.0, **_TE_MARGINS)
        for p in _EC_PROPS:
            s = _style(p)
            a1.plot(
                Te,
                _q(results, p)["I_sp"],
                color=s["color"],
                ls=s["linestyle"],
                label=s["label"],
            )
        a1.set_ylim(0, 5500)
        a1.set_yticks([0, 1000, 2000, 3000, 4000, 5000])
        a1.set_ylabel(r"$I_{sp,\mathrm{ion}}$ [s]")
        h, _ = a1.get_legend_handles_labels()
        a1.legend(
            handles=h + te_line_handles(te_boundary, te_ref),
            loc="upper left",
            ncol=2,
            columnspacing=1.0,
        )
        xe = _q(results, "Xe")["I_sp"]
        for p in ("N2O", "Ar"):
            s = _style(p)
            a2.plot(
                Te, _q(results, p)["I_sp"] / xe, color=s["color"], ls=s["linestyle"]
            )
        _te_axis(a2)
        a2.set_ylim(1.6, 1.86)
        a2.set_ylabel(r"$I_{sp,\mathrm{ion}}/I_{sp,\mathrm{ion,Xe}}$ [-]")
        mark_te_stacked((a1, a2), te_boundary=te_boundary, te_ref=te_ref)
    return fig


def fig_ep_res_pt(baseline, *, te_boundary, te_ref, xlim=TE_XLIM, ylim=(1e1, 1e3)):
    """Absorbed power per unit thrust of Xe and N2O against T_e, log y."""
    base = baseline.get("baseline", baseline)
    Te = base["Te_grid"]
    with print_style():
        fig, (ax,) = new_figure(HALF_W, 6.4, **_TE_MARGINS)
        for p in ("Xe", "N2O"):
            s = _style(p)
            pt = base[p]["P_in_arr"] / (base[p]["T_N_arr"] * 1e3)
            ax.semilogy(Te, pt, color=s["color"], ls=s["linestyle"], label=s["label"])
        ax.set_ylim(*ylim)
        log_ticks(ax)
        ax.grid(True, which="major", alpha=0.3)
        _te_axis(ax, xlim=xlim)
        ax.set_ylabel(r"$P_{\mathrm{abs}}/T$ [W/mN]")
        h, _ = ax.get_legend_handles_labels()
        ax.legend(
            handles=h + te_line_handles(te_boundary, te_ref),
            loc="upper right",
            **LEGEND_IN,
        )
        mark_te(ax, te_boundary=te_boundary, te_ref=te_ref)
    return fig


def fig_ep_res_decomposition(results):
    """Absorbed power per unit thrust relative to Xe for Xe, Ar and N2O,
    with the Xe -> Ar and Ar -> N2O factors and the net factor, and each
    value in W/mN above its bar."""
    tp = _thrust_to_power_vs_xe(results)  # (Xe, Ar, N2O)
    pt = tuple(1.0 / v for v in tp)  # power per thrust vs Xe
    wmn = _watts_per_mn(results)
    labels = ("Xe", "Ar", "N2O")
    c = [_style(p)["color"] for p in labels]
    with print_style():
        fig, (ax,) = new_figure(HALF_W, 7.6, l=1.45, r=0.15, b=0.75, t=0.2)
        x = (0.0, 1.0, 2.15)
        w = 0.56
        xa = 0.5 * (x[1] + x[2])  # Ar -> N2O arrow
        for xi, v, col in zip(x, pt, c):
            ax.bar(
                xi,
                v,
                width=w,
                color=col,
                alpha=0.85,
                edgecolor="black",
                linewidth=0.6,
                zorder=3,
            )
        for xi, v, p in zip(x, pt, wmn):
            ax.text(
                xi,
                v + 0.025,
                f"{p:.2f}\nW/mN",
                ha="center",
                va="bottom",
                color=INK,
                linespacing=1.0,
                zorder=5,
            )
        guide = dict(color="0.65", lw=0.7, ls=(0, (4, 3)), zorder=2)
        ax.plot([w / 2, 3.15], [pt[0]] * 2, **guide)
        ax.plot([x[1] + w / 2, xa + 0.06], [pt[1]] * 2, **guide)
        ax.plot([x[2] + w / 2, 3.15], [pt[2]] * 2, **guide)

        def arrow(xa, y0, y1, col):
            ax.annotate(
                "",
                xy=(xa, y1),
                xytext=(xa, y0),
                arrowprops=dict(
                    arrowstyle="<->",
                    color=col,
                    lw=1.1,
                    shrinkA=0,
                    shrinkB=0,
                    mutation_scale=6,
                ),
                zorder=4,
            )

        def box(xt, yt, txt):
            ax.text(
                xt,
                yt,
                txt,
                ha="center",
                va="center",
                color=INK,
                linespacing=1.05,
                zorder=6,
                bbox=dict(
                    boxstyle="round,pad=0.2", fc="white", ec="0.8", lw=0.5, alpha=0.95
                ),
            )

        f1, f2, fn = pt[1] / pt[0], pt[2] / pt[1], pt[2] / pt[0]
        arrow(0.5, pt[0], pt[1], "0.2")
        box(0.5, 0.5 * (pt[0] + pt[1]), f"$\\times${f1:.3f}")
        arrow(xa, pt[1], pt[2], c[2])
        box(xa, pt[1] - 0.2, f"$\\times${f2:.3f}")
        ax.plot([xa, xa], [pt[1] - 0.13, pt[1]], color=c[2], lw=0.6, zorder=4)
        arrow(2.97, pt[0], pt[2], "0.2")
        box(2.97, 0.5 * (pt[0] + pt[2]), f"net\n$\\times${fn:.3f}")
        ax.set_xticks(x)
        ax.set_xticklabels([_style(p)["label"] for p in labels])
        ax.set_xlim(-0.42, 3.45)
        ax.set_ylim(0.0, 2.45)
        ax.set_yticks([0.0, 0.5, 1.0, 1.5, 2.0])
        ax.set_ylabel(
            r"$(P_\mathrm{abs}/T)$ / $(P_\mathrm{abs}/T)_\mathrm{Xe}$ [-]", fontsize=9
        )
        ax.grid(axis="x", visible=False)
        ax.tick_params(axis="x", length=0)
    return fig


# Chapter 6: spacecraft level

ARCH_NAMES = ARCHITECTURE_NAMES
MISSION_NAMES = {"leo": "LEO mission", "geo": "GEO transfer mission"}
WET_LABEL = "Propulsion wet mass [kg]"


def _arch_colour(key):
    return ARCHITECTURE_COLOURS[key]


def _arch_ls(key):
    return ARCHITECTURE_LINESTYLES[key]


def _finite(records):
    return [r for r in records if not r.get("diverged") and not math.isnan(r["wet"])]


# one-at-a-time sensitivity (tornado)

# Row labels of the tornado, by parameter name
TORNADO_ROW_LABELS = {
    "pf_metal": "Performance factor, metal",
    "m_min_tank": "Minimum tank mass",
    "m_fixed_tank": "Fixed mass per tank",
    "m_head_ecr": "Head mass, electrodeless",
    "m_head_mono": "Monopropellant head mass",
    "isp_mono": "Monopropellant specific impulse",
    "p_heater_mono": "Catalyst heater power",
    "eta_dcrf": "PPU efficiency, electrodeless",
    "sm_ppu": "PPU specific mass, Hall",
    "sm_rf": "PPU specific mass, electrodeless",
    "sm_eps_inc": "Specific mass, added power system",
    "power_margin_150": "Power margin, LEO bus",
    "power_margin_300": "Power margin, 300 kg bus",
    "bus_mass_leo": "Bus mass, LEO mission",
    "bus_mass_geo": "Bus mass, GEO transfer mission",
    "mission_life_leo": "Mission duration, LEO",
    "dv_ca_budget_1d": "Avoidance budget per event",
    "dv_eol_leo": "Disposal budget, LEO",
    "eol_leo_ep_loss_ratio": "Disposal loss ratio, LEO",
    "dv_transfer_chem": "Transfer budget, chemical",
    "eta_planning_transfer": "Transfer planning efficiency",
    "dv_sk_leo_rate": "Station-keeping rate, LEO",
    "w_duty_transfer": "Duty factor, transfer",
    "w_duty_leo": "Duty factor, LEO classes",
    "eps_cap_frac": "Power-system limit",
    "envelope_frac": "Bus mass-fraction limit",
    "thruster_R_Isp_ue": "Exhaust-velocity band",
    "thruster_R_T": "Thrust band",
    "thruster_R_Pin": "Absorbed-power band",
}


def _rows(tor, mission_key, metric, base):
    """Tornado rows sorted by their largest deviation from the baseline."""
    rows = []
    for bar in tor["bars"]:
        devs, err = [], None
        for e in ("lo", "hi"):
            d = bar["edges"][e]
            if d["error"]:
                err = d["error"]
                continue
            v = d[mission_key][metric]
            if not math.isnan(v):
                devs.append(v)
        rows.append(
            dict(
                name=bar["name"],
                group=bar["group"],
                devs=devs,
                err=err,
                span=max((abs(v - base) for v in devs), default=0.0),
                flips=any(v * base < 0 for v in devs),
            )
        )
    rows.sort(key=lambda r: r["span"])
    return rows


def tornado_rows(tor, mission_key, *, drop=()):
    """(baseline, rows) in plotting order (bottom first), without the rows
    in `drop`."""
    metric = "d_mm_se" if mission_key == "leo" else "d_mm_hy"
    base = tor["baseline"][mission_key][metric]
    rows = [r for r in _rows(tor, mission_key, metric, base) if r["name"] not in drop]
    return base, rows


def fig_res_tornado(tor, mission_key, *, row_labels, drop=()):
    """One-at-a-time sensitivity of the wet-mass difference between the
    two leading architectures of a mission.

    row_labels maps parameter names to row labels; rows in `drop` are left
    out. Rows are numbered from the top; a bar crossing zero is hatched
    and its label marked *.
    """
    base, rows = tornado_rows(tor, mission_key, drop=drop)
    missing = [r["name"] for r in rows if r["name"] not in row_labels]
    if missing:
        raise KeyError(f"no report-language label for {missing}")
    unevaluated = [r["name"] for r in rows if r["err"]]
    if unevaluated:
        raise ValueError(f"rows with an unevaluated bound: {unevaluated}")
    c_sys, c_thr = "#1f77b4", "#9467bd"
    b_cm = below_margin(2)
    h_cm = 0.36 * len(rows) + 0.7 + b_cm
    l_cm = 7.3
    with print_style():
        fig, (ax,) = new_figure(FULL_W, h_cm, l=l_cm, r=0.3, b=b_cm, t=0.15)
        lo = base
        for i, r in enumerate(rows):
            col = c_thr if r["group"] == "thruster" else c_sys
            for v in r["devs"]:
                flip = v * base < 0
                lo = min(lo, v)
                ax.barh(
                    i,
                    v - base,
                    left=base,
                    height=0.62,
                    color=col,
                    alpha=0.75 if v >= base else 0.45,
                    edgecolor=("#d62728" if flip else "none"),
                    hatch=("////" if flip else None),
                    lw=0.8,
                )
        ax.axvline(base, color="k", lw=0.9, zorder=3)
        ax.axvline(0.0, color="#d62728", ls="--", lw=0.9, zorder=3)
        n = len(rows)
        ax.set_yticks(range(n))
        ax.set_yticklabels(
            [row_labels[r["name"]] + (" *" if r["flips"] else "") for r in rows]
        )
        # Row numbers, 1 at the top
        x_num = -(l_cm - 0.15) / (FULL_W - l_cm - 0.3)
        for i in range(n):
            ax.text(
                x_num,
                i,
                f"{n - i}",
                transform=ax.get_yaxis_transform(),
                ha="left",
                va="center",
                color=INK,
            )
        ax.set_ylim(-0.7, n - 0.3)
        ax.grid(axis="y", visible=False)
        ax.tick_params(axis="y", length=0)
        pair = (
            ("multimode", "single_elec")
            if mission_key == "leo"
            else ("multimode", "hybrid")
        )
        ax.set_xlabel(
            f"{ARCH_NAMES[pair[0]]} minus {ARCH_NAMES[pair[1]].lower()} wet mass [kg]"
        )
        handles = [
            Patch(color=c_sys, alpha=0.75, label="System-level parameter"),
            Patch(color=c_thr, alpha=0.75, label="Band on the electrodeless unit"),
            Patch(
                facecolor="none",
                edgecolor="#d62728",
                hatch="////",
                lw=0.8,
                label="Crosses zero (marked *)",
            ),
            Line2D([], [], color="k", lw=0.9, label="Baseline"),
            Line2D([], [], color="#d62728", ls="--", lw=0.9, label="Zero difference"),
        ]
        hi = max([base, 0.0] + [v for r in rows for v in r["devs"]])
        span = hi - lo
        ax.set_xlim(lo - 0.06 * span, hi + 0.06 * span)
        # Legend under the x label, from the left edge of the figure
        fig.legend(
            handles=handles,
            loc="upper left",
            ncol=3,
            columnspacing=1.4,
            frameon=False,
            borderaxespad=0.0,
            bbox_to_anchor=(0.15 / FULL_W, ax.get_position().y0 - 1.0 / h_cm),
        )
    return fig


# LEO warning-time pivot

_TWARN_ORDER = ("multimode", "mono_multimode", "single_elec", "hybrid", "single_chem")


def fig_res_twarn_pivot(sw, flip_band, *, nominal_h, open_label):
    """Wet mass against collision-avoidance warning time, five
    architectures; filled markers feasible, open markers infeasible
    (`open_label`); the single-electric limit dashed with its sensitivity
    range shaded; the nominal warning time dotted."""
    H = 3600.0
    grid_h = [t / H for t in sw["grid"]]
    with print_style():
        fig, (ax,) = new_figure(FULL_W, 8.2, l=1.45, r=0.3, b=1.0, t=0.2)
        for key in _TWARN_ORDER:
            if key not in sw["records"]:
                continue
            recs = sw["records"][key]
            col = _arch_colour(key)
            ax.plot(
                grid_h,
                [r["wet"] for r in recs],
                color=col,
                ls=_arch_ls(key),
                lw=1.2,
                label=ARCH_NAMES[key],
                zorder=2,
            )
            ok = [(h, r["wet"]) for h, r in zip(grid_h, recs) if r["feasible"]]
            bad = [(h, r["wet"]) for h, r in zip(grid_h, recs) if not r["feasible"]]
            if ok:
                ax.plot(*zip(*ok), "o", color=col, ms=2.6, zorder=3)
            if bad:
                ax.plot(
                    *zip(*bad), "o", mfc="white", mec=col, mew=0.8, ms=3.2, zorder=3
                )
        handles, labels = ax.get_legend_handles_labels()
        handles += [
            Line2D(
                [], [], ls="none", marker="o", color="0.35", ms=3.2, label="feasible"
            ),
            Line2D(
                [],
                [],
                ls="none",
                marker="o",
                mfc="white",
                mec="0.35",
                mew=0.8,
                ms=3.6,
                label=open_label,
            ),
        ]
        if sw["flip_s"] is not None:
            ax.axvline(sw["flip_s"] / H, color="k", ls="--", lw=0.9, zorder=1)
            handles.append(
                Line2D(
                    [], [], color="k", ls="--", lw=0.9, label="single-electric limit"
                )
            )
        if flip_band is not None and flip_band["band_lo_s"] is not None:
            ax.axvspan(
                flip_band["band_lo_s"] / H,
                flip_band["band_hi_s"] / H,
                color="0.55",
                alpha=0.16,
                lw=0,
                zorder=1,
            )
            handles.append(
                Patch(
                    color="0.55",
                    alpha=0.16,
                    lw=0,
                    label="single-electric limit, sensitivity range",
                )
            )
        ax.axvline(nominal_h, color="k", ls=":", lw=0.9, zorder=1)
        ax.set_xscale("log")
        ticks = [t for t in (2, 4, 6, 12, 24, 48, 72) if grid_h[0] <= t <= grid_h[-1]]
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"{t:g}" for t in ticks])
        ax.minorticks_off()
        ax.set_xlim(grid_h[0] * 0.95, grid_h[-1] * 1.05)
        ax.set_xlabel(r"$t_{\mathrm{warn}}$ [h]")
        ax.set_ylabel(WET_LABEL)
        lo_w = min(min(r["wet"] for r in sw["records"][k]) for k in sw["records"])
        hi_w = max(max(r["wet"] for r in sw["records"][k]) for k in sw["records"])
        ax.set_ylim(math.floor(lo_w) - 0.5, hi_w + 0.42 * (hi_w - lo_w) + 0.5)
        # 'nominal' label in the widest gap between curves
        at = {
            k: float(
                np.interp(
                    math.log(nominal_h),
                    np.log(grid_h),
                    [r["wet"] for r in sw["records"][k]],
                )
            )
            for k in sw["records"]
        }
        ys = sorted(at.values())
        gaps = [
            (ys[i + 1] - ys[i], 0.5 * (ys[i] + ys[i + 1])) for i in range(len(ys) - 1)
        ]
        y_nom = max(gaps)[1]
        ax.text(nominal_h * 1.04, y_nom, "nominal", ha="left", va="center", color=INK)
        ax.legend(
            handles=handles,
            loc="upper right",
            ncol=3,
            columnspacing=0.7,
            handlelength=1.6,
            handletextpad=0.5,
            fontsize=8.1,
        )
    return fig


# GEO allocation sweep

_FEP_ORDER = ("multimode", "mono_multimode", "hybrid", "single_elec", "single_chem")


def fig_res_fep_geo(sweep, *, power_label="over power-system limit"):
    """Propulsion wet mass against the electric allocation fraction, GEO
    transfer mission, log y; open squares mark points over the power-system
    limit."""
    with print_style():
        fig, (ax,) = new_figure(FULL_W, 8.6, l=1.85, r=0.3, b=1.0, t=0.2)
        derived = {
            k: sweep[k]["f_ep_derived"]
            for k in _FEP_ORDER
            if k in sweep
            and len(_finite(sweep[k]["records"])) > 1
            and sweep[k].get("f_ep_derived") is not None
        }
        shared = (
            min(derived.values())
            if derived and max(derived.values()) - min(derived.values()) < 1e-6
            else None
        )
        all_wet, drawn_pw = [], False
        for key in _FEP_ORDER:
            if key not in sweep:
                continue
            recs = _finite(sweep[key]["records"])
            col = _arch_colour(key)
            name = ARCH_NAMES[key]
            if not recs:
                ax.plot([], [], "-", color=col, label=f"{name} (no solution)")
                continue
            all_wet += [r["wet"] for r in recs]
            pw = [
                (r["f_ep"], r["wet"])
                for r in recs
                if r["excluded"] and any(x.startswith("power") for x in r["reasons"])
            ]
            rest = [
                (r["f_ep"], r["wet"]) for r in recs if (r["f_ep"], r["wet"]) not in pw
            ]
            if len(recs) > 1:
                ax.plot(
                    [r["f_ep"] for r in recs],
                    [r["wet"] for r in recs],
                    color=col,
                    ls=_arch_ls(key),
                    lw=1.2,
                    marker="o",
                    ms=0,
                    label=name,
                    zorder=2,
                )
                ms = 2.8
            else:
                ax.plot([], [], "o", ms=4.5, color=col, label=name)
                ms = 4.5
            if rest:
                ax.plot(*zip(*rest), "o", color=col, ms=ms, zorder=3)
            if pw:
                ax.plot(*zip(*pw), "s", mfc="white", mec=col, mew=0.9, ms=3.6, zorder=3)
                drawn_pw = True
            if shared is None and key in derived:
                ax.axvline(derived[key], color=col, ls=":", lw=0.9)
        if shared is not None:
            ax.axvline(shared, color="0.35", ls=":", lw=1.0, zorder=1)
        ax.set_yscale("log")
        lo, hi = min(all_wet), max(all_wet)
        ax.set_ylim(
            10 ** math.floor(math.log10(lo)), 10 ** (math.ceil(math.log10(hi)) + 0.25)
        )
        log_ticks(ax)
        ax.set_xlim(-0.03, 1.03)
        ax.set_xlabel(r"$f_{EP}$ (impulsive-equivalent) [-]")
        ax.set_ylabel(WET_LABEL)
        ax.grid(True, which="major", alpha=0.3)
        handles, _ = ax.get_legend_handles_labels()
        if drawn_pw:
            handles.append(
                Line2D(
                    [],
                    [],
                    ls="none",
                    marker="s",
                    mfc="white",
                    mec="0.35",
                    mew=0.9,
                    ms=3.6,
                    label=power_label,
                )
            )
        if derived:
            handles.append(
                Line2D(
                    [],
                    [],
                    ls=":",
                    color="0.35",
                    lw=1.0,
                    label="nominal"
                    + (f" ($f_{{EP}}$ = {shared:.3f})" if shared is not None else ""),
                )
            )
        ax.legend(handles=handles, loc="upper left", ncol=2, columnspacing=1.4)
    return fig


# minimum viable eta_u


def _eta_panel(
    ax,
    sw,
    eta_range,
    eta_measured,
    eta_hi,
    *,
    log_y,
    ylim,
    water_range,
    notes,
    water_label_y,
):
    refs = sw["references"]
    ref_style = {
        "multimode_f0": (
            "0.35",
            (0, (4, 1.5, 1, 1.5)),
            "Biprop. multimode, chemical mode only",
        ),
        "single_chem": (
            _arch_colour("single_chem"),
            _arch_ls("single_chem"),
            ARCH_NAMES["single_chem"],
        ),
        "single_elec": (
            _arch_colour("single_elec"),
            _arch_ls("single_elec"),
            ARCH_NAMES["single_elec"],
        ),
        "hybrid": (_arch_colour("hybrid"), _arch_ls("hybrid"), ARCH_NAMES["hybrid"]),
    }
    if log_y:
        ax.set_yscale("log")
    ax.set_ylim(*ylim)
    if log_y:
        log_ticks(ax)
    ax.set_xlim(0.0, eta_hi)
    ax.axvspan(*water_range, color="#1f77b4", alpha=0.10, lw=0, zorder=0)
    ax.text(
        water_range[0] + 0.005,
        water_label_y,
        "measured water\nutilisation",
        transform=ax.get_xaxis_transform(),
        ha="left",
        va="center",
        color="#1f5f9a",
        fontsize=8,
        linespacing=1.0,
        zorder=6,
        bbox=dict(fc="white", ec="none", alpha=0.7, pad=1.0),
    )
    shade_beyond_measured(ax, eta_measured, hi=eta_hi)
    for key in ("multimode_f0", "single_chem", "single_elec", "hybrid"):
        w = refs.get(key, float("nan"))
        if math.isnan(w):
            continue
        col, ls, lab = ref_style[key]
        ax.axhline(w, color=col, ls=ls, lw=1.1, label=lab, zorder=2)
    recs = [r for r in sw["records"] if not math.isnan(r["wet"])]
    com = [r for r in recs if r["in_range"]]
    ext = [r for r in recs if not r["in_range"] and r["eta_u"] > eta_range[1]]
    if com and ext:
        ext = [com[-1]] + ext
    col = _arch_colour("multimode")
    ax.plot(
        [r["eta_u"] for r in com],
        [r["wet"] for r in com],
        "-",
        color=col,
        lw=1.8,
        zorder=4,
        label=(
            f"{ARCH_NAMES['multimode']}, ${eta_range[0]:g} \\leq "
            f"\\eta_u \\leq {eta_range[1]:g}$"
        ),
    )
    ax.plot(
        [r["eta_u"] for r in ext],
        [r["wet"] for r in ext],
        "--",
        color=col,
        lw=1.4,
        zorder=4,
        label=f"{ARCH_NAMES['multimode']}, $\\eta_u > {eta_range[1]:g}$",
    )
    ref_of = {"own_chem": refs["multimode_f0"], "single_chem": refs["single_chem"]}
    for crit, x_c in sw["min_viable"].items():
        if x_c is None or crit not in notes or notes[crit] is None:
            continue
        text, xy_text = notes[crit]
        y_c = ref_of[crit]
        ax.plot(x_c, y_c, "o", color="k", ms=4, zorder=6)
        ax.vlines(x_c, ylim[0], y_c, color="k", lw=0.6, alpha=0.5, zorder=3)
        ax.annotate(
            f"{x_c:.3f} {text}",
            xy=(x_c, y_c),
            xytext=xy_text,
            textcoords="data",
            ha="left",
            va="center",
            color=INK,
            zorder=6,
            arrowprops=dict(
                arrowstyle="-", lw=0.6, color="k", alpha=0.7, shrinkA=1, shrinkB=2
            ),
        )
    ax.set_ylabel(WET_LABEL)


def fig_res_min_eta_u(
    sw_a,
    sw_b,
    *,
    eta_range,
    eta_measured,
    water_range=(0.05, 0.20),
    eta_hi=1.0,
    ylim_a=(0.0, 30.0),
    ylim_b=(1e2, 3e3),
):
    """Bipropellant multimode wet mass against eta_u of N2O, (a) LEO and
    (b) GEO transfer, with the eta_u-independent references, the threshold
    against single-chemical and the measured water utilisation band; shaded
    above the measured xenon utilisation."""
    with print_style():
        fig, (a1, a2) = new_figure(
            FULL_W, 13.5, nrows=2, sharex=True, hspace=0.6, l=1.75, r=0.3, b=2.6, t=0.45
        )
        a1.tick_params(labelbottom=False)
        _eta_panel(
            a1,
            sw_a,
            eta_range,
            eta_measured,
            eta_hi,
            log_y=False,
            ylim=ylim_a,
            water_range=water_range,
            water_label_y=0.80,
            notes={
                "own_chem": None,
                "single_chem": ("vs. single-chemical", (0.36, 16.0)),
            },
        )
        _eta_panel(
            a2,
            sw_b,
            eta_range,
            eta_measured,
            eta_hi,
            log_y=True,
            ylim=ylim_b,
            water_range=water_range,
            water_label_y=0.30,
            notes={
                "own_chem": None,
                "single_chem": ("vs. single-chemical", (0.25, 370.0)),
            },
        )
        a2.set_xlabel(r"$\eta_u$ [-]")
        h, l = a1.get_legend_handles_labels()
        fig.legend(
            h,
            l,
            loc="lower center",
            ncol=2,
            frameon=False,
            bbox_to_anchor=(0.5, 0.0),
            columnspacing=1.4,
        )
    return fig


# LEO collision-avoidance feasibility map


def fig_res_ca_map(
    boundary,
    *,
    t_range_h=(2.0, 72.0),
    dv_range=(0.01, 1.0),
    nominal=(24.0, 0.1),
    t_ref_h=24.0,
    crossing_h=None,
    grid_pts=None,
):
    """Per-event avoidance delta-v against warning time: where the
    electrodeless unit can fly the avoidance burn within the warning time
    and where the chemical mode is needed.

    boundary : [(dv_event, t_min_h)], dv_event ascending; t_min_h is the
        shortest warning time [h] at which the electric burn fits (nan if
        none within t_range_h, -inf if it fits across all of it).
    nominal : (t_warn_h, dv_event) of the nominal point.
    crossing_h : where the budget line dv_event = nominal[1] t_ref_h / t_warn
        meets the boundary.
    grid_pts : optional [(t_h, dv_event, feasible)] spot checks.
    """
    lo, hi = t_range_h
    dv = np.array([b[0] for b in boundary], float)
    t = np.array([b[1] for b in boundary], float)
    t_plot = np.where(np.isnan(t), hi, np.where(np.isinf(t), lo, np.clip(t, lo, hi)))
    with print_style():
        fig, (ax,) = new_figure(
            FULL_W, 9.1, l=1.9, r=0.3, b=below_margin(2) + 0.2, t=0.2
        )
        ax.fill_betweenx(dv, lo, t_plot, color="#ff7f0e", alpha=0.16, lw=0, zorder=0)
        ax.fill_betweenx(dv, t_plot, hi, color="#2ca02c", alpha=0.16, lw=0, zorder=0)
        inside = ~np.isnan(t) & ~np.isinf(t) & (t >= lo) & (t <= hi)
        ax.plot(
            t[inside], dv[inside], color="k", lw=1.2, zorder=3, label="electric limit"
        )
        if grid_pts:
            for th, d, ok in grid_pts:
                ax.plot(
                    th,
                    d,
                    "o",
                    ms=2.2,
                    zorder=4,
                    color=("#2ca02c" if ok else "#ff7f0e"),
                    mec="none",
                    alpha=0.9,
                )
        tt = np.geomspace(lo, hi, 200)
        ax.plot(
            tt,
            nominal[1] * t_ref_h / tt,
            color="0.25",
            ls="--",
            lw=1.0,
            zorder=3,
            label=(
                f"budget line, {nominal[1]:g} m/s"
                f" $\\times$ ({t_ref_h:g} h / $t_{{\\mathrm{{warn}}}}$)"
            ),
        )
        ax.plot(
            *nominal,
            marker="o",
            ms=5,
            mfc="k",
            mec="k",
            ls="none",
            zorder=6,
            label="nominal point",
        )
        if crossing_h is not None:
            ax.axvline(
                crossing_h,
                color="k",
                ls="--",
                lw=0.9,
                zorder=4,
                label="minimum warning time",
            )
            ax.text(
                crossing_h / 1.03,
                dv_range[1] / 1.25,
                f"{crossing_h:.1f} h",
                ha="right",
                va="top",
                color=INK,
                zorder=6,
                bbox=dict(fc="white", ec="none", alpha=0.8, pad=1.0),
            )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(lo, hi)
        ax.set_ylim(*dv_range)
        xt = [x for x in (2, 4, 6, 12, 24, 48, 72) if lo <= x <= hi]
        ax.set_xticks(xt)
        ax.set_xticklabels([f"{x:g}" for x in xt])
        ax.xaxis.set_minor_locator(FixedLocator([]))
        log_ticks(ax)
        ax.set_xlabel(r"$t_{\mathrm{warn}}$ [h]")
        ax.set_ylabel("Per-event avoidance $\\Delta v$\nat the warning time [m/s]")
        ax.text(
            0.03,
            0.62,
            "chemical mode required",
            transform=ax.transAxes,
            ha="left",
            va="top",
            color="#a04f00",
        )
        ax.text(
            0.97,
            0.05,
            "electric burn fits\nthe warning window",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            color="#1b6e1b",
            linespacing=1.0,
        )
        h, _ = ax.get_legend_handles_labels()
        legend_below(fig, ax, h, ncol=2, columnspacing=1.6, label_cm=1.15)
    return fig


# wet-mass bars, one mission

# Row order of both missions, top to bottom
WET_ORDER = ("multimode", "mono_multimode", "single_elec", "hybrid", "single_chem")


def fig_res_wet_bars(rows, *, truncate_above=None, h_cm=6.2):
    """Horizontal stacked bars, dry mass then loaded propellant, one per
    architecture in WET_ORDER from the top; values on the segments, or in
    brackets after the wet mass where a segment is too narrow.

    rows : [{key, dry, prop, wet, converged}, ...].
    truncate_above : x-axis limit [kg]; a longer bar is cut at the limit
                     and labelled with its full value.
    """
    by_key = {r["key"]: r for r in rows if r["converged"]}
    rows = [by_key[k] for k in WET_ORDER if k in by_key]
    xmax = (
        truncate_above
        if truncate_above is not None
        else max(r["wet"] for r in rows) * 1.22
    )
    fmt = "{:.2f}" if xmax < 100 else "{:.1f}"
    with print_style():
        fig, (ax,) = new_figure(FULL_W, h_cm, l=3.6, r=0.3, b=1.0, t=0.6)
        y = list(range(len(rows)))[::-1]
        for yi, r in zip(y, rows):
            col = _arch_colour(r["key"])
            cut = r["wet"] > xmax
            dry = min(r["dry"], xmax)
            prop = min(r["prop"], xmax - dry)
            ax.barh(yi, dry, color=col, height=0.6, zorder=3)
            ax.barh(
                yi,
                prop,
                left=dry,
                color=col,
                alpha=0.40,
                height=0.6,
                edgecolor="white",
                lw=1.0,
                zorder=3,
            )
            if cut:
                # break mark at the cut
                for dx in (-0.018, -0.008):
                    ax.plot(
                        [
                            xmax * (1 + dx) - 0.004 * xmax,
                            xmax * (1 + dx) + 0.004 * xmax,
                        ],
                        [yi - 0.38, yi + 0.38],
                        color="white",
                        lw=2.2,
                        zorder=4,
                        clip_on=False,
                    )
                ax.text(
                    xmax * 0.97,
                    yi,
                    f"{r['wet']:.1f} kg",
                    ha="right",
                    va="center",
                    color=INK,
                    zorder=5,
                    bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.0),
                )
                continue
            t_wet = ax.text(
                r["wet"] + 0.01 * xmax,
                yi,
                f"{fmt.format(r['wet'])} kg",
                ha="left",
                va="center",
                color=INK,
            )
            if dry / xmax > 0.07 and prop / xmax > 0.07:
                ax.text(
                    dry / 2,
                    yi,
                    fmt.format(r["dry"]),
                    ha="center",
                    va="center",
                    color="white",
                    fontsize=8,
                    zorder=5,
                )
                ax.text(
                    dry + prop / 2,
                    yi,
                    fmt.format(r["prop"]),
                    ha="center",
                    va="center",
                    color=INK,
                    fontsize=8,
                    zorder=5,
                )
            else:
                ax.annotate(
                    f" ({fmt.format(r['dry'])} dry / {fmt.format(r['prop'])} prop.)",
                    xy=(1.0, 0.5),
                    xycoords=t_wet,
                    ha="left",
                    va="center",
                    color=INK,
                    fontsize=7,
                    zorder=5,
                )
        ax.set_yticks(y)
        ax.set_yticklabels([ARCH_NAMES[r["key"]] for r in rows])
        ax.set_xlim(0, xmax)
        ax.set_ylim(-0.6, len(rows) - 0.4)
        ax.set_xlabel(WET_LABEL)
        ax.grid(axis="y", visible=False)
        ax.tick_params(axis="y", length=0)
        ax.legend(
            handles=[
                Patch(color="0.35", label="dry mass"),
                Patch(color="0.35", alpha=0.4, label="loaded propellant"),
            ],
            loc="lower right",
            bbox_to_anchor=(1.0, 1.0),
            ncol=2,
            frameon=False,
            borderaxespad=0.1,
        )
    return fig


# LEO response variant

_RESP_ORDER = ("multimode", "mono_multimode", "single_elec", "hybrid", "single_chem")


def fig_res_response(rows, *, se_limit_total, x_max=100.0, h_cm=9.0):
    """Propulsion wet mass against the total response delta-v of the LEO
    response variant, five architectures, each drawn where feasible.

    rows : [{key, total, wet, feasible}, ...].
    se_limit_total : total response [m/s] at which single-electric is
                     first excluded, drawn as a vertical line.
    Points with total > x_max are not drawn.
    """
    with print_style():
        fig, (ax,) = new_figure(
            FULL_W, h_cm, l=1.45, r=0.3, b=below_margin(4) + 0.15, t=0.2
        )
        curves = {}
        for key in _RESP_ORDER:
            pts = sorted(
                (r["total"], r["wet"])
                for r in rows
                if r["key"] == key and r["feasible"] and r["total"] <= x_max + 1e-9
            )
            if not pts:
                continue
            curves[key] = pts
            ax.plot(
                *zip(*pts),
                color=_arch_colour(key),
                ls=_arch_ls(key),
                lw=1.2,
                marker="o",
                ms=2.6,
                label=ARCH_NAMES[key],
                zorder=3,
            )
        handles, _ = ax.get_legend_handles_labels()
        ax.axvline(se_limit_total, color="k", ls="--", lw=0.9, zorder=2)
        handles.append(
            Line2D(
                [],
                [],
                color="k",
                ls="--",
                lw=0.9,
                label=f"single-electric limit, {se_limit_total:.1f} m/s",
            )
        )
        ax.set_xlim(0.0, x_max)
        # y axis from 1 kg below the lightest design
        top = max(w for pts in curves.values() for _, w in pts)
        low = min(w for pts in curves.values() for _, w in pts)
        ax.set_ylim(math.floor(low) - 1.0, math.ceil(top / 5.0) * 5.0)
        ax.set_xlabel(r"Total response $\Delta v^{\mathrm{imp}}_{\mathrm{resp}}$ [m/s]")
        ax.set_ylabel(WET_LABEL)
        legend_below(fig, ax, handles, ncol=2, columnspacing=1.6, label_cm=1.15)
    return fig


# Chapter 3: Sheppard and Little (2022) exit-velocity check

SL_MARKERS = {"Ar": "o", "H2O": "s", "Xe": "D"}


def fig_sl_validation(
    points,
    per_propellant,
    *,
    order=("Ar", "H2O", "Xe"),
    ylim=(0.0, 1.8),
    xlim=(5.0, 15.0),
):
    """Per-point exit-velocity factor k = u_e,meas / u_e,model against the
    interpolated T_e, one marker per point, each propellant's mean k as a
    horizontal line in its colour and linestyle (PROPELLANT_STYLE; Xe green).
    `points` maps propellant -> (Te array, k array); `per_propellant` maps
    propellant -> {'n', 'mean', 'se'} (computed by the caller)."""
    with print_style():
        n_leg = len(order)
        fig, (ax,) = new_figure(
            FULL_W, 5.6 + below_margin(n_leg) - 0.95, b=below_margin(n_leg), l=1.45
        )
        handles = []
        for p in order:
            st = PROPELLANT_STYLE[p]
            te, k = points[p]
            s = per_propellant[p]
            ax.plot(
                te,
                k,
                ls="none",
                marker=SL_MARKERS[p],
                ms=5.0,
                mfc=st["color"],
                mec="white",
                mew=0.6,
                zorder=3,
            )
            ax.axhline(
                s["mean"], color=st["color"], lw=1.1, ls=st["linestyle"], zorder=2
            )
            handles.append(
                Line2D(
                    [],
                    [],
                    color=st["color"],
                    ls=st["linestyle"],
                    lw=1.1,
                    marker=SL_MARKERS[p],
                    ms=5.0,
                    mfc=st["color"],
                    mec="white",
                    mew=0.6,
                    label=(
                        f"{st['label']}: {s['n']} points, mean $k$ = "
                        f"{s['mean']:.3f} (standard error {s['se']:.3f})"
                    ),
                )
            )
        ax.set_xlabel(r"Electron temperature at the exit-velocity point, $T_e$ [eV]")
        ax.set_ylabel(r"$k = u_{e,\mathrm{meas}} / u_{e,\mathrm{model}}$")
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        legend_below(fig, ax, handles, ncol=1)
    return fig
