"""
Page geometry and markers shared by the thesis figures: vector PDF at the
print width (FULL_W 16 cm, HALF_W 7.8 cm), margins in centimetres, 9 pt text.
"""

from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from matplotlib.lines import Line2D
from matplotlib.ticker import LogFormatterMathtext, LogLocator, NullFormatter


CM = 1.0 / 2.54
FULL_W = 16.0  # cm
HALF_W = 7.8  # cm

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
    "xtick.major.size": 3.0,
    "ytick.major.size": 3.0,
    "xtick.minor.size": 1.8,
    "ytick.minor.size": 1.8,
    "axes.grid": True,
    "grid.linewidth": 0.4,
    "grid.alpha": 0.35,
    "grid.linestyle": "-",
    "savefig.bbox": None,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.dpi": 300,
    "legend.framealpha": 0.92,
    "legend.edgecolor": "0.8",
    "legend.borderpad": 0.3,
    "legend.handlelength": 1.6,
    "legend.labelspacing": 0.3,
    "legend.borderaxespad": 0.4,
    "axes.titlepad": 3.0,
}

INK = "0.15"  # annotation text
MARK = "0.40"  # range-marker text and lines (unobtrusive)
RULER = "0.72"  # representative-range ruler bars

# Representative electron-temperature range of each thruster class [eV]
REPRESENTATIVE_TE = (("helicon", 1.0, 10.0), ("ECR", 5.0, 23.0))


@contextmanager
def print_style():
    """Context in which every thesis figure is built."""
    with plt.rc_context(PRINT_RC):
        yield


def new_figure(
    w_cm,
    h_cm,
    *,
    ncols=1,
    nrows=1,
    l=1.35,
    r=0.3,
    b=0.95,
    t=0.25,
    wspace=1.4,
    hspace=0.9,
    width_ratios=None,
    height_ratios=None,
    sharex=False,
    sharey=False,
):
    """Figure of exact size with margins and gaps fixed in centimetres.

    `wspace` / `hspace` are the gaps between panels in cm. Returns the
    figure and a flat list of axes in row-major order.
    """
    fig = plt.figure(figsize=(w_cm * CM, h_cm * CM))
    inner_w = w_cm - l - r
    inner_h = h_cm - b - t
    gs = fig.add_gridspec(
        nrows,
        ncols,
        left=l / w_cm,
        right=1.0 - r / w_cm,
        bottom=b / h_cm,
        top=1.0 - t / h_cm,
        wspace=(
            wspace * ncols / max(inner_w - wspace * (ncols - 1), 1e-6)
            if ncols > 1
            else 0.0
        ),
        hspace=(
            hspace * nrows / max(inner_h - hspace * (nrows - 1), 1e-6)
            if nrows > 1
            else 0.0
        ),
        width_ratios=width_ratios,
        height_ratios=height_ratios,
    )
    axes = []
    for i in range(nrows):
        for j in range(ncols):
            kw = {}
            if sharex and axes:
                kw["sharex"] = axes[0]
            if sharey and axes:
                kw["sharey"] = axes[0]
            axes.append(fig.add_subplot(gs[i, j], **kw))
    return fig, axes


def _offset(ax, dy_pt):
    """Blended transform: x in data, y in axes fraction, shifted dy points."""
    return mtransforms.blended_transform_factory(
        ax.transData, ax.transAxes
    ) + mtransforms.ScaledTranslation(0.0, dy_pt / 72.0, ax.figure.dpi_scale_trans)


def mark_representative_range(ax, *, labels=True, rows_pt=(4.0, 12.5)):
    """Ruler above the frame, one bar per thruster class, clipped to the
    x limits; needs about 0.55 cm of top margin."""
    lo_ax, hi_ax = ax.get_xlim()
    span = hi_ax - lo_ax
    for (name, lo, hi), dy in zip(REPRESENTATIVE_TE, rows_pt):
        a, b = max(lo, lo_ax), min(hi, hi_ax)
        if b <= a:
            continue
        tr = _offset(ax, dy)
        ax.plot(
            [a, b],
            [1.0, 1.0],
            transform=tr,
            color=RULER,
            lw=2.6,
            solid_capstyle="butt",
            clip_on=False,
            zorder=1,
        )
        if not labels:
            continue
        if (hi_ax - b) / span > 0.16:
            ax.text(
                b + 0.012 * span,
                1.0,
                name,
                transform=tr,
                ha="left",
                va="center",
                color=MARK,
                fontsize=8,
            )
        else:
            ax.text(
                a - 0.012 * span,
                1.0,
                name,
                transform=tr,
                ha="right",
                va="center",
                color=MARK,
                fontsize=8,
            )


TE_LINE_STYLES = {"boundary": (0, (4.0, 2.5)), "reference": (0, (1.0, 1.6))}


def te_line_labels(te_boundary, te_ref):
    """Label text of the two T_e lines, each carrying its value."""
    return (
        f"electropositive boundary, {te_boundary:.1f} eV"
        if te_boundary is not None
        else None,
        f"reference, {te_ref:g} eV" if te_ref is not None else None,
    )


def te_line_handles(te_boundary, te_ref):
    """Legend entries for the two T_e lines, labelled with their values."""
    lab_b, lab_r = te_line_labels(te_boundary, te_ref)
    out = []
    for lab, key in ((lab_b, "boundary"), (lab_r, "reference")):
        if lab is not None:
            out.append(
                Line2D([], [], color=MARK, lw=0.8, ls=TE_LINE_STYLES[key], label=lab)
            )
    return out


def mark_te_lines(ax, *, te_boundary=None, te_ref=None):
    """Electropositive boundary (dashed) and reference T_e (dotted), where
    inside the x limits."""
    lo, hi = ax.get_xlim()
    for x, key in ((te_boundary, "boundary"), (te_ref, "reference")):
        if x is None or not (lo <= x <= hi):
            continue
        ax.axvline(x, color=MARK, lw=0.8, ls=TE_LINE_STYLES[key], zorder=1.5)


def mark_te(ax, *, te_boundary, te_ref, labels=True, ruler=True):
    """Ruler and T_e lines; call after the x limits are final."""
    if ruler:
        mark_representative_range(ax, labels=labels)
    mark_te_lines(ax, te_boundary=te_boundary, te_ref=te_ref)


def mark_te_stacked(axes, *, te_boundary, te_ref):
    """T_e markers on a stacked figure: the ruler on the top panel, the
    lines on every panel."""
    mark_te(axes[0], te_boundary=te_boundary, te_ref=te_ref)
    for ax in axes[1:]:
        mark_te_lines(ax, te_boundary=te_boundary, te_ref=te_ref)


def stacked_figure(
    h_cm, *, height_ratios=None, hspace=0.5, l=1.75, r=0.3, b=1.05, t=0.62
):
    """Two panels stacked at full width sharing x. Returns (fig, axes)."""
    fig, axes = new_figure(
        FULL_W,
        h_cm,
        nrows=2,
        sharex=True,
        height_ratios=height_ratios,
        hspace=hspace,
        l=l,
        r=r,
        b=b,
        t=t,
    )
    axes[0].tick_params(labelbottom=False)
    return fig, axes


LEGEND_ROW_CM = 0.42  # height of one 9 pt legend row


def below_margin(n_rows, b=1.05):
    """Bottom margin [cm] for the x label and a legend of n_rows below it."""
    return b + 0.15 + LEGEND_ROW_CM * n_rows


def legend_below(fig, ax, handles, *, ncol=1, label_cm=0.95, **kw):
    """Frameless legend under the x label, left-aligned with the frame;
    `label_cm` is the depth of the tick labels and x label [cm]."""
    w_cm = fig.get_figwidth() / CM
    h_cm = fig.get_figheight() / CM
    pos = ax.get_position()
    return fig.legend(
        handles=handles,
        loc="upper left",
        ncol=ncol,
        bbox_to_anchor=(pos.x0 - 0.9 / w_cm, pos.y0 - label_cm / h_cm),
        frameon=False,
        borderaxespad=0.0,
        **kw,
    )


def log_ticks(ax, which="y", *, minor=True):
    """Decade ticks labelled in 10^n form, unlabelled minor ticks."""
    axis = ax.yaxis if which == "y" else ax.xaxis
    axis.set_major_locator(LogLocator(base=10.0))
    axis.set_major_formatter(LogFormatterMathtext(base=10.0))
    axis.set_minor_locator(
        LogLocator(base=10.0, subs=tuple(range(2, 10)))
        if minor
        else LogLocator(base=10.0, subs=())
    )
    axis.set_minor_formatter(NullFormatter())


def shade_beyond_measured(
    ax, eta_measured, *, hi=1.0, label=True, label_y=0.97, va="top"
):
    """Shade eta_u above the measured xenon utilisation, up to `hi`."""
    ax.axvspan(eta_measured, hi, color="0.5", alpha=0.10, lw=0, zorder=0)
    if label:
        ax.text(
            0.5 * (eta_measured + hi),
            label_y,
            "beyond measured\nxenon utilisation",
            transform=ax.get_xaxis_transform(),
            ha="center",
            va=va,
            color="0.30",
            fontsize=8,
            linespacing=1.0,
            zorder=6,
        )


def _panel_groups(fig):
    """Axes of `fig` grouped into panels (twins, insets and colorbars join
    their host), in reading order."""
    axes = [a for a in fig.axes if a.get_visible()]
    hosts, extra = [], []
    for a in axes:
        if getattr(a, "_colorbar", None) is not None:
            extra.append(a)
            continue
        p = a.get_position()
        for g in hosts:
            q = g[0].get_position()
            same = all(abs(u - v) < 1e-6 for u, v in zip(p.bounds, q.bounds))
            inside = (
                p.x0 >= q.x0 - 1e-6
                and p.x1 <= q.x1 + 1e-6
                and p.y0 >= q.y0 - 1e-6
                and p.y1 <= q.y1 + 1e-6
            )
            if same or inside:
                g.append(a)
                break
        else:
            hosts.append([a])
    for a in extra:  # colorbars: to the nearest panel
        c = a.get_position()
        g = min(
            hosts,
            key=lambda g: (
                (g[0].get_position().x1 - c.x0) ** 2
                + (
                    0.5 * (g[0].get_position().y0 + g[0].get_position().y1)
                    - 0.5 * (c.y0 + c.y1)
                )
                ** 2
            ),
        )
        g.append(a)
    hosts.sort(
        key=lambda g: (-round(g[0].get_position().y1, 3), g[0].get_position().x0)
    )
    return hosts


def export_panels(fig, path):
    """Write each panel of a multi-panel figure as `<stem>_a.pdf`,
    `<stem>_b.pdf`, ..., all the same page size with the frames aligned,
    and a figure legend as `<stem>_legend.pdf`. Returns the paths."""
    groups = _panel_groups(fig)
    if len(groups) < 2:
        return []
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    inch = fig.dpi_scale_trans.inverted()
    tight = [
        mtransforms.Bbox.union([a.get_tightbbox(r) for a in g]).transformed(inch)
        for g in groups
    ]
    frame = [g[0].get_window_extent(r).transformed(inch) for g in groups]
    left = max(f.x0 - t.x0 for f, t in zip(frame, tight))
    right = max(t.x1 - f.x0 for f, t in zip(frame, tight))
    below = max(f.y0 - t.y0 for f, t in zip(frame, tight))
    above = max(t.y1 - f.y0 for f, t in zip(frame, tight))
    pad = 0.03
    out = []
    path = Path(path)
    legends = list(fig.legends)
    for k, f in enumerate(frame):
        box = mtransforms.Bbox.from_extents(
            f.x0 - left - pad,
            f.y0 - below - pad,
            f.x0 + right + pad,
            f.y0 + above + pad,
        )
        # hide the other panels and figure legends while writing this one
        others = [a for j, g in enumerate(groups) if j != k for a in g]
        for a in others + legends:
            a.set_visible(False)
        p = path.with_name(f"{path.stem}_{chr(ord('a') + k)}.pdf")
        fig.savefig(p, bbox_inches=box)
        for a in others + legends:
            a.set_visible(True)
        out.append(p)
    for k, leg in enumerate(legends):
        bb = leg.get_window_extent(r).transformed(inch).padded(pad)
        hide = [a for g in groups for a in g] + [x for x in legends if x is not leg]
        for a in hide:
            a.set_visible(False)
        p = path.with_name(f"{path.stem}_legend{'' if k == 0 else k + 1}.pdf")
        fig.savefig(p, bbox_inches=bb)
        for a in hide:
            a.set_visible(True)
        out.append(p)
    return out


def save_pdf(fig, out_dir, name, *, panels=True):
    """Write `name` (a .pdf) into `out_dir` at the exact figure size, and
    each panel of a multi-panel figure."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / name
    # exact page size: no tight bounding box
    with plt.rc_context({"savefig.bbox": None, "pdf.fonttype": 42}):
        fig.savefig(path)
        extra = export_panels(fig, path) if panels else []
    plt.close(fig)
    print(f"Saved: {path}" + (f" (+ {len(extra)} panel files)" if extra else ""))
    return path
