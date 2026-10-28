"""Convert the digitised Song et al. (2019) Fig. 28 curves to the loader CSVs."""

# Reads the WebPlotDigitizer exports in
# data/cross_sections/n2o/digitised_raw/Song2019_fig28/ and writes, next to
# the other N2O tables,
#
#     electronic_swarm_Qex1_5p6eV.csv     swarm component, threshold 5.6 eV
#     electronic_swarm_Qex2_8p5eV.csv     swarm component, threshold 8.5 eV
#     electronic_swarm_Qex3_9p6eV.csv     swarm component, threshold 9.6 eV
#     electronic_swarm_total_fig28.csv    their sum (thick curve), check only
#
# columns E_eV,sigma_1e16_cm2, with a '#' provenance header.
#
# Processing, per curve (the header of each output repeats what was applied):
#     1. read two comma-separated columns (x, y), no header, '.' decimal;
#     2. stable sort by x;
#     3. drop exact duplicate (x, y) rows;
#     4. repeated x with different y (vertical runs where the curve rises
#        steeply near threshold on the log energy axis, one digitiser pixel
#        column; x equal to within X_RTOL relative): keep the point with the
#        largest y, so the energy grid is strictly increasing as the log-log
#        interpolant needs;
#     5. drop points below the curve's threshold and add an explicit
#        (threshold, 0) row;
#     6. flag, without correcting: non-positive y, and every local extremum of
#        y after processing other than the curve's main maximum.
#
# Usage:
#     python scripts/convert_song_fig28.py            # convert, print the log
#     python scripts/convert_song_fig28.py --check    # also run the consistency
#                                                     # checks (components vs
#                                                     # digitised total vs
#                                                     # Song Table XIII sum)

import argparse
import hashlib
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
N2O_DIR = REPO / "data" / "cross_sections" / "n2o"
RAW_DIR = N2O_DIR / "digitised_raw" / "Song2019_fig28"
X_RTOL = 1e-9  # x values closer than this (relative) are one pixel column

SOURCE = (
    "Song et al. (2019), J. Phys. Chem. Ref. Data 48, 043104, Fig. 28 "
    "(electronic excitation cross sections of N2O)"
)

# raw file, output file, curve, threshold [eV], description
CURVES = (
    (
        "Qex1.csv",
        "electronic_swarm_Qex1_5p6eV.csv",
        "Qex1",
        5.6,
        "thin curve, component after Nakamura (2007)",
    ),
    (
        "Qex2.csv",
        "electronic_swarm_Qex2_8p5eV.csv",
        "Qex2",
        8.5,
        "thin curve, component after Nakamura (2007)",
    ),
    (
        "Qex3.csv",
        "electronic_swarm_Qex3_9p6eV.csv",
        "Qex3",
        9.6,
        "thin curve, component after Nakamura (2007)",
    ),
    (
        "Total.csv",
        "electronic_swarm_total_fig28.csv",
        "total",
        5.6,
        "thick curve, sum of the three components",
    ),
)


def read_raw(path):
    """Two comma-separated float columns, no header, '.' decimal."""
    rows = []
    for n, line in enumerate(path.read_text(encoding="ascii").splitlines(), 1):
        if not line.strip():
            continue
        parts = line.split(",")
        if len(parts) != 2:
            raise ValueError(f"{path.name} line {n}: expected 2 fields, got {line!r}")
        rows.append((float(parts[0]), float(parts[1])))
    return np.array(rows, dtype=float)


def process(raw, threshold):
    """Steps 2-6 of the module docstring. Returns (E, sigma, log lines)."""
    log = [
        f"raw rows: {len(raw)}; raw x order monotone: "
        f"{bool(np.all(np.diff(raw[:, 0]) >= 0))}"
    ]
    order = np.argsort(raw[:, 0], kind="stable")
    xy = raw[order]

    keep = np.ones(len(xy), dtype=bool)
    keep[1:] = ~((xy[1:, 0] == xy[:-1, 0]) & (xy[1:, 1] == xy[:-1, 1]))
    log.append(f"exact duplicate (x, y) rows dropped: {int((~keep).sum())}")
    xy = xy[keep]

    # one group per digitiser pixel column: x equal to within X_RTOL
    gid = np.concatenate(([0], np.cumsum(np.diff(xy[:, 0]) > X_RTOL * xy[1:, 0])))
    ux, uy, counts, glog = [], [], [], []
    for g in range(gid[-1] + 1):
        grp = xy[gid == g]
        k = int(np.argmax(grp[:, 1]))  # first row with the largest y
        ux.append(grp[k, 0])
        uy.append(grp[k, 1])
        counts.append(len(grp))
        if len(grp) > 1:
            xs = np.unique(grp[:, 0])
            xtxt = (
                repr(float(xs[0]))
                if len(xs) == 1
                else f"{float(xs[0])!r} to {float(xs[-1])!r}"
            )
            glog.append(
                f"    x = {xtxt}: {len(grp)} points, y {grp[:, 1].min():.6g} to "
                f"{grp[:, 1].max():.6g}, kept ({float(grp[k, 0])!r}, "
                f"{float(grp[k, 1])!r})"
            )
    ux, uy, counts = np.array(ux), np.array(uy), np.array(counts)
    rep = counts > 1
    log.append(
        f"repeated-x groups (x equal to within {X_RTOL:g} relative): {int(rep.sum())} "
        f"({int(counts[rep].sum())} points); kept the largest y in each, dropped "
        f"{int(counts[rep].sum() - rep.sum())}"
    )
    log += glog

    nonpos = uy <= 0
    log.append(f"non-positive y (flagged, kept): {int(nonpos.sum())}")

    below = ux < threshold
    log.append(f"points below threshold {threshold} eV dropped: {int(below.sum())}")
    for x, y in zip(ux[below], uy[below]):
        log.append(f"    ({float(x)!r}, {float(y)!r})")
    if np.any(ux == threshold):
        raise ValueError(f"a digitised point lies exactly at the threshold {threshold}")
    E = np.concatenate(([threshold], ux[~below]))
    s = np.concatenate(([0.0], uy[~below]))

    d = np.sign(np.diff(s[1:]))
    turns = [
        i + 1
        for i in range(len(d) - 1)
        if d[i] != 0 and d[i + 1] != 0 and d[i] != d[i + 1]
    ]
    i_max = int(np.argmax(s))
    extra = [i + 1 for i in turns if i + 1 != i_max]
    log.append(
        f"main maximum: {float(s[i_max])!r} at {float(E[i_max])!r} eV; "
        f"other local extrema (flagged, kept): {len(extra)}"
    )
    for i in extra:
        kind = "max" if s[i] > s[i - 1] else "min"
        log.append(
            f"    local {kind} at E = {float(E[i])!r} eV, sigma = {float(s[i])!r} "
            f"(neighbours {float(s[i - 1])!r}, {float(s[i + 1])!r})"
        )
    log.append(
        f"output rows: {len(E)} (incl. threshold row); E {float(E[1])!r} to "
        f"{float(E[-1])!r} eV"
    )
    return E, s, log


def header(raw_name, curve, threshold, desc):
    """Citation and units header of an output CSV."""
    return (
        "\n".join(
            (
                f"# N2O electronic excitation, swarm-derived {curve}: {desc}.",
                f"# Source: {SOURCE}, digitised",
                f"#   (digitised_raw/Song2019_fig28/{raw_name}); threshold {threshold} eV.",
                "# Units: E_eV in eV; sigma in 1e-16 cm^2.",
                "E_eV,sigma_1e16_cm2",
            )
        )
        + "\n"
    )


def convert():
    out = {}
    for raw_name, out_name, curve, thr, desc in CURVES:
        path = RAW_DIR / raw_name
        sha = hashlib.sha256(path.read_bytes()).hexdigest()
        E, s, log = process(read_raw(path), thr)
        body = "".join(f"{float(e)!r},{float(v)!r}\n" for e, v in zip(E, s))
        text = header(raw_name, curve, thr, desc) + body
        (N2O_DIR / out_name).write_bytes(text.encode("ascii"))
        print(f"== {raw_name} -> {out_name} (sha256 {sha[:12]}...)")
        for line in log:
            print("  " + line)
        out[curve] = (E, s, thr)
    return out


# Consistency checks (report only; nothing is corrected)


def loglog(E, s, x):
    """Log-log interpolation over the positive points, as the model's
    interpolant (cross_sections._build_loglog_interpolant): 0 below the
    first positive point; NaN above the last one (outside digitised range)."""
    pos = s > 0
    e, v = E[pos], s[pos]
    y = np.zeros_like(x)
    inside = (x >= e[0]) & (x <= e[-1])
    y[inside] = np.exp(np.interp(np.log(x[inside]), np.log(e), np.log(v)))
    y[x > e[-1]] = np.nan
    return y


def _stats(name, x, num, ref, lines, flag=0.15):
    with np.errstate(divide="ignore", invalid="ignore"):
        dev = (num - ref) / ref
    ok = np.isfinite(dev)
    lines.append(
        f"{name}: n = {int(ok.sum())} energies with both values defined "
        f"({int((~ok).sum())} excluded: outside a digitised range or 0/0)"
    )
    if ok.any():
        a = np.abs(dev[ok])
        k = int(np.argmax(a))
        lines.append(
            f"    max |dev| {100 * a[k]:.2f} % at {x[ok][k]:.4g} eV; "
            f"median |dev| {100 * np.median(a):.2f} %; "
            f"> {100 * flag:.0f} %: {int((a > flag).sum())} energies"
        )
    return dev


def check(curves):
    lines = []
    comp = [curves[c] for c in ("Qex1", "Qex2", "Qex3")]
    Et, st, _ = curves["total"]

    # (a) sum of components vs digitised total, log grid 5.6-1000 eV
    x = np.geomspace(5.6, 1000.0, 100)
    parts = [loglog(E, s, x) for E, s, _ in comp]
    ssum = sum(parts)
    tot = loglog(Et, st, x)
    lines.append(
        "(a) Qex1 + Qex2 + Qex3 vs digitised total, 100 log-spaced energies, 5.6-1000 eV"
    )
    dev = _stats("    sum vs total", x, ssum, tot, lines)
    lines.append("    E_eV, Qex1, Qex2, Qex3, sum, total, dev_pct  [1e-16 cm^2]")
    for i in range(len(x)):
        d = f"{100 * dev[i]:+.2f}" if np.isfinite(dev[i]) else "n/a"
        fl = "  FLAG >15 %" if np.isfinite(dev[i]) and abs(dev[i]) > 0.15 else ""
        lines.append(
            f"    {x[i]:9.4f}, {parts[0][i]:.5g}, {parts[1][i]:.5g}, {parts[2][i]:.5g}, "
            f"{ssum[i]:.5g}, {tot[i]:.5g}, {d}{fl}"
        )

    # (b) both vs Song Table XIII swarm sum at its tabulated energies
    rows = [
        ln
        for ln in (N2O_DIR / "electronic_swarm_sum.csv").read_text().splitlines()
        if ln.strip() and not ln.startswith("#")
    ]
    assert rows[0] == "E_eV,swarm_sum", rows[0]
    tab = np.array([[float(v) for v in ln.split(",")] for ln in rows[1:]])
    xe, ref = tab[:, 0], tab[:, 1] * 0.01  # 1e-18 -> 1e-16 cm^2
    psum = sum(loglog(E, s, xe) for E, s, _ in comp)
    ptot = loglog(Et, st, xe)
    lines.append("")
    lines.append(
        "(b) vs Song Table XIII swarm-derived sum (electronic_swarm_sum.csv, "
        "1e-18 cm^2 -> 1e-16), at its tabulated energies"
    )
    d1 = _stats("    component sum vs Table XIII", xe, psum, ref, lines)
    d2 = _stats("    digitised total vs Table XIII", xe, ptot, ref, lines)
    lines.append("    E_eV, TableXIII, comp_sum, dev_pct, total, dev_pct  [1e-16 cm^2]")

    def f(d):
        if not np.isfinite(d):
            return "n/a"
        return f"{100 * d:+.2f}" + (" FLAG" if abs(d) > 0.15 else "")

    for i in range(len(xe)):
        lines.append(
            f"    {xe[i]:7.2f}, {ref[i]:.5g}, {psum[i]:.5g}, {f(d1[i])}, "
            f"{ptot[i]:.5g}, {f(d2[i])}"
        )
    print("\n".join(lines))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()
    curves = convert()
    if args.check:
        print()
        check(curves)
