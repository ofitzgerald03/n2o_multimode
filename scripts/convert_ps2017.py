"""Convert the digitised Petro and Sedwick (2017) curves to the loader CSVs."""

# Reads the WebPlotDigitizer exports in data/ps2017/digitised_raw/ (one
# subfolder per figure) and writes data/ps2017/<name>.csv, one file per
# curve, and data/ps2017/PROVENANCE.md.
#
# Processing, per curve:
#     1. read two comma-separated columns (x, y), no header;
#     2. stable sort by x;
#     3. drop exact duplicate (x, y) rows;
#     4. repeated x (equal to within X_RTOL relative, one digitiser pixel
#        column): one point at the mean y, taken in axis space (arithmetic on
#        a linear y axis, geometric on a log y axis);
#     4a. drop the stray clicks listed in REMOVED (matched exactly on the raw
#        floats);
#     5. negative y on a linear y axis: clamped to 0 and flagged; non-positive
#        y on a log y axis: dropped and flagged;
#     6. flag, without correcting, local extrema other than the global maximum
#        and minimum that reverse by more than EXTREMUM_K x the RMS residual,
#        and leave-one-out residuals above OUTLIER_K x RMS;
#     7. digitising error: RMS of the leave-one-out residuals about a local
#        quadratic fitted to up to N_SIDE neighbours each side (at least
#        N_SIDE_MIN), in axis space (decades on a log y axis).
#
# Usage: python scripts/convert_ps2017.py

import hashlib
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
PS_DIR = REPO / "data" / "ps2017"
RAW_DIR = PS_DIR / "digitised_raw"
X_RTOL = 1e-9  # x values closer than this (relative) are one pixel column
N_SIDE = 3  # neighbours per side in each leave-one-out quadratic fit
N_SIDE_MIN = 2  # points with fewer neighbours on a side are not evaluated
OUTLIER_K = 4.0  # residual / RMS above which a point is flagged
EXTREMUM_K = 3.0  # reversal / RMS above which a local extremum is flagged

SOURCE = (
    'Petro, E. M. and Sedwick, R. J. (2017), "Effects of Water-Vapor '
    'Propellant on Electrodeless Thruster Performance", J. Propul. Power '
    "33(6), 1410-1417, doi:10.2514/1.B36389"
)

# Per figure: raw subfolder, caption, x and y (column, symbol,
# unit) and the axis scales of the published figure
FIGS = {
    "fig1": dict(
        raw_dir="fig 1",
        caption="rate constants for ions generated in a water-vapour plasma",
        x=("Te_eV", "T_e", "eV"),
        y=("k_m3_s", "rate constant kappa", "m^3/s"),
        xscale="linear",
        yscale="log",
    ),
    "fig2": dict(
        raw_dir="fig2",
        caption="steady-state ratio of negative ions to neutral particles, Eq. (3)",
        x=("Te_eV", "T_e", "eV"),
        y=("negion_ratio", "n_i-/n_n", "-"),
        xscale="linear",
        yscale="linear",
    ),
    "fig3": dict(
        raw_dir="fig3",
        caption=(
            "distribution of positive ions in a H2O plasma; H2O+ from the "
            "main axes (0-1), H+, O+, OH+ from the inset (0-0.02)"
        ),
        x=("Te_eV", "T_e", "eV"),
        y=("ion_fraction", "positive-ion fraction", "-"),
        xscale="linear",
        yscale="linear",
    ),
    "fig5": dict(
        raw_dir="fig 5",
        caption=(
            "fractional distribution of energy loss by process, chi_j; "
            "electronic excitation includes radiating and nonradiating "
            "dissociation"
        ),
        x=("Te_eV", "T_e", "eV"),
        y=("chi", "fractional energy loss chi_j", "-"),
        xscale="linear",
        yscale="linear",
    ),
    "fig7a": dict(
        raw_dir="fig7a",
        caption="thrust efficiency vs electron temperature",
        x=("Te_eV", "T_e", "eV"),
        y=("eta_T", "thrust efficiency eta_T", "-"),
        xscale="linear",
        yscale="linear",
    ),
    "fig7b": dict(
        raw_dir="fig7b",
        caption="thrust efficiency vs upstream magnetic mirror ratio",
        x=("R_m", "R_m", "-"),
        y=("eta_T", "thrust efficiency eta_T", "-"),
        xscale="log",
        yscale="linear",
    ),
}

# Step 4a: (figure, raw file) -> ((raw x, raw y, reason), ...)
REMOVED = {
    ("fig5", "rotation.csv"): (
        (
            4.286632014291518,
            0.05737715325280757,
            "stray click (raw row 35; 0.013 below both neighbours)",
        ),
    ),
}

# figure, raw file, output name, curve

CURVES = (
    ("fig1", "H2Oplus.csv", "PS2017_fig1_k_H2Oplus", "H2O+ (ionisation)"),
    ("fig1", "H2plus.csv", "PS2017_fig1_k_H2plus", "H2+ (dissociative ionisation)"),
    ("fig1", "Hplus.csv", "PS2017_fig1_k_Hplus", "H+ (dissociative ionisation)"),
    ("fig1", "OHplus.csv", "PS2017_fig1_k_OHplus", "OH+ (dissociative ionisation)"),
    ("fig1", "Oplus.csv", "PS2017_fig1_k_Oplus", "O+ (dissociative ionisation)"),
    (
        "fig1",
        "Oplusplus.csv",
        "PS2017_fig1_k_Oplusplus",
        "O++ (dissociative ionisation)",
    ),
    ("fig1", "Hminus.csv", "PS2017_fig1_k_Hminus", "H- (dissociative attachment)"),
    ("fig1", "OHminus.csv", "PS2017_fig1_k_OHminus", "OH- (dissociative attachment)"),
    ("fig1", "Ominus.csv", "PS2017_fig1_k_Ominus", "O- (dissociative attachment)"),
    (
        "fig2",
        "Petro2017_fig2_line.csv",
        "PS2017_fig2_negion_total",
        "total negative-ion-to-neutral ratio (the single published curve)",
    ),
    ("fig3", "H2Oplus.csv", "PS2017_fig3_frac_H2Oplus", "H2O+ (main axes)"),
    ("fig3", "Hplus.csv", "PS2017_fig3_frac_Hplus", "H+ (inset)"),
    ("fig3", "OHplus.csv", "PS2017_fig3_frac_OHplus", "OH+ (inset)"),
    ("fig3", "Oplus.csv", "PS2017_fig3_frac_Oplus", "O+ (inset)"),
    ("fig5", "ionization.csv", "PS2017_fig5_chi_ionisation", "ionisation"),
    (
        "fig5",
        "electronic_excitation.csv",
        "PS2017_fig5_chi_electronic_excitation",
        "electronic excitation (incl. radiating and nonradiating dissociation)",
    ),
    ("fig5", "vibration.csv", "PS2017_fig5_chi_vibration", "vibration"),
    ("fig5", "rotation.csv", "PS2017_fig5_chi_rotation", "rotation"),
    ("fig5", "elastic.csv", "PS2017_fig5_chi_elastic", "elastic"),
    ("fig5", "attachment.csv", "PS2017_fig5_chi_attachment", "attachment"),
    ("fig7a", "Ar Rm_1.csv", "PS2017_fig7a_etaT_Ar_Rm1", "Ar, R_m = 1"),
    ("fig7a", "Ar Rm_10.csv", "PS2017_fig7a_etaT_Ar_Rm10", "Ar, R_m = 10"),
    ("fig7a", "Ar Rm_100.csv", "PS2017_fig7a_etaT_Ar_Rm100", "Ar, R_m = 100"),
    ("fig7a", "H2O Rm_1.csv", "PS2017_fig7a_etaT_H2O_Rm1", "H2O, R_m = 1"),
    ("fig7a", "H2O Rm_10.csv", "PS2017_fig7a_etaT_H2O_Rm10", "H2O, R_m = 10"),
    ("fig7a", "H2O Rm_100.csv", "PS2017_fig7a_etaT_H2O_Rm100", "H2O, R_m = 100"),
    ("fig7b", "Ar Te_2.csv", "PS2017_fig7b_etaT_Ar_Te2", "Ar, T_e = 2 eV"),
    ("fig7b", "Ar Te_5.csv", "PS2017_fig7b_etaT_Ar_Te5", "Ar, T_e = 5 eV"),
    ("fig7b", "Ar Te_10.csv", "PS2017_fig7b_etaT_Ar_Te10", "Ar, T_e = 10 eV"),
    ("fig7b", "H2O Te_2.csv", "PS2017_fig7b_etaT_H2O_Te2", "H2O, T_e = 2 eV"),
    ("fig7b", "H2O Te_5.csv", "PS2017_fig7b_etaT_H2O_Te5", "H2O, T_e = 5 eV"),
    ("fig7b", "H2O Te_10.csv", "PS2017_fig7b_etaT_H2O_Te10", "H2O, T_e = 10 eV"),
)


# Raw files


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


# Processing


def collapse_repeated_x(xy, ylog, rtol=X_RTOL):
    """Step 4: x values equal to within rtol (relative) become one point.

    Returns (xy_collapsed, list of (x, n, y_lo, y_hi, y_kept)).
    """
    if len(xy) == 0:
        return xy, []
    gid = np.concatenate(([0], np.cumsum(np.diff(xy[:, 0]) > rtol * np.abs(xy[1:, 0]))))
    out, groups = [], []
    for g in range(gid[-1] + 1):
        grp = xy[gid == g]
        if len(grp) == 1:
            out.append(grp[0])
            continue
        y = grp[:, 1]
        if ylog and np.all(y > 0):
            yk = float(np.exp(np.mean(np.log(y))))
        else:
            yk = float(np.mean(y))
        xk = float(np.mean(grp[:, 0]))
        out.append((xk, yk))
        groups.append((xk, len(grp), float(y.min()), float(y.max()), yk))
    return np.array(out, dtype=float), groups


def loo_residuals(x, y, xlog, ylog):
    """Step 7. Leave-one-out residuals about a local quadratic, axis space."""
    u = np.log10(x) if xlog else x.copy()
    v = np.log10(y) if ylog else y.copy()
    n = len(u)
    r = np.full(n, np.nan)
    for i in range(N_SIDE_MIN, n - N_SIDE_MIN):
        nb = np.r_[max(0, i - N_SIDE) : i, i + 1 : min(n, i + 1 + N_SIDE)]
        if np.any(~np.isfinite(v[nb])) or not np.isfinite(v[i]):
            continue
        coef = np.polyfit(u[nb] - u[i], v[nb], 2)
        r[i] = v[i] - coef[-1]
    return r


def process(raw, fig, removals=()):
    """Steps 2-7 (with 4a). Returns (x, y, info dict, log lines)."""
    meta = FIGS[fig]
    xlog, ylog = meta["xscale"] == "log", meta["yscale"] == "log"
    info = {"n_raw": len(raw), "removed": []}
    sorted_raw = bool(np.all(np.diff(raw[:, 0]) >= 0))
    order = np.argsort(raw[:, 0], kind="stable")
    moved = int(np.sum(order != np.arange(len(raw))))
    info["sorted_raw"] = sorted_raw
    log = [
        f"raw rows: {len(raw)}; raw x order monotone: {sorted_raw}"
        + ("" if sorted_raw else f" (sorted: {moved} rows change position)")
    ]
    xy = raw[order]

    for xr, yr, why in removals:  # step 4a
        hit = (xy[:, 0] == xr) & (xy[:, 1] == yr)
        if hit.sum() != 1:
            raise ValueError(
                f"removal ({xr!r}, {yr!r}) matched {int(hit.sum())} rows, expected 1"
            )
        xy = xy[~hit]
        info["removed"].append(f"removed ({xr!r}, {yr!r}): {why}")
        log.append(f"removed ({xr!r}, {yr!r}): {why}")

    keep = np.ones(len(xy), dtype=bool)
    keep[1:] = ~((xy[1:, 0] == xy[:-1, 0]) & (xy[1:, 1] == xy[:-1, 1]))
    info["n_dup"] = int((~keep).sum())
    log.append(f"exact duplicate (x, y) rows dropped: {info['n_dup']}")
    xy = xy[keep]

    xy, groups = collapse_repeated_x(xy, ylog)
    info["n_groups"] = len(groups)
    log.append(
        f"repeated-x groups (x equal to within {X_RTOL:g} relative): "
        f"{len(groups)}; each replaced by its mean y "
        f"({'geometric' if ylog else 'arithmetic'})"
    )
    for xk, n, lo, hi, yk in groups:
        log.append(f"    x = {xk!r}: {n} points, y {lo:.6g} to {hi:.6g}, kept {yk!r}")

    flags = []
    neg = xy[:, 1] < 0 if not ylog else xy[:, 1] <= 0
    for x0, y0 in xy[neg]:
        if ylog:
            flags.append(
                f"non-positive y on a log axis at x = {float(x0)!r} "
                f"(y = {float(y0)!r}): dropped"
            )
        else:
            flags.append(
                f"negative y at x = {float(x0)!r} (y = {float(y0)!r}): clamped to 0"
            )
    if ylog:
        xy = xy[~neg]
    else:
        xy[neg, 1] = 0.0
    x, y = xy[:, 0], xy[:, 1]

    # step 7: digitising error (computed first: it scales the step 6 flags)
    yy = np.where(y > 0, y, np.nan) if ylog else y
    r = loo_residuals(x, yy, xlog, ylog)
    ok = np.isfinite(r)
    rms = float(np.sqrt(np.mean(r[ok] ** 2))) if ok.any() else float("nan")
    info["err"] = rms
    info["n_err"] = int(ok.sum())

    # step 6: significant local extrema other than the global max and min
    v = np.log10(yy) if ylog else yy
    d = np.sign(np.diff(y))
    i_max, i_min = int(np.argmax(y)), int(np.argmin(y))
    n_small = 0
    for i in range(1, len(y) - 1):
        if d[i - 1] != 0 and d[i] != 0 and d[i - 1] != d[i] and i not in (i_max, i_min):
            rev = abs(v[i] - 0.5 * (v[i - 1] + v[i + 1]))
            if not rev > EXTREMUM_K * rms:
                n_small += 1
                continue
            kind = "max" if d[i - 1] > 0 else "min"
            flags.append(
                f"local {kind} at x = {x[i]:.6g} (y = {y[i]:.6g}; "
                f"neighbours {y[i - 1]:.6g}, {y[i + 1]:.6g}; "
                f"reversal {rev / rms:.1f} x RMS)"
            )
    info["n_small_extrema"] = n_small
    for i in np.where(ok & (np.abs(r) > OUTLIER_K * rms))[0]:
        flags.append(
            f"residual outlier at x = {x[i]:.6g} (y = {y[i]:.6g}): "
            f"{r[i] / rms:+.1f} x RMS"
        )
    if ylog:
        rel = 100 * (10**rms - 1)
        info["err_txt"] = f"{rms:.3g} decades (about {rel:.2g} % relative)"
    else:
        info["err_txt"] = f"{rms:.3g} (absolute, {meta['y'][2]})"
    log.append(
        f"digitising error (RMS leave-one-out residual over {info['n_err']} "
        f"interior points): {info['err_txt']}"
    )
    log.append(
        f"local extrema reversing by less than {EXTREMUM_K:g} x RMS "
        f"(within scatter, not flagged): {n_small}"
    )

    info["flags"] = flags
    log.append(f"flags: {len(flags)}")
    log += ["    " + f for f in flags]
    info["n_out"] = len(x)
    info["x_range"] = (float(x[0]), float(x[-1]))
    log.append(f"output rows: {len(x)}; x {float(x[0])!r} to {float(x[-1])!r}")
    return x, y, info, log


# Output


def header(fig, raw_name, curve, info):
    """Citation, units and digitising-error header of an output CSV."""
    meta = FIGS[fig]
    xc, xs, xu = meta["x"]
    yc, ys, yu = meta["y"]
    lines = [
        f"# Petro and Sedwick (2017) {fig.replace('fig', 'Fig. ')}: {curve}.",
        f"# Source: {SOURCE}, digitised",
        f"#   (digitised_raw/{meta['raw_dir']}/{raw_name}).",
        f"# Units: {xc} = {xs} [{xu}]; {yc} = {ys} [{yu}]; "
        f"axes x {meta['xscale']}, y {meta['yscale']}.",
        f"# Digitising error (RMS leave-one-out residual, {info['n_err']} "
        f"interior points): {info['err_txt']}.",
        f"{xc},{yc}",
    ]
    return "\n".join(lines) + "\n"


def write_provenance(done, missing, unmapped):
    """Write PROVENANCE.md: axes, raw-to-output map, removals and flags."""
    L = [
        "# Petro and Sedwick (2017) digitised curves",
        "",
        f"Source: {SOURCE}.",
        "",
        "Digitised with WebPlotDigitizer, calibrated in data units. The raw "
        "exports are in `digitised_raw/`; every other file in this folder is "
        "written by `scripts/convert_ps2017.py`. Used by "
        "`scripts/ps2017_reproduction.py`.",
        "",
        "## Axes",
        "",
        "| figure | raw subfolder | x | y | axis scales |",
        "|---|---|---|---|---|",
    ]
    for fig, meta in FIGS.items():
        L.append(
            f"| {fig.replace('fig', 'Fig. ')} | `{meta['raw_dir']}/` | "
            f"{meta['x'][1]} [{meta['x'][2]}] | {meta['y'][1]} [{meta['y'][2]}] | "
            f"x {meta['xscale']}, y {meta['yscale']} |"
        )
    L += [
        "",
        "Fig. 6 (E_c) is not digitised: Petro and Sedwick give analytic fits, "
        "Eqs. (9) and (10).",
        "",
        "## Raw-to-output map",
        "",
        "| raw file | output file | curve | raw rows | output rows | "
        "digitising error | flags | sha256 (raw) |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for fig, raw_name, canon, curve, sha, info in done:
        L.append(
            f"| `{FIGS[fig]['raw_dir']}/{raw_name}` | `{canon}.csv` | {curve} | "
            f"{info['n_raw']} | {info['n_out']} | "
            f"{info['err_txt']} | {len(info['flags'])} | `{sha}` |"
        )
    L += ["", "## Removed points", ""]
    rem = [
        (canon, r)
        for fig, raw_name, canon, curve, sha, info in done
        for r in info["removed"]
    ]
    L += [f"- `{canon}`: {r}" for canon, r in rem] or ["None."]
    L += ["", "## Flags (kept as digitised unless stated)", ""]
    any_flag = False
    for fig, raw_name, canon, curve, sha, info in done:
        if info["flags"]:
            any_flag = True
            L.append(f"- `{canon}`:")
            L += [f"    - {f}" for f in info["flags"]]
    if not any_flag:
        L.append("None.")
    if missing or unmapped:
        L += [
            "",
            "## Missing and unmapped",
            "",
            "Mapped curves without a raw file: "
            + (", ".join(f"`{c}`" for c in missing) or "none")
            + ".",
            "",
            "Raw CSVs not in the map: "
            + (", ".join(f"`{u}`" for u in unmapped) or "none")
            + ".",
        ]
    L += [
        "",
        "## Fridman detachment curves",
        "",
        "Electron-impact detachment cross sections of H- and O-, Fridman "
        "(2021), Fig. 2.9 (cited by Petro and Sedwick for detachment), "
        "digitised; in "
        "`data/cross_sections/n2o/digitised_raw/Fridman2021_fig2_9/`, two "
        "columns without a header, E [eV] and sigma [cm^2]. Used by "
        "`scripts/ps2017_reproduction.py`.",
        "",
    ]
    (PS_DIR / "PROVENANCE.md").write_bytes("\n".join(L).encode("ascii"))


def convert():
    PS_DIR.mkdir(parents=True, exist_ok=True)
    done, missing, curves = [], [], {}
    mapped = set()
    for fig, raw_name, canon, curve in CURVES:
        path = RAW_DIR / FIGS[fig]["raw_dir"] / raw_name
        mapped.add(path.resolve())
        if not path.exists():
            missing.append(canon)
            print(f"== MISSING {FIGS[fig]['raw_dir']}/{raw_name} -> {canon}: skipped")
            continue
        sha = hashlib.sha256(path.read_bytes()).hexdigest()
        x, y, info, log = process(read_raw(path), fig, REMOVED.get((fig, raw_name), ()))
        body = "".join(f"{float(a)!r},{float(b)!r}\n" for a, b in zip(x, y))
        text = header(fig, raw_name, curve, info) + body
        (PS_DIR / f"{canon}.csv").write_bytes(text.encode("ascii"))
        print(
            f"== {FIGS[fig]['raw_dir']}/{raw_name} -> {canon}.csv (sha256 {sha[:12]}...)"
        )
        for line in log:
            print("  " + line)
        done.append((fig, raw_name, canon, curve, sha, info))
        curves[canon] = (x, y, info)
    unmapped = sorted(
        p.relative_to(RAW_DIR).as_posix()
        for p in RAW_DIR.rglob("*.csv")
        if p.resolve() not in mapped
    )
    for u in unmapped:
        print(f"== UNMAPPED raw file {u}: not in CURVES, skipped")
    write_provenance(done, missing, unmapped)
    print(
        f"\n{len(done)} curves converted, {len(missing)} missing, "
        f"{len(unmapped)} unmapped; PROVENANCE.md written."
    )
    return curves


if __name__ == "__main__":
    convert()
