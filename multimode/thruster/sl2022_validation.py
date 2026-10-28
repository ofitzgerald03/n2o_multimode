"""
Exit-velocity check against Sheppard and Little (2022): per-point factor
k = u_e,meas / u_e,model(T_e) for Ar, H2O and Xe, and the offsets between
propellants, from the digitised data in data/sl2022/.
"""

# For each measured u_e point, T_e is the measured T_e curve of the same
# mass-flow group, interpolated linearly in mdot/P; points outside that
# curve's range are dropped. offset = mean(k_num) / mean(k_den) - 1, with a
# first-order standard error, a per-point bootstrap and a bootstrap over
# mass-flow groups.

from pathlib import Path

import numpy as np
import pandas as pd

from ..config import PARAMS, PATHS
from ..constants import AMU
from .exhaust import exit_velocity_PS
from .pipeline import run_pipeline

SEED = 20261007
X = "mdot_over_P_ug_per_J"

PROPELLANTS = ("Ar", "H2O", "Xe")
PAIRS = (("Xe", "Ar"), ("H2O", "Xe"), ("H2O", "Ar"))  # (numerator, denominator)


# Digitised data


def _load_csv(path):
    return pd.read_csv(path, comment="#", skipinitialspace=True)


def load_digitised(input_dir, pair=("Ar", "H2O")):
    """T_e and u_e tables of pair = (den, num) from SL2022_{prop}_{Te,ue}.csv;
    columns mdot_mg_s, mdot_over_P_ug_per_J and the quantity."""
    d = Path(input_dir)
    den, num = pair
    return {
        "den_Te": _load_csv(d / f"SL2022_{den}_Te.csv"),
        "den_ue": _load_csv(d / f"SL2022_{den}_ue.csv"),
        "num_Te": _load_csv(d / f"SL2022_{num}_Te.csv"),
        "num_ue": _load_csv(d / f"SL2022_{num}_ue.csv"),
    }


def _pair_within_group(Te_g, ue_g):
    """Pair each T_e point with the nearest unused u_e point in mdot/P."""
    Te = Te_g.reset_index(drop=True)
    ue = ue_g.reset_index(drop=True)
    used = np.zeros(len(ue), dtype=bool)
    cand_x = ue["mdot_over_P_ug_per_J"].to_numpy()
    rows = []
    for _, te in Te.iterrows():
        dist = np.abs(cand_x - te["mdot_over_P_ug_per_J"])
        dist[used] = np.inf
        if not np.isfinite(dist).any():
            continue
        j = int(np.argmin(dist))
        used[j] = True
        rows.append(
            {
                "Te_eV": float(te["Te_eV"]),
                "ue_m_s": float(ue.iloc[j]["ue_m_s"]),
                "mdot_over_P_diff": float(dist[j]),
            }
        )
    return pd.DataFrame(rows)


def pair_within_propellant(Te_df, ue_df):
    """One-to-one nearest (T_e, u_e) pairs over all mass-flow groups."""
    chunks = []
    groups = sorted(set(Te_df["mdot_mg_s"]) | set(ue_df["mdot_mg_s"]))
    for mdot in groups:
        Te_g = Te_df[Te_df["mdot_mg_s"] == mdot]
        ue_g = ue_df[ue_df["mdot_mg_s"] == mdot]
        if len(Te_g) == 0 or len(ue_g) == 0:
            continue
        paired = _pair_within_group(Te_g, ue_g)
        paired.insert(0, "mdot_mg_s", mdot)
        chunks.append(paired)
    return (
        pd.concat(chunks, ignore_index=True)
        if chunks
        else pd.DataFrame(columns=["mdot_mg_s", "Te_eV", "ue_m_s"])
    )


def overlap_range(den_df, num_df):
    """T_e range covered by both propellants."""
    lo = max(den_df["Te_eV"].min(), num_df["Te_eV"].min())
    hi = min(den_df["Te_eV"].max(), num_df["Te_eV"].max())
    return float(lo), float(hi)


# Model exit velocity


def model_ue():
    """{name: callable(Te) -> u_e [m/s]}: Ar and Xe single species, H2O with
    the common Bohm velocity, 'H2O_ss' pure H2O+ (Eq. (7) of Petro and
    Sedwick 2017); u_e / sqrt(Te) is interpolated linearly in Te."""
    r = run_pipeline(PARAMS, propellants=("Xe", "Ar", "H2O"), anchor_propellant="Xe")
    Te = r["Te"]

    def on_grid(ue):
        c = ue / np.sqrt(Te)
        return lambda t: np.interp(t, Te, c, left=np.nan, right=np.nan) * np.sqrt(t)

    m_ar = float(r["Ar"]["q_A"]["m_i_eff"][0])
    m_xe = float(r["Xe"]["q_A"]["m_i_eff"][0])
    return {
        "Ar": on_grid(r["Ar"]["q_A"]["u_e"]),
        "Xe": on_grid(r["Xe"]["q_A"]["u_e"]),
        "H2O": on_grid(r["H2O"]["q_B"]["u_e"]),
        "H2O_ss": lambda t: exit_velocity_PS(np.asarray(t, float), 18.0 * AMU),
        "_check": dict(
            Ar_qA_equals_PS=float(
                np.max(np.abs(r["Ar"]["q_A"]["u_e"] / exit_velocity_PS(Te, m_ar) - 1))
            ),
            Xe_qA_equals_PS=float(
                np.max(np.abs(r["Xe"]["q_A"]["u_e"] / exit_velocity_PS(Te, m_xe) - 1))
            ),
            ratio_B_anchor=float(
                np.interp(
                    PARAMS["Te_ref"], Te, r["H2O"]["q_B"]["u_e"] / r["Ar"]["q_A"]["u_e"]
                )
            ),
            ratio_B_range=[
                float(v)
                for v in (lambda q: (q.min(), q.max()))(
                    r["H2O"]["q_B"]["u_e"] / r["Ar"]["q_A"]["u_e"]
                )
            ],
            ratio_ss=float(
                exit_velocity_PS(10.0, 18.0 * AMU) / exit_velocity_PS(10.0, m_ar)
            ),
        ),
    }


# Per-point estimator


def te_at_ue_points(Te_df, ue_df, mode="interp", logx=False):
    """T_e at each u_e point from the T_e curve of its mass-flow group.

    mode 'interp' drops points outside the curve's range, 'clamp' takes the
    nearest end point; logx interpolates in log(mdot/P).
    """
    rows = []
    for g, ue_g in ue_df.groupby("mdot_mg_s", sort=True):
        te_g = Te_df[Te_df["mdot_mg_s"] == g].sort_values(X)
        xs, ts = te_g[X].to_numpy(), te_g["Te_eV"].to_numpy()
        fx = np.log if logx else (lambda v: v)
        for _, u in ue_g.iterrows():
            x = float(u[X])
            inside = xs[0] <= x <= xs[-1]
            if inside:
                te, rule = float(np.interp(fx(x), fx(xs), ts)), "interp"
            elif mode == "clamp":
                te, rule = float(ts[0] if x < xs[0] else ts[-1]), "clamped"
            else:
                te, rule = np.nan, "dropped (outside T_e x-range)"
            rows.append(
                dict(
                    mdot_mg_s=g,
                    x=x,
                    ue_m_s=float(u["ue_m_s"]),
                    Te_eV=te,
                    rule=rule,
                    x_range=(float(xs[0]), float(xs[-1])),
                )
            )
    return pd.DataFrame(rows)


def nearest_pairing(Te_df, ue_df):
    """T_e at each u_e point by one-to-one nearest pairing."""
    p = pair_within_propellant(Te_df, ue_df)
    return pd.DataFrame(
        dict(
            mdot_mg_s=p["mdot_mg_s"],
            x=np.nan,
            ue_m_s=p["ue_m_s"],
            Te_eV=p["Te_eV"],
            rule="nearest pairing",
        )
    )


def stats(k):
    k = np.asarray(k, float)
    n = len(k)
    return dict(
        n=n,
        mean=float(k.mean()),
        sd=float(k.std(ddof=1)),
        se=float(k.std(ddof=1) / np.sqrt(n)),
    )


def estimate(kA, kW, gA=None, gW=None, B=20000, seed=SEED):
    """Offset mean(kW) / mean(kA) - 1, its first-order SE and 95 % bootstrap
    intervals (per point, and by mass-flow group when gA, gW are given)."""
    sA, sW = stats(kA), stats(kW)
    R = sW["mean"] / sA["mean"]
    se = R * np.sqrt((sW["se"] / sW["mean"]) ** 2 + (sA["se"] / sA["mean"]) ** 2)
    rng = np.random.default_rng(seed)
    kA, kW = np.asarray(kA), np.asarray(kW)
    bs = rng.choice(kW, (B, len(kW))).mean(1) / rng.choice(kA, (B, len(kA))).mean(1)
    out = dict(
        Ar=sA,
        H2O=sW,
        ratio_of_means=float(R),
        offset=float(R - 1),
        se=float(se),
        se_interval95=[float(R - 1 - 1.96 * se), float(R - 1 + 1.96 * se)],
        boot_B=B,
        boot_seed=seed,
        boot_interval95=[
            float(np.percentile(bs, 2.5) - 1),
            float(np.percentile(bs, 97.5) - 1),
        ],
    )
    if gA is not None:
        # resample whole mass-flow groups within each propellant
        def clus(k, g):
            labs = np.unique(g)
            idx = [np.flatnonzero(g == l) for l in labs]
            sums = np.array([k[i].sum() for i in idx])
            cnt = np.array([len(i) for i in idx])
            pick = rng.integers(0, len(labs), (B, len(labs)))
            return sums[pick].sum(1) / cnt[pick].sum(1)

        cb = clus(kW, np.asarray(gW)) / clus(kA, np.asarray(gA))
        out["cluster_boot_interval95"] = [
            float(np.percentile(cb, 2.5) - 1),
            float(np.percentile(cb, 97.5) - 1),
        ]
        out["n_groups"] = dict(Ar=int(len(np.unique(gA))), H2O=int(len(np.unique(gW))))
    return out


def sensitivity_variants(B=20000):
    """H2O / Ar offset under each variant of the method.

    Returns ({variant: estimate}, {variant: tables}, model checks).
    """
    raw = load_digitised(PATHS["sl2022"])
    M = model_ue()
    res, tables = {}, {}

    def k_of(df, model):
        d = df.dropna(subset=["Te_eV"]).copy()
        d["ue_model"] = model(d["Te_eV"].to_numpy())
        d["k"] = d["ue_m_s"] / d["ue_model"]
        return d

    variants = {
        "baseline (interp in x, drop outside, multi-species H2O)": (
            "interp",
            False,
            "H2O",
            None,
        ),
        "single-species H2O model": ("interp", False, "H2O_ss", None),
        "end-point clamp instead of drop": ("clamp", False, "H2O", None),
        "interpolation in log x": ("interp", True, "H2O", None),
        "one-to-one nearest pairing": ("nearest", False, "H2O", None),
        "only points inside the T_e overlap": ("interp", False, "H2O", "overlap"),
    }
    for name, (mode, logx, wmodel, restrict) in variants.items():
        if mode == "nearest":
            A0 = nearest_pairing(raw["den_Te"], raw["den_ue"])
            W0 = nearest_pairing(raw["num_Te"], raw["num_ue"])
        else:
            A0 = te_at_ue_points(raw["den_Te"], raw["den_ue"], mode, logx)
            W0 = te_at_ue_points(raw["num_Te"], raw["num_ue"], mode, logx)
        A, W = k_of(A0, M["Ar"]), k_of(W0, M[wmodel])
        if restrict == "overlap":
            lo, hi = overlap_range(A, W)
            A = A[(A.Te_eV >= lo) & (A.Te_eV <= hi)]
            W = W[(W.Te_eV >= lo) & (W.Te_eV <= hi)]
        e = estimate(
            A["k"].to_numpy(),
            W["k"].to_numpy(),
            A["mdot_mg_s"].to_numpy(),
            W["mdot_mg_s"].to_numpy(),
            B=B,
        )
        e["Te_range"] = dict(
            Ar=[float(A.Te_eV.min()), float(A.Te_eV.max())],
            H2O=[float(W.Te_eV.min()), float(W.Te_eV.max())],
        )
        e["dropped"] = dict(
            Ar=int(A0["Te_eV"].isna().sum()), H2O=int(W0["Te_eV"].isna().sum())
        )
        # least-squares slope of k against T_e, relative to mean k [1/eV]
        e["k_trend_rel_per_eV"] = {
            lab: float(np.polyfit(d.Te_eV, d.k, 1)[0] / d.k.mean())
            for lab, d in (("Ar", A), ("H2O", W))
        }
        res[name] = e
        tables[name] = (A0, W0, A, W)
    return res, tables, M["_check"]


def per_point_k(B=20000):
    """Per-point k of each propellant and the offset of each pair in PAIRS.

    Returns {'points', 'per_propellant', 'pairs', 'model_checks'}.
    """
    M = model_ue()
    pts, dropped = {}, {}
    for p in PROPELLANTS:
        raw = load_digitised(PATHS["sl2022"], pair=(p, p))
        d0 = te_at_ue_points(raw["den_Te"], raw["den_ue"], "interp", False)
        dropped[p] = int(d0["Te_eV"].isna().sum())
        d = d0.dropna(subset=["Te_eV"]).copy()
        d["ue_model"] = M[p](d["Te_eV"].to_numpy())
        d["k"] = d["ue_m_s"] / d["ue_model"]
        pts[p] = d
    per = {
        p: dict(
            stats(d["k"]),
            Te_range=[float(d.Te_eV.min()), float(d.Te_eV.max())],
            dropped=dropped[p],
            groups=sorted(float(g) for g in set(d["mdot_mg_s"])),
        )
        for p, d in pts.items()
    }
    pairs = {}
    for num, den in PAIRS:
        e = estimate(
            pts[den]["k"].to_numpy(),
            pts[num]["k"].to_numpy(),
            pts[den]["mdot_mg_s"].to_numpy(),
            pts[num]["mdot_mg_s"].to_numpy(),
            B=B,
        )
        e = {k: v for k, v in e.items() if k not in ("Ar", "H2O")}
        # estimate() labels the denominator 'Ar' and the numerator 'H2O'
        e["n_groups"] = {den: e["n_groups"]["Ar"], num: e["n_groups"]["H2O"]}
        e["n"] = {den: len(pts[den]), num: len(pts[num])}
        pairs[f"{num}/{den}"] = e
    return dict(points=pts, per_propellant=per, pairs=pairs, model_checks=M["_check"])
