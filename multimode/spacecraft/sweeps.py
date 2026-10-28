"""
Sweeps of the spacecraft-level model over f_EP, eta_u, warning time and
transfer window, the one-at-a-time sensitivity and the corner grid.
"""

import math

T_WARN_SET_S = (7200.0, 21600.0, 86400.0, 259200.0)  # 2 h, 6 h, 1 d, 3 d
# Six-hourly points between 24 h and 72 h
T_WARN_TAIL_S = tuple(3600.0 * h for h in (30, 36, 42, 48, 54, 60, 66))
# GEO transfer windows [yr]: 1 and 2 weeks, 1, 3, 6, 12 and 24 months
T_TRANSFER_MAX_SET_YR = (1.0 / 52.0, 2.0 / 52.0, 1.0 / 12.0, 0.25, 0.5, 1.0, 2.0)
ETA_U_RANGE = (0.05, 0.61)
ETA_U_OPTIMISTIC = (0.61, 0.85)


def sweep_f_ep(arch_keys, mission, n_points=21):
    """Wet mass and verdict against f_EP for each architecture; single-mode
    architectures give one point. Diverged closures carry wet = nan.

    Returns {arch_key: {'records': [{f_ep, wet, converged, diverged,
    excluded, reasons}, ...], 'f_ep_derived': f_EP of the nominal
    allocation}}.
    """
    from math import nan
    from .allocation import nominal_allocation, sweep_allocation, derived_f_ep
    from .architectures import get_architecture
    from .closure import close
    from .feasibility import assess

    out = {}
    for key in arch_keys:
        arch = get_architecture(key)
        if arch.has_ep and arch.has_chem:
            fs = [i / (n_points - 1) for i in range(n_points)]
        elif arch.has_ep:
            fs = [1.0]  # single-electric: a point
        else:
            fs = [0.0]  # single-chemical
        records = []
        for f in fs:
            alloc = sweep_allocation(key, mission, f)
            r = close(key, mission, alloc)
            if r.diverged or not r.converged:
                records.append(
                    dict(
                        f_ep=f,
                        wet=nan,
                        converged=r.converged,
                        diverged=r.diverged,
                        excluded=True,
                        reasons=[
                            "closure diverged (no fixed point)"
                            if r.diverged
                            else "closure not converged"
                        ],
                    )
                )
                continue
            rep = assess(key, mission, alloc, r)
            records.append(
                dict(
                    f_ep=f,
                    wet=r.wet,
                    converged=True,
                    diverged=False,
                    excluded=rep.excluded,
                    reasons=rep.reasons,
                )
            )
        try:
            f_der = derived_f_ep(mission, nominal_allocation(key, mission))
        except KeyError:
            f_der = None
        out[key] = {"records": records, "f_ep_derived": f_der}
    return out


def _eta_u_default_grid():
    """eta_u from 0.05 to 0.85 in steps of 0.02."""
    return tuple(i / 100.0 for i in range(5, 86, 2))


def _close_arch_at_eta_u(mission, eta_u, allocation, arch_key="multimode"):
    """Closure and assessment of an electrodeless architecture with its
    electric unit at eta_u. Returns (arch, result, report or None, wet or
    nan).
    """
    from math import nan
    from .architectures import with_ecr_eta_u
    from .closure import close
    from .feasibility import assess

    arch = with_ecr_eta_u(arch_key, eta_u)
    r = close(arch_key, mission, allocation, arch=arch)
    if r.diverged or not r.converged:
        return arch, r, None, nan
    rep = assess(arch_key, mission, allocation, r, arch=arch)
    return arch, r, rep, r.wet


def min_viable_eta_u(
    mission, criterion="own_chem", lo=0.02, hi=None, tol=1e-5, arch_key="multimode"
):
    """Smallest eta_u at which the architecture at its nominal allocation is
    no heavier than its chemical alternative, by bisection.

    criterion 'own_chem': the same architecture at f_EP = 0;
    'single_chem': the single-chemical architecture.
    Returns the threshold, lo if viable everywhere, or None if not viable up
    to hi (default: eta_u parity, 0.61) or if the reference does not
    converge.
    """
    from math import inf, isnan
    from .allocation import nominal_allocation, sweep_allocation
    from .closure import close

    if hi is None:
        from ..thruster.ecr_unit import eta_u_parity

        hi = eta_u_parity()

    alloc_ep = nominal_allocation(arch_key, mission)
    if criterion == "own_chem":
        r_ref = close(arch_key, mission, sweep_allocation(arch_key, mission, 0.0))
        ref = r_ref.wet if r_ref.converged else float("nan")
    elif criterion == "single_chem":
        r_ref = close(
            "single_chem", mission, nominal_allocation("single_chem", mission)
        )
        ref = r_ref.wet if r_ref.converged else float("nan")
    else:
        raise ValueError(f"unknown criterion {criterion!r}")
    if isnan(ref):
        return None

    def wet(eta_u):
        _, _, _, w = _close_arch_at_eta_u(mission, eta_u, alloc_ep, arch_key)
        return inf if isnan(w) else w

    if wet(hi) > ref:
        return None  # never viable in the range
    if wet(lo) <= ref:
        return lo  # viable everywhere probed
    a, b = lo, hi  # wet(a) > ref >= wet(b)
    while b - a > tol:
        m = 0.5 * (a + b)
        if wet(m) <= ref:
            b = m
        else:
            a = m
    return 0.5 * (a + b)


def sweep_eta_u(mission, eta_u_grid=None, arch_key="multimode"):
    """Wet mass of an electrodeless architecture at its nominal allocation
    against eta_u of N2O, with the eta_u-independent references.

    Returns {'grid', 'records': [{eta_u, in_range, wet, converged,
    diverged, excluded, reasons, isp, n_ep}, ...], 'references':
    {'own_f0', 'single_chem', 'single_elec', 'hybrid'} (and 'multimode_f0'
    for the multimode), 'min_viable': {'own_chem', 'single_chem'}}.
    """
    from math import nan
    from .allocation import nominal_allocation, sweep_allocation
    from .architectures import get_architecture
    from .closure import close

    grid = _eta_u_default_grid() if eta_u_grid is None else tuple(eta_u_grid)
    alloc_ep = nominal_allocation(arch_key, mission)

    records = []
    for eu in grid:
        arch, r, rep, wet = _close_arch_at_eta_u(mission, eu, alloc_ep, arch_key)
        if rep is None:
            reasons = (
                ["closure diverged (no fixed point)"]
                if r.diverged
                else [
                    "closure not converged within max_iter "
                    "(contraction ~ 1 at this eta_u; stack "
                    "unbounded-for-purpose)"
                ]
            )
        else:
            reasons = rep.reasons
        records.append(
            dict(
                eta_u=eu,
                in_range=(ETA_U_RANGE[0] <= eu <= ETA_U_RANGE[1]),
                wet=wet,
                converged=r.converged,
                diverged=r.diverged,
                excluded=(True if rep is None else rep.excluded),
                reasons=reasons,
                isp=arch.ep_unit.isp,
                n_ep=(r.n_units["ep"] if r.converged else None),
            )
        )

    # References independent of eta_u
    r_f0 = close(arch_key, mission, sweep_allocation(arch_key, mission, 0.0))
    r_sc = close("single_chem", mission, nominal_allocation("single_chem", mission))
    refs = {
        "own_f0": r_f0.wet if r_f0.converged else nan,
        "single_chem": r_sc.wet if r_sc.converged else nan,
    }
    if arch_key == "multimode":
        refs["multimode_f0"] = refs["own_f0"]
    for key in ("single_elec", "hybrid"):
        arch = get_architecture(key)
        assert arch.ep_propellant != "n2o", (
            f"{key} draws N2O for EP — the eta_u-invariance assumption "
            f"of this sweep no longer holds; recompute per point."
        )
        r = close(key, mission, nominal_allocation(key, mission))
        refs[key] = r.wet if r.converged else nan

    return {
        "grid": grid,
        "records": records,
        "references": refs,
        "min_viable": {
            "own_chem": min_viable_eta_u(mission, "own_chem", arch_key=arch_key),
            "single_chem": min_viable_eta_u(mission, "single_chem", arch_key=arch_key),
        },
    }


def ep_max_units(unit, mission, arr_degr=0.0):
    """Largest electric unit count whose added power-system mass stays within
    the power-system limit of the bus, at array power loss arr_degr.
    """
    from .mass_models import eps_increment, eps_mass_cap

    cap = eps_mass_cap(mission.bus_mass)
    n = 0
    while (
        n < 100000
        and eps_increment((n + 1) * unit.p_dc, mission.bus_class, arr_degr) <= cap
    ):
        n += 1
    return n


def ep_transfer_allocation(arch_key, mission, s_tol=1e-5, arch=None):
    """Nominal allocation with the largest electric share s of the transfer
    that passes the power-system limit (maximum-electric transfer); a chemical
    mode flies the rest. s is found by bisection on the full closure; where
    s = 0 fails (no fixed point), the largest passing share on a 0.05 grid
    is refined. The single-electric architecture flies the whole transfer
    electrically.

    Returns (allocation, {'s_ep', 'forced_all_ep', 'n_max',
    'no_feasible_s'}).
    """
    from .allocation import nominal_allocation
    from .architectures import get_architecture
    from .closure import close
    from .feasibility import assess

    arch = get_architecture(arch_key) if arch is None else arch
    if not arch.has_ep:
        raise ValueError(
            f"{arch_key} has no EP mode; the EP-transfer sweep mode does not apply"
        )
    tr = mission.get_class("transfer")
    base = nominal_allocation(arch_key, mission)
    # unit ceiling of a fully electric transfer
    from .mass_models import transfer_array_derate

    n_max = ep_max_units(arch.ep_unit, mission, transfer_array_derate(mission, 1.0))

    def alloc_at(s):
        a = dict(base)
        a[tr.key] = {
            "chem": (1.0 - s) * tr.dv_impulsive,
            "ep": s * tr.dv_impulsive * tr.ep_loss_ratio,
        }
        return a

    if not arch.has_chem:  # single_elec: no top-up
        return alloc_at(1.0), dict(
            s_ep=1.0, forced_all_ep=True, n_max=n_max, no_feasible_s=False
        )

    def power_ok(s):
        a = alloc_at(s)
        r = close(arch_key, mission, a, arch=arch)
        if r.diverged or not r.converged:
            return False
        return assess(arch_key, mission, a, r, arch=arch).power_ok

    def _bisect_up(lo, hi):
        """Largest passing s in [lo, hi] given power_ok(lo) True,
        power_ok(hi) False."""
        while hi - lo > s_tol:
            mid = 0.5 * (lo + hi)
            if power_ok(mid):
                lo = mid
            else:
                hi = mid
        return lo  # the feasible side of the boundary

    no_feasible_s = False
    if power_ok(1.0):
        s_star = 1.0
    elif power_ok(0.0):
        s_star = _bisect_up(0.0, 1.0)
    else:
        # s = 0 fails: scan for the largest passing share, then bisect
        n_scan = 20
        grid = [i / n_scan for i in range(n_scan + 1)]
        passing = [s for s in grid if power_ok(s)]
        if not passing:
            s_star = 0.0
            no_feasible_s = True
        else:
            s_lo = max(passing)
            hi = min(s_lo + 1.0 / n_scan, 1.0)
            s_star = s_lo if hi <= s_lo or power_ok(hi) else _bisect_up(s_lo, hi)
    return alloc_at(s_star), dict(
        s_ep=s_star, forced_all_ep=False, n_max=n_max, no_feasible_s=no_feasible_s
    )


def mono_corner_grid():
    """Joint 3 x 3 grid over the monopropellant head mass and Isp (bracket
    low, nominal, bracket high), LEO mission, all five architectures.

    Returns {'levels', 'records': [{head_level, isp_level, m_head_mono,
    isp_mono, corner, wets, converged, d_mmm_mm, lightest}, ...], 'corner':
    {'computed', 'superposed', 'interaction'} at the unfavourable corner
    (head high, Isp low), 'lightest_everywhere'}.
    """
    from math import nan, isnan
    from .allocation import nominal_allocation
    from .architectures import ARCHITECTURE_ORDER, invalidate_cache
    from .closure import close
    from .mission import leo_mission
    from .parameters import REGISTER, override

    levels = {}
    for key in ("m_head_mono", "isp_mono"):
        p = REGISTER[key]
        levels[key] = {"lo": p.range[0], "nom": p.value, "hi": p.range[1]}

    corner_names = {
        ("hi", "lo"): "unfavourable (head hi, isp lo)",
        ("lo", "hi"): "favourable (head lo, isp hi)",
    }
    records = []
    for hl in ("lo", "nom", "hi"):
        for il in ("lo", "nom", "hi"):
            ov = {
                "m_head_mono": levels["m_head_mono"][hl],
                "isp_mono": levels["isp_mono"][il],
            }
            try:
                with override(**ov):
                    invalidate_cache()
                    m = leo_mission()
                    wets, conv = {}, {}
                    for k in ARCHITECTURE_ORDER:
                        r = close(k, m, nominal_allocation(k, m))
                        conv[k] = r.converged
                        wets[k] = r.wet if r.converged else nan
            finally:
                invalidate_cache()
            finite = {k: w for k, w in wets.items() if not isnan(w)}
            records.append(
                dict(
                    head_level=hl,
                    isp_level=il,
                    m_head_mono=ov["m_head_mono"],
                    isp_mono=ov["isp_mono"],
                    corner=corner_names.get((hl, il), ""),
                    wets=wets,
                    converged=conv,
                    d_mmm_mm=wets["mono_multimode"] - wets["multimode"],
                    lightest=(min(finite, key=finite.get) if finite else None),
                )
            )

    d = {(r["head_level"], r["isp_level"]): r["d_mmm_mm"] for r in records}
    d_nom = d[("nom", "nom")]
    superposed = d_nom + (d[("hi", "nom")] - d_nom) + (d[("nom", "lo")] - d_nom)
    corner = {
        "computed": d[("hi", "lo")],
        "superposed": superposed,
        "interaction": d[("hi", "lo")] - superposed,
    }
    return {
        "levels": levels,
        "records": records,
        "corner": corner,
        "lightest_everywhere": all(r["lightest"] == "mono_multimode" for r in records),
    }


def _t_warn_default_grid(t_warn_set=T_WARN_SET_S, step_s=1800.0):
    """Warning-time grid [s]: the set t_warn_set, half-hourly from 2 h to
    24 h, and six-hourly from 30 h to 66 h.
    """
    lo, hi = 7200.0, 86400.0
    n = int(round((hi - lo) / step_s))
    fine = {lo + i * step_s for i in range(n + 1)}
    return tuple(sorted(fine | set(t_warn_set) | set(T_WARN_TAIL_S)))


def sweep_t_warn(
    arch_keys,
    mission=None,
    t_warn_set=T_WARN_SET_S,
    grid_s=None,
    locate_flip=True,
    flip_tol_s=60.0,
):
    """Wet mass, verdict and collision-avoidance quantities against warning
    time, LEO mission; the mission and nominal allocation are rebuilt at each
    point. With locate_flip the warning time at which single-electric
    becomes feasible is refined by bisection.

    Returns {'grid', 'records': {arch: [{t_warn_s, wet, converged, diverged,
    excluded, feasible, reasons, ca_factor, drift_ok, ca_prop, n_ep}, ...]},
    'optimal': [{t_warn_s, arch, wet}] (lightest feasible design), 'flip_s'}.
    """
    from math import nan
    from .allocation import nominal_allocation
    from .architectures import CONTENDERS
    from .closure import close
    from .feasibility import assess
    from .mission import leo_mission

    if mission is not None and mission.key != "leo":
        raise ValueError(f"sweep_t_warn needs the LEO mission; got {mission.key!r}")
    grid = _t_warn_default_grid(t_warn_set) if grid_s is None else tuple(sorted(grid_s))

    def _point(arch_key, tw):
        """One (architecture, t_warn) evaluation — mission and
        allocation rebuilt per point."""
        m_tw = leo_mission(t_warn_s=tw)
        alloc = nominal_allocation(arch_key, m_tw)
        r = close(arch_key, m_tw, alloc)
        if r.diverged or not r.converged:
            return dict(
                t_warn_s=tw,
                wet=nan,
                converged=r.converged,
                diverged=r.diverged,
                excluded=True,
                feasible=False,
                reasons=[
                    "closure diverged (no fixed point)"
                    if r.diverged
                    else "closure not converged"
                ],
                ca_factor=None,
                drift_ok=None,
                ca_prop=nan,
                n_ep=r.n_units["ep"],
            )
        rep = assess(arch_key, m_tw, alloc, r)
        ca = r.ca_charge
        return dict(
            t_warn_s=tw,
            wet=r.wet,
            converged=True,
            diverged=False,
            excluded=rep.excluded,
            feasible=not rep.excluded,
            reasons=rep.reasons,
            ca_factor=(None if ca is None else ca["factor"]),
            drift_ok=(None if ca is None else ca["drift_ok"]),
            ca_prop=sum(r.per_class_prop["ca"].values()),
            n_ep=r.n_units["ep"],
        )

    records = {k: [_point(k, tw) for tw in grid] for k in arch_keys}

    # Lightest feasible design per grid point
    optimal = []
    for i, tw in enumerate(grid):
        best = None
        for k in arch_keys:
            if k not in CONTENDERS:
                continue
            rec = records[k][i]
            if rec["feasible"] and (best is None or rec["wet"] < best[1]):
                best = (k, rec["wet"])
        optimal.append(
            {
                "t_warn_s": tw,
                "arch": None if best is None else best[0],
                "wet": nan if best is None else best[1],
            }
        )

    # Single-electric limit, by bisection
    flip = None
    if locate_flip and "single_elec" in arch_keys:
        feas = [rec["feasible"] for rec in records["single_elec"]]
        for i in range(1, len(grid)):
            if not feas[i - 1] and feas[i]:
                lo, hi = grid[i - 1], grid[i]
                while hi - lo > flip_tol_s:
                    mid = 0.5 * (lo + hi)
                    if _point("single_elec", mid)["feasible"]:
                        hi = mid
                    else:
                        lo = mid
                flip = 0.5 * (lo + hi)
                break

    return {"grid": grid, "records": records, "optimal": optimal, "flip_s": flip}


def ca_flip_band(flip_tol_s=60.0, grid_s=None):
    """Single-electric warning-time limit at the nominal inputs and at the
    bracket ends of dv_ca_budget_1d (the band) and of w_duty_leo, with the
    square-root scaling flip ~ sqrt(budget / w) for comparison.

    Returns {'nominal', 'band', 'w_sensitivity', 'band_lo_s', 'band_hi_s'}.
    """
    from .architectures import invalidate_cache
    from .parameters import REGISTER, override, value

    if grid_s is None:
        # hourly, 2 h to 24 h
        grid_s = tuple(3600.0 * h for h in range(2, 25))

    def _flip():
        sw = sweep_t_warn(
            ("single_elec",), grid_s=grid_s, locate_flip=True, flip_tol_s=flip_tol_s
        )
        return sw["flip_s"]

    flip_nom = _flip()
    b_nom = value("dv_ca_budget_1d")
    w_nom = value("w_duty_leo")

    def _edges(key):
        rng = REGISTER[key].range
        out = []
        for v in rng:
            with override(**{key: v}):
                invalidate_cache()
                f = _flip()
            invalidate_cache()
            if key == "dv_ca_budget_1d":
                sqrt_s = flip_nom * math.sqrt(v / b_nom)
            else:  # w_duty_leo
                sqrt_s = flip_nom * math.sqrt(w_nom / v)
            out.append(
                {"param": key, "value": v, "flip_s": f, "sqrt_scaling_s": sqrt_s}
            )
        return out

    band = _edges("dv_ca_budget_1d")
    w_sens = _edges("w_duty_leo")
    flips = [e["flip_s"] for e in band if e["flip_s"] is not None]
    return {
        "nominal": flip_nom,
        "band": band,
        "w_sensitivity": w_sens,
        "band_lo_s": min(flips) if flips else None,
        "band_hi_s": max(flips) if flips else None,
    }


def _transfer_commission_s(arch, mission, result):
    """Wall-clock time of the transfer [s], chemical then electric."""
    cls = mission.get_class("transfer")
    t = 0.0
    for mode in ("chem", "ep"):
        mp = result.per_class_prop[cls.key][mode]
        if mp <= 0.0:
            continue
        unit = arch.chem_unit if mode == "chem" else arch.ep_unit
        n = result.n_units[mode]
        t += (mp / (max(n, 1) * unit.mdot)) / cls.w_duty
    return t


def sweep_t_transfer_max(
    mission=None, t_set=T_TRANSFER_MAX_SET_YR, arch_keys=None, ep_transfer=True
):
    """Designs of the GEO mission against the transfer window; the mission is
    rebuilt at each window.

    Returns {'t_set_yr', 'records' (time-optimal transfer: wet, n_ep, peak
    DC power, added power-system mass, electronics, transfer time, verdict),
    'ep_transfer' (maximum-electric transfer, electric architectures; empty with
    ep_transfer False), 'asymmetry' (closed-form check per window)}.
    """
    from math import nan
    from .allocation import nominal_allocation
    from .architectures import ARCHITECTURE_ORDER, get_architecture
    from .closure import close
    from .feasibility import assess, transfer_asymmetry_check
    from .mission import geo_mission

    if mission is not None and mission.key != "geo":
        raise ValueError(
            f"sweep_t_transfer_max needs the GEO mission; got {mission.key!r}"
        )
    keys = ARCHITECTURE_ORDER if arch_keys is None else arch_keys

    records = {k: [] for k in keys}
    ep_keys = (
        tuple(k for k in keys if get_architecture(k).has_ep) if ep_transfer else ()
    )
    ep_records = {k: [] for k in ep_keys}
    asymmetry = {}
    for t_yr in t_set:
        m_t = geo_mission(t_transfer_max_yr=t_yr)
        for k in keys:
            arch = get_architecture(k)
            alloc = nominal_allocation(k, m_t)
            r = close(k, m_t, alloc)
            if r.diverged or not r.converged:
                records[k].append(
                    dict(
                        t_yr=t_yr,
                        wet=nan,
                        converged=r.converged,
                        diverged=r.diverged,
                        excluded=True,
                        feasible=False,
                        reasons=[
                            "closure diverged (no fixed point)"
                            if r.diverged
                            else "closure not converged"
                        ],
                        n_ep=None,
                        p_dc_w=nan,
                        eps_inc=nan,
                        electronics=nan,
                        t_comm_s=nan,
                    )
                )
                continue
            rep = assess(k, m_t, alloc, r)
            records[k].append(
                dict(
                    t_yr=t_yr,
                    wet=r.wet,
                    converged=True,
                    diverged=False,
                    excluded=rep.excluded,
                    feasible=not rep.excluded,
                    reasons=rep.reasons,
                    n_ep=r.n_units["ep"],
                    p_dc_w=r.p_dc_peak,
                    eps_inc=r.eps_inc,
                    electronics=r.electronics,
                    t_comm_s=_transfer_commission_s(arch, m_t, r),
                )
            )

        for k in ep_keys:
            arch = get_architecture(k)
            alloc, dg = ep_transfer_allocation(k, m_t)
            r = close(k, m_t, alloc)
            base = dict(
                t_yr=t_yr,
                s_ep=dg["s_ep"],
                n_max=dg["n_max"],
                forced_all_ep=dg["forced_all_ep"],
                no_feasible_s=dg["no_feasible_s"],
                dv_ep=alloc["transfer"]["ep"],
                dv_chem=alloc["transfer"]["chem"],
            )
            if r.diverged or not r.converged:
                if r.diverged and dg["no_feasible_s"]:
                    reason = (
                        "no power-feasible EP share exists at this "
                        "window; the all-chemical (s = 0) fallback "
                        "has no fixed point (160 s transfer)"
                    )
                elif r.diverged:
                    reason = (
                        "closure diverged (no fixed point: "
                        "EPS-feedback slope > 1 at this window)"
                    )
                else:
                    reason = "closure not converged"
                ep_records[k].append(
                    dict(
                        base,
                        wet=nan,
                        converged=r.converged,
                        diverged=r.diverged,
                        excluded=True,
                        feasible=False,
                        power_ok=False,
                        reasons=[reason],
                        n_ep=None,
                        p_dc_w=nan,
                        eps_inc=nan,
                        tr_prop_ep=nan,
                        tr_prop_chem=nan,
                        t_comm_s=nan,
                    )
                )
                continue
            rep = assess(k, m_t, alloc, r)
            ep_records[k].append(
                dict(
                    base,
                    wet=r.wet,
                    converged=True,
                    diverged=False,
                    excluded=rep.excluded,
                    feasible=not rep.excluded,
                    power_ok=rep.power_ok,
                    reasons=rep.reasons,
                    n_ep=r.n_units["ep"],
                    p_dc_w=r.p_dc_peak,
                    eps_inc=r.eps_inc,
                    tr_prop_ep=r.per_class_prop["transfer"]["ep"],
                    tr_prop_chem=r.per_class_prop["transfer"]["chem"],
                    t_comm_s=_transfer_commission_s(arch, m_t, r),
                )
            )

        asymmetry[t_yr] = transfer_asymmetry_check(t_max_yr=t_yr)

    return {
        "t_set_yr": tuple(t_set),
        "records": records,
        "ep_transfer": ep_records,
        "asymmetry": asymmetry,
    }


def _tornado_bars():
    """Sensitivity cases: each 'oat' register row at its bracket ends, and
    three bands on the electrodeless unit: exhaust velocity +/-15 % (Isp
    and mass flow, thrust and power held), thrust +/-10 % (thrust and mass
    flow) and absorbed power +/-20 %.
    """
    from .parameters import REGISTER, value

    bars = []
    for key, p in REGISTER.items():
        if p.sweep == "oat":
            bars.append(
                dict(
                    name=key,
                    group="register",
                    flag=p.flag,
                    lo={key: p.range[0]},
                    hi={key: p.range[1]},
                    lo_label=f"{p.range[0]:g}",
                    hi_label=f"{p.range[1]:g}",
                )
            )
    isp0 = value("isp_ecr")
    mdot0 = value("mdot_unit_ecr")
    t0 = value("t_unit_ecr")
    p0 = value("p_abs_ecr")
    bars.append(
        dict(
            name="thruster_R_Isp_ue",
            group="thruster",
            flag="T",
            lo={"isp_ecr": 0.85 * isp0, "mdot_unit_ecr": mdot0 / 0.85},
            hi={"isp_ecr": 1.15 * isp0, "mdot_unit_ecr": mdot0 / 1.15},
            lo_label="-15% (probable edge)",
            hi_label="+15%",
        )
    )
    bars.append(
        dict(
            name="thruster_R_T",
            group="thruster",
            flag="T",
            lo={"t_unit_ecr": 0.9 * t0, "mdot_unit_ecr": 0.9 * mdot0},
            hi={"t_unit_ecr": 1.1 * t0, "mdot_unit_ecr": 1.1 * mdot0},
            lo_label="-10% (probable edge)",
            hi_label="+10%",
        )
    )
    bars.append(
        dict(
            name="thruster_R_Pin",
            group="thruster",
            flag="T",
            lo={"p_abs_ecr": 0.8 * p0},
            hi={"p_abs_ecr": 1.2 * p0},
            lo_label="-20%",
            hi_label="+20%",
        )
    )
    return bars


def tornado(
    arch_keys=("multimode", "single_elec", "hybrid", "mono_multimode"), missions=None
):
    """One-at-a-time sensitivity of the wet masses and of the differences
    multimode - single-electric (LEO) and multimode - hybrid (GEO). An
    input the model cannot evaluate at a bracket end is recorded with its
    error.

    Returns {'baseline': {mission: metrics}, 'bars': [{name, group, flag,
    edges: {'lo'/'hi': {value_label, error, leo: {...}, geo: {...}}}}]};
    metrics wet_mm, wet_se, wet_hy, wet_mmm, d_mm_se, d_mm_hy, d_mmm_mm.
    """
    from math import nan
    from .allocation import nominal_allocation
    from .architectures import invalidate_cache
    from .closure import close
    from .mission import leo_mission, geo_mission
    from .parameters import override

    def _metrics():
        out = {}
        for build, mkey in ((leo_mission, "leo"), (geo_mission, "geo")):
            m = build()
            wets = {}
            for k in arch_keys:
                r = close(k, m, nominal_allocation(k, m))
                wets[k] = r.wet if r.converged else nan
            out[mkey] = {
                "wet_mm": wets.get("multimode", nan),
                "wet_se": wets.get("single_elec", nan),
                "wet_hy": wets.get("hybrid", nan),
                "wet_mmm": wets.get("mono_multimode", nan),
                "d_mm_se": wets.get("multimode", nan) - wets.get("single_elec", nan),
                "d_mm_hy": wets.get("multimode", nan) - wets.get("hybrid", nan),
                "d_mmm_mm": wets.get("mono_multimode", nan)
                - wets.get("multimode", nan),
            }
        return out

    baseline = _metrics()

    bars = []
    for bar in _tornado_bars():
        edges = {}
        for edge in ("lo", "hi"):
            try:
                with override(**bar[edge]):
                    invalidate_cache()
                    met = _metrics()
                err = None
            except Exception as e:  # input not evaluable
                met = None
                err = f"{type(e).__name__}: {e}"
            finally:
                invalidate_cache()
            edges[edge] = dict(
                value_label=bar[f"{edge}_label"],
                error=err,
                **(met if met is not None else {}),
            )
        bars.append(
            dict(name=bar["name"], group=bar["group"], flag=bar["flag"], edges=edges)
        )
    return {"baseline": baseline, "bars": bars}


# Planning efficiency and array power loss, GEO mission

ETA_PLANNING_SET = (0.75, 0.70, 0.65, 0.62)
ARR_DEGR_REPORT_SET = (0.0, 0.24)
STRATEGIES = ("time-optimal", "maximum-electric")

STRATEGY_LABEL = {
    "time-optimal": "time-optimal transfer: chemical wherever a chemical "
    "mode exists (nominal)",
    "maximum-electric": "maximum-electric transfer: largest power-feasible electric share",
}


def _eta_planning_alloc(arch_key, mission, strategy):
    """Allocation of one architecture under one transfer strategy."""
    from .allocation import nominal_allocation
    from .architectures import get_architecture

    arch = get_architecture(arch_key)
    if strategy == "maximum-electric" and arch.has_ep:
        alloc, _dg = ep_transfer_allocation(arch_key, mission)
        return alloc
    return nominal_allocation(arch_key, mission)


def _eta_planning_point(eta, arr_degr, strategies=STRATEGIES, arch_keys=None):
    """Every architecture of the GEO mission at one planning efficiency and
    array power loss.
    """
    from .allocation import transfer_ep_share
    from .architectures import ARCHITECTURE_ORDER, invalidate_cache
    from .closure import close
    from .feasibility import assess
    from .mission import geo_mission
    from .parameters import override, value

    keys = ARCHITECTURE_ORDER if arch_keys is None else arch_keys
    dv_chem = value("dv_transfer_chem")
    dv_ep = dv_chem / eta
    out = {}
    with override(eta_planning_transfer=eta, arr_degr_transfer=arr_degr):
        invalidate_cache()
        mission = geo_mission()
        for strategy in strategies:
            rows = []
            for key in keys:
                alloc = _eta_planning_alloc(key, mission, strategy)
                res = close(key, mission, alloc)
                ok = res.converged and not res.diverged
                rep = assess(key, mission, alloc, res) if ok else None
                rows.append(
                    dict(
                        strategy=strategy,
                        eta=eta,
                        arr_degr=arr_degr,
                        arch=key,
                        dv_ep=dv_ep,
                        loss_ratio=dv_ep / dv_chem,
                        s_ep=transfer_ep_share(mission, alloc),
                        d_eff=res.arr_degr,
                        p_dc_w=res.p_dc_peak if ok else float("nan"),
                        n_ep=res.n_units["ep"] if ok else None,
                        eps_inc=res.eps_inc if ok else float("nan"),
                        wet=res.wet if ok else float("nan"),
                        env_frac=(res.wet / mission.bus_mass) if ok else float("nan"),
                        envelope_cap=(rep.envelope_cap if rep else float("nan")),
                        power_ok=rep.power_ok if rep else None,
                        envelope_ok=rep.envelope_ok if rep else None,
                        excluded=rep.excluded if rep else None,
                        converged=res.converged,
                        diverged=res.diverged,
                    )
                )
            out[strategy] = rows
        invalidate_cache()
    invalidate_cache()
    return mission, out


def sweep_eta_planning(
    eta_set=None, degr_set=None, strategies=STRATEGIES, arch_keys=None
):
    """GEO designs over planning efficiency and array power loss, for both
    transfer strategies.

    Returns {'grid': {(eta, d): {strategy: [row, ...]}}, 'combos',
    'survivors': {strategy: {(eta, d): tuple}}, 'first_change',
    'mission', 'eta_set', 'degr_set', 'strategies'}.
    """
    eta_set = ETA_PLANNING_SET if eta_set is None else tuple(eta_set)
    degr_set = ARR_DEGR_REPORT_SET if degr_set is None else tuple(degr_set)
    combos = [(e, d) for e in eta_set for d in degr_set]

    grid = {}
    mission = None
    for eta, d in combos:
        mission, out = _eta_planning_point(eta, d, strategies, arch_keys)
        grid[(eta, d)] = out

    survivors, first_change = {}, {}
    for strategy in strategies:
        base = None
        survivors[strategy] = {}
        for combo in combos:
            surv = tuple(
                sorted(
                    r["arch"] for r in grid[combo][strategy] if r["envelope_ok"] is True
                )
            )
            survivors[strategy][combo] = surv
            if base is None:
                base = surv
            elif surv != base and strategy not in first_change:
                first_change[strategy] = (combo[0], combo[1], surv)
    return dict(
        grid=grid,
        combos=combos,
        survivors=survivors,
        first_change=first_change,
        mission=mission,
        eta_set=eta_set,
        degr_set=degr_set,
        strategies=strategies,
    )


def eta_planning_crossover(arch_key, strategy, arr_degr, lo=0.50, hi=0.80, tol=1e-4):
    """Lowest planning efficiency at which arch_key stays within the mass
    limit, by bisection; None if [lo, hi] does not straddle the boundary.
    """

    def ok_at(eta):
        _m, out = _eta_planning_point(
            eta, arr_degr, strategies=(strategy,), arch_keys=(arch_key,)
        )
        return out[strategy][0]["envelope_ok"] is True

    if not ok_at(hi):
        return None
    if ok_at(lo):
        return None
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if ok_at(mid):
            hi = mid
        else:
            lo = mid
    return hi  # the passing side of the boundary
