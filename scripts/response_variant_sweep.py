"""LEO response variant: propulsion wet mass against the total response delta-v."""

# Each grid point sets the response burn and the return burn to the same
# per-event delta-v (total / n_ca_events) at a 1 d warning time and closes
# every architecture at its nominal allocation (and the two multimode
# designs also at f_EP = 1). Writes, to output/response_variant/:
#
#     response_variant.csv   one row per grid point and architecture
#     response_variant.json  single-electric limit and lightest-design changes

import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from multimode.config import PATHS  # noqa: E402
from multimode.spacecraft import closure as C  # noqa: E402
from multimode.spacecraft.allocation import (  # noqa: E402
    derived_f_ep,
    nominal_allocation,
    sweep_allocation,
)
from multimode.spacecraft.architectures import (  # noqa: E402
    ARCHITECTURE_ORDER,
    get_architecture,
    invalidate_cache,
)
from multimode.spacecraft.feasibility import assess  # noqa: E402
from multimode.spacecraft.mission import leo_mission  # noqa: E402
from multimode.spacecraft.parameters import override, value  # noqa: E402

# Total response delta-v grid [m/s]
TOTALS = (1.5, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0)
N_EV = int(value("n_ca_events"))
T_WARN = 86400.0
MM = ("multimode", "mono_multimode")
CSV_ORDER = ("multimode", "mono_multimode", "single_elec", "hybrid", "single_chem")
OUT_DIR = PATHS["output"] / "response_variant"

COLUMNS = [
    "total_response_dv_mps",
    "per_event_mps",
    "response_total_mps",
    "return_total_mps",
    "mission_total_impulsive_mps",
    "architecture",
    "wet_kg",
    "mass_reported",
    "wet_kg_as_closed",
    "converged",
    "feasible",
    "failing_check",
    "reasons",
    "n_ep",
    "n_chem",
    "f_ep",
    "response_mode",
    "return_mode",
    "response_burn_h",
    "response_limit_h",
    "response_charge_factor",
    "drift_ok",
    "prop_response_kg",
    "prop_return_kg",
    "m_at_response_kg",
    "return_ep_wall_per_event_d",
    "event_interval_d",
    "sk_plus_return_ep_wall_share",
    "f1_wet_kg",
    "f1_mass_reported",
    "f1_wet_kg_as_closed",
    "f1_feasible",
    "f1_failing_check",
    "f1_n_ep",
    "f1_response_mode",
    "f1_return_mode",
    "f1_response_burn_h",
    "f1_drift_ok",
]


def check_label(reasons):
    """Short name of each failing check."""
    out = []
    for r in reasons:
        if r.startswith("ca/") and "drift unreachable" in r:
            out.append("CA drift (no root)")
        elif r.startswith("ca/"):
            out.append("CA duty limit")
        elif r.startswith("power"):
            out.append("power ceiling")
        elif r.startswith("envelope"):
            out.append("envelope")
        else:
            out.append(r.split(":")[0])
    return out


def mode_of(entry):
    modes = [m for m in ("chem", "ep") if entry[m] > 0.0]
    return "+".join(modes) if modes else "none"


def evaluate(k, m, al):
    """Closure and assessment of one design at one grid point."""
    arch = get_architecture(k)
    r = C.close(k, m, al)
    rec = dict(
        converged=r.converged,
        f_ep=derived_f_ep(m, al),
        response_mode=mode_of(al["ca"]),
        return_mode=mode_of(al["ca_restore"]),
    )
    if not r.converged:
        return rec
    rep = assess(k, m, al, r)
    v = [x for x in rep.verdicts if x.class_key == "ca"]
    v = v[0] if v else None
    pc = r.per_class_prop
    rec.update(
        wet=r.wet,
        n_ep=r.n_units["ep"],
        n_chem=r.n_units["chem"],
        feasible=not rep.excluded,
        reasons=rep.reasons,
        failing=check_label(rep.reasons),
        ca_charge=r.ca_charge,
        drift_ok=(r.ca_charge is None or r.ca_charge["drift_ok"]),
        ca_burn_h=(v.t_burn_s / 3600 if v else None),
        ca_limit_h=(v.t_limit_s / 3600 if v else None),
        prop_ca=sum(pc["ca"].values()),
        prop_restore=sum(pc["ca_restore"].values()),
        # spacecraft mass when the response burn starts (SK, then CA)
        m_at_ca=m.bus_mass + r.wet - sum(pc["sk"].values()),
    )
    # Electric wall-clock time of the return burn against the mean
    # interval between events
    if arch.has_ep and r.n_units["ep"] > 0:
        mdot = r.n_units["ep"] * arch.ep_unit.mdot
        w = m.get_class("ca_restore").w_duty
        t_on_rs = pc["ca_restore"]["ep"] / mdot
        t_on_sk = pc["sk"]["ep"] / mdot
        rec.update(
            restore_ep_wall_per_event_d=t_on_rs / w / N_EV / 86400,
            event_interval_d=m.duration_s / N_EV / 86400,
            sk_plus_restore_wall_share=(t_on_sk + t_on_rs) / w / m.duration_s,
        )
    return rec


def point(per_event, archs=ARCHITECTURE_ORDER, f1=MM):
    """Grid point at this per-event delta-v: nominal allocation for `archs`,
    f_EP = 1 for `f1`."""
    try:
        with override(dv_ca_budget_1d=per_event, dv_ca_restore=per_event):
            invalidate_cache()
            m = leo_mission(t_warn_s=T_WARN)
            meta = dict(
                per_event=per_event,
                ca_total=m.get_class("ca").dv_impulsive,
                restore_total=m.get_class("ca_restore").dv_impulsive,
                mission_total=m.total_dv_impulsive,
            )
            nom = {k: evaluate(k, m, nominal_allocation(k, m)) for k in archs}
            ae = {k: evaluate(k, m, sweep_allocation(k, m, 1.0)) for k in f1}
    finally:
        invalidate_cache()
    return meta, nom, ae


def bisect(pred, lo, hi, tol=1e-12):
    """pred(lo) False, pred(hi) True; returns (lo, hi)."""
    assert not pred(lo) and pred(hi), (lo, hi)
    while hi - lo > tol * max(1.0, hi):
        mid = 0.5 * (lo + hi)
        if pred(mid):
            hi = mid
        else:
            lo = mid
    return lo, hi


def transitions(pairs):
    """Indices i where the value changes between pairs[i] and pairs[i+1]."""
    return [i for i in range(len(pairs) - 1) if pairs[i][1] != pairs[i + 1][1]]


def reported(v):
    return v.get("converged") and v.get("drift_ok", True)


def single_elec_limit():
    """Total response [m/s] at which single-electric is first excluded."""

    def feasible(b):
        return point(b, archs=("single_elec",), f1=())[1]["single_elec"]["feasible"]

    xs = [10 ** (-2 + 3.2 * i / 400) for i in range(401)]
    sc = [(x, feasible(x)) for x in xs]
    i = transitions(sc)[0]
    _lo, hi = bisect(lambda b: not feasible(b), sc[i][0], sc[i + 1][0])
    return hi * N_EV


def lightest_changes():
    """Where the lightest feasible design (nominal allocation) changes
    along the total-response axis: 0.25 m/s scan, each change bisected."""

    def lightest_at(tot):
        n = point(tot / N_EV, f1=())[1]
        cand = {
            k: v["wet"]
            for k, v in n.items()
            if v.get("feasible") and v.get("drift_ok", True)
        }
        return min(cand, key=cand.get)

    ts = [1.5 + 0.25 * i for i in range(int((200 - 1.5) / 0.25) + 1)]
    sc = [(t, lightest_at(t)) for t in ts]
    out = []
    for i in transitions(sc):
        after = sc[i + 1][1]
        lo, hi = bisect(lambda t: lightest_at(t) == after, sc[i][0], sc[i + 1][0])
        out.append(dict(total_lo=lo, total_hi=hi, before=sc[i][1], after=after))
    return out


def write_csv(path):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(COLUMNS)
        for tot in TOTALS:
            mt, nom, ae = point(tot / N_EV)
            for k in CSV_ORDER:
                v = nom[k]
                a = ae.get(k, {})
                cf = (
                    (v.get("ca_charge") or {}).get("factor", 1.0)
                    if v.get("converged")
                    else ""
                )
                row = [
                    tot,
                    mt["per_event"],
                    mt["ca_total"],
                    mt["restore_total"],
                    mt["mission_total"],
                    k,
                    v["wet"] if reported(v) else "",
                    bool(reported(v)),
                    v.get("wet", ""),
                    v["converged"],
                    v.get("feasible", ""),
                    "; ".join(v.get("failing", [])),
                    " | ".join(v.get("reasons", [])),
                    v.get("n_ep", ""),
                    v.get("n_chem", ""),
                    v["f_ep"],
                    v["response_mode"],
                    v["return_mode"],
                    v.get("ca_burn_h", ""),
                    v.get("ca_limit_h", ""),
                    cf,
                    v.get("drift_ok", ""),
                    v.get("prop_ca", ""),
                    v.get("prop_restore", ""),
                    v.get("m_at_ca", ""),
                    v.get("restore_ep_wall_per_event_d", ""),
                    v.get("event_interval_d", ""),
                    v.get("sk_plus_restore_wall_share", ""),
                ]
                if a:
                    row += [
                        a["wet"] if reported(a) else "",
                        bool(reported(a)),
                        a.get("wet", ""),
                        a.get("feasible", ""),
                        "; ".join(a.get("failing", [])),
                        a.get("n_ep", ""),
                        a["response_mode"],
                        a["return_mode"],
                        a.get("ca_burn_h", ""),
                        a.get("drift_ok", ""),
                    ]
                else:
                    row += [""] * 10
                w.writerow(row)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    write_csv(OUT_DIR / "response_variant.csv")
    summary = dict(
        single_elec_first_excluded_total=single_elec_limit(),
        lightest_changes=lightest_changes(),
    )
    (OUT_DIR / "response_variant.json").write_text(
        json.dumps(summary, indent=1), encoding="utf-8"
    )
    print(f"written {OUT_DIR}")


if __name__ == "__main__":
    main()
