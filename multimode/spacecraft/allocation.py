"""
Allocation of each manoeuvre class between the chemical and electric modes,
as class key -> {'chem': dv, 'ep': dv}; the electric delta-v includes the
class loss ratio.
"""

# nominal_allocation gives the default table below; sweep_allocation
# assigns a fraction f_EP of the impulsive budget to the electric mode in
# order of schedule tolerance (disposal, station keeping, transfer,
# collision avoidance). Station keeping and transfer can be split; collision
# avoidance and disposal go to one mode whole.
#
# Default table:
#
#     Class            Multimode  Mono-mm  Hybrid  Single-chem  Single-elec
#     SK               EP         EP       EP      chem         EP
#     CA               chem       chem     chem    chem         EP
#     CA return        EP         EP       EP      chem         EP
#     Transfer         chem       chem     chem    chem         EP
#     Disposal         EP         EP       EP      chem         EP
#
# The transfer is flown chemically by every architecture with a chemical
# mode (time-optimal), and the return burn electrically by every
# architecture with an electric mode.

from .architectures import get_architecture

# Default mode per class and architecture
_NOMINAL_TABLE = {
    "leo": {
        "single_elec": {"sk": "ep", "ca": "ep", "ca_restore": "ep", "eol": "ep"},
        "multimode": {"sk": "ep", "ca": "chem", "ca_restore": "ep", "eol": "ep"},
        "hybrid": {"sk": "ep", "ca": "chem", "ca_restore": "ep", "eol": "ep"},
        "single_chem": {
            "sk": "chem",
            "ca": "chem",
            "ca_restore": "chem",
            "eol": "chem",
        },
        "mono_multimode": {"sk": "ep", "ca": "chem", "ca_restore": "ep", "eol": "ep"},
    },
    "geo": {
        "single_elec": {"transfer": "ep", "sk": "ep", "eol": "ep"},
        "multimode": {"transfer": "chem", "sk": "ep", "eol": "ep"},
        "hybrid": {"transfer": "chem", "sk": "ep", "eol": "ep"},
        "single_chem": {"transfer": "chem", "sk": "chem", "eol": "chem"},
        "mono_multimode": {"transfer": "chem", "sk": "ep", "eol": "ep"},
    },
}


def _entry(cls, share):
    """{'chem': dv, 'ep': dv including losses} for a class at an electric
    share of its impulsive budget."""
    dv_ep_imp = share * cls.dv_impulsive
    return {"chem": cls.dv_impulsive - dv_ep_imp, "ep": dv_ep_imp * cls.ep_loss_ratio}


def transfer_ep_share(mission, allocation):
    """Electric share s of the transfer, on the impulsive basis:
    s = (ep / loss ratio) / dv_impulsive. 0 without a transfer class."""
    for cls in mission.classes:
        if cls.kind != "transfer":
            continue
        entry = allocation[cls.key]
        if cls.dv_impulsive <= 0.0 or cls.ep_loss_ratio <= 0.0:
            return 0.0
        return (entry["ep"] / cls.ep_loss_ratio) / cls.dv_impulsive
    return 0.0


def nominal_allocation(arch_key, mission):
    """Default allocation of an architecture (table above)."""
    arch = get_architecture(arch_key)
    table = _NOMINAL_TABLE[mission.key][arch_key]
    alloc = {}
    for cls in mission.classes:
        mode = table[cls.key]
        if cls.kind == "transfer":
            mode = "chem" if arch.has_chem else "ep"
            assert mode == table[cls.key], (
                f"_NOMINAL_TABLE transfer entry for {arch_key!r} does not "
                f"match the has-chem rule"
            )
        elif cls.kind == "ca_restore":
            mode = "ep" if arch.has_ep else "chem"
            assert mode == table[cls.key], (
                f"_NOMINAL_TABLE ca_restore entry for {arch_key!r} does "
                f"not match the has-EP rule"
            )
        alloc[cls.key] = _entry(cls, 1.0 if mode == "ep" else 0.0)
    return alloc


def sweep_allocation(arch_key, mission, f_ep, tol=1e-12):
    """Allocation with a fraction f_ep of the impulsive budget electric,
    assigned in order of schedule tolerance. Single-mode architectures
    admit only f_ep = 0 or 1."""
    arch = get_architecture(arch_key)
    if not 0.0 <= f_ep <= 1.0:
        raise ValueError(f"f_ep = {f_ep} outside [0, 1]")
    if not arch.has_ep and f_ep > tol:
        raise ValueError(f"{arch_key} has no EP mode; f_ep must be 0")
    if not arch.has_chem and f_ep < 1.0 - tol:
        raise ValueError(f"{arch_key} has no chemical mode; f_ep must be 1")

    remaining = f_ep * mission.total_dv_impulsive
    alloc = {}
    for cls in sorted(mission.classes, key=lambda c: c.assignment_rank):
        if cls.divisible:
            take = min(remaining, cls.dv_impulsive)
        else:
            take = cls.dv_impulsive if remaining >= cls.dv_impulsive - tol else 0.0
        remaining -= take
        alloc[cls.key] = _entry(cls, take / cls.dv_impulsive)
    return {c.key: alloc[c.key] for c in mission.classes}


def derived_f_ep(mission, allocation):
    """Electric fraction of the impulsive budget of an allocation."""
    ep_imp = 0.0
    for cls in mission.classes:
        ep_imp += allocation[cls.key]["ep"] / cls.ep_loss_ratio
    return ep_imp / mission.total_dv_impulsive
