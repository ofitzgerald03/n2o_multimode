"""
Parameter register of the spacecraft-level model: every input value, with
units, sensitivity range, provenance flag, source and sweep role. Modules
read values through `value(key)`; `export_register_csv` writes the table.
"""

# Flags: P primary source; S secondary source; E estimated from comparable
# hardware or a stated construction; D declared for this study; V derived
# from other entries or geometry; T composed by the thruster model.
# Sweep roles: 'fixed', 'oat' (one-at-a-time sensitivity), 'axis' (swept
# axis).

import csv
from contextlib import contextmanager
from dataclasses import dataclass, asdict
from typing import Optional, Tuple, Union

from ..config import PATHS

# Seconds per year
YEAR_S = 3.156e7

Value = Union[float, int, str]


@dataclass(frozen=True)
class Parameter:
    """One register row; range is the (lo, hi) sensitivity bracket or
    None."""

    key: str
    value: Value
    units: str
    range: Optional[Tuple[float, float]]
    flag: str
    source: str
    sweep: str
    note: str = ""


REGISTER: dict = {}


def _add(key, value, units, rng, flag, source, sweep, note=""):
    if key in REGISTER:
        raise ValueError(f"Duplicate register key: {key!r}")
    REGISTER[key] = Parameter(key, value, units, rng, flag, source, sweep, note)


def value(key):
    """Value of a register row."""
    return REGISTER[key].value


@contextmanager
def override(**overrides):
    """Temporarily replace the values of register rows.

    Objects built from register values and cached (the architectures) do
    not see an override: pair it with architectures.invalidate_cache() on
    entry and on exit.
    """
    import dataclasses

    saved = {}
    try:
        for key, val in overrides.items():
            if key not in REGISTER:
                raise KeyError(f"override: no register row {key!r}")
            saved[key] = REGISTER[key]
            REGISTER[key] = dataclasses.replace(saved[key], value=val)
        yield
    finally:
        REGISTER.update(saved)


FLAG_LEGEND = (
    ("P", "primary: from a datasheet or primary literature at the condition stated"),
    ("S", "secondary: from a source reporting another's work"),
    (
        "E",
        "estimated: constructed from comparable hardware or by a stated construction",
    ),
    (
        "D",
        "declared: selected for this study within a range the "
        "literature supports, and swept across it",
    ),
    ("V", "derived: by calculation from other entries or from geometry"),
    ("T", "composed by the thruster model at the anchor (eta_u parity)"),
)


def export_register_csv(path=None, legend=True):
    """Write the register as CSV (default output/parameter_register.csv);
    legend appends the flag legend as comment rows."""
    if path is None:
        out_dir = PATHS["output"]
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / "parameter_register.csv"
    cols = ("key", "value", "units", "range", "flag", "source", "sweep", "note")
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for p in REGISTER.values():
            d = asdict(p)
            d["range"] = "" if p.range is None else f"{p.range[0]}..{p.range[1]}"
            w.writerow([d[c] for c in cols])
        if legend:
            w.writerow([])
            w.writerow(["# flag legend", "", "", "", "flag", "meaning", "", ""])
            for flag, meaning in FLAG_LEGEND:
                w.writerow(["#", "", "", "", flag, meaning, "", ""])
    return path


# Tanks and storage

_add(
    "t_max_storage",
    30.0,
    "degC",
    (30.0, 35.0),
    "P",
    "B20 and SatDrive qualified storage range, upper limit 30 C (datasheet)",
    "fixed",
    "Storage pressures and densities below are evaluated at this "
    "temperature; N2O critical temperature 36.4 C.",
)

_add(
    "p_meop_n2o",
    63.075e5,
    "Pa",
    None,
    "P",
    "N2O saturation pressure at 30 C, 63.075 bar: NIST Chemistry WebBook "
    "(SRD 69), saturation table",
    "fixed",
)
_add(
    "rho_min_n2o",
    688.11,
    "kg/m3",
    None,
    "P",
    "N2O saturated-liquid density at 30 C, 688.11 kg/m3: NIST Chemistry "
    "WebBook (SRD 69), saturation table",
    "fixed",
)
_add(
    "p_meop_propene",
    13.050e5,
    "Pa",
    None,
    "P",
    "Propene saturation pressure at 30 C, 13.050 bar: NIST Chemistry "
    "WebBook (SRD 69), saturation table",
    "fixed",
    "Fuel is propene, per the B20 datasheet.",
)
_add(
    "rho_min_propene",
    497.45,
    "kg/m3",
    None,
    "P",
    "Propene saturated-liquid density at 30 C, 497.45 kg/m3: NIST "
    "Chemistry WebBook (SRD 69), saturation table",
    "fixed",
)
_add(
    "p_meop_xe", 150e5, "Pa", None, "S", "Xe supercritical storage at 150 bar", "fixed"
)
_add("rho_min_xe", 1600.0, "kg/m3", None, "S", "Xe density at about 150 bar", "fixed")

_add("ullage_factor", 1.10, "-", None, "S", "Falcone et al. (10 % ullage)", "fixed")
_add(
    "residual_factor",
    1.05,
    "-",
    None,
    "S",
    "Falcone et al. (5 % unusable propellant)",
    "fixed",
)

_add(
    "pf_metal",
    5080.0,
    "m",
    (3500.0, 5080.0),
    "V",
    "Berg and Rovey 2014, Sec. C.1: titanium tank sizing parameter "
    "6350 m on a burst basis, divided by their burst factor 1.25 to the "
    "MEOP basis: 5080 m. Lower edge: factor-six titanium relation of "
    "Humble, used by Donius and Rovey 2011, about 3500 m",
    "oat",
    "Tank performance factor PF = P V / (m_tank g) of metal (Type 1) tanks.",
)
_add(
    "alpha_xe_copv",
    0.045,
    "-",
    None,
    "P",
    "Dankanich and Woodcock: Xe COPV tank fraction 4.5 %",
    "fixed",
    "Sets the COPV performance factor.",
)
_add("m_min_tank", 0.5, "kg", (0.3, 1.0), "D", "Minimum tank mass", "oat")
_add(
    "m_fixed_tank",
    0.3,
    "kg",
    (0.15, 0.45),
    "E",
    "Propellant management device, bosses and mounts per tank",
    "oat",
    "Items not carried by the performance factor.",
)
_add(
    "chem_tank_technology",
    "metal",
    "-",
    None,
    "P",
    "SatDrive datasheet: Type 1 all-metal flown configurations",
    "fixed",
    "Tank type of N2O and propene (metal or copv); Xe is always COPV.",
)

# Thruster heads and units

_add(
    "m_head_ecr",
    1.0,
    "kg",
    (0.5, 1.5),
    "E",
    "N2O electrodeless head. Comparators: Berg and Rovey 2014 Table 4, "
    "helicon-type thruster 0.5 kg; Donius and Rovey 2011, SPT-50 0.8 kg; "
    "Falcone et al. 2025 Table 2, electric thruster 1.20 kg",
    "oat",
)
_add(
    "m_head_b20",
    0.695,
    "kg",
    None,
    "P",
    "Dawn Aerospace B20 datasheet: dry mass 695 g",
    "fixed",
)
_add(
    "m_head_hall",
    1.10,
    "kg",
    None,
    "S",
    "BHT-200 class: 0.98 kg thruster and 0.12 kg cathode",
    "fixed",
)
_add(
    "m_head_mono",
    0.5,
    "kg",
    (0.25, 0.75),
    "E",
    "N2O catalytic monopropellant head. Comparators: Berg and Rovey 2014 "
    "Table 4, monopropellant thruster 0.5 kg; Donius and Rovey 2011, "
    "surveyed mean 0.21 kg",
    "oat",
)

_add(
    "isp_chem",
    271.0,
    "s",
    None,
    "P",
    "Dawn Aerospace B20 datasheet: delivered Isp of at least 271 s",
    "fixed",
)
_add(
    "isp_mono",
    160.0,
    "s",
    (130.0, 182.0),
    "S",
    "N2O monopropellant delivered Isp",
    "oat",
)
_add(
    "p_heater_mono",
    25.0,
    "W",
    (10.0, 50.0),
    "S",
    "Catalyst-heater power of the N2O monopropellant head: Zakirov et al. "
    "2000, 24 W for catalyst activation",
    "oat",
)
_add("isp_hall", 1390.0, "s", None, "S", "BHT-200 class Isp", "fixed")
_add("isp_ecr", 1808.0, "s", None, "T", "Thruster model, N2O at eta_u = 0.61", "fixed")

_add("of_ratio", 7.0, "-", None, "S", "N2O/hydrocarbon O/F of 7", "fixed")

_add(
    "t_unit_ecr",
    0.98e-3,
    "N",
    None,
    "T",
    "Thruster model, N2O at eta_u = 0.61",
    "fixed",
)
_add(
    "p_abs_ecr",
    61.3,
    "W",
    None,
    "T",
    "Thruster model, N2O at eta_u = 0.61: absorbed microwave power",
    "fixed",
    "DC power follows through eta_dcrf.",
)
_add(
    "mdot_unit_ecr",
    0.0553e-6,
    "kg/s",
    None,
    "T",
    "Thruster model, N2O at eta_u = 0.61",
    "fixed",
)

_add("t_unit_hall", 13e-3, "N", None, "S", "BHT-200 class thrust", "fixed")
_add("p_disch_hall", 200.0, "W", None, "S", "BHT-200 class discharge power", "fixed")

_add(
    "t_b20_nominal",
    12.3,
    "N",
    (6.46, 18.11),
    "V",
    "Midpoint of the B20 datasheet thrust range 6.46-18.11 N",
    "fixed",
    "Chemical burn times only.",
)
_add(
    "t_mono_nominal",
    10.0,
    "N",
    None,
    "D",
    "Monopropellant thrust in the B20 thrust class",
    "fixed",
    "Chemical burn times only.",
)

# Power chain

_add(
    "eta_dcrf",
    0.65,
    "-",
    (0.5, 0.7),
    "S",
    "DC-to-RF conversion efficiency at 2.45 GHz",
    "oat",
)
_add("eta_ppu", 0.93, "-", (0.92, 0.95), "S", "Hall PPU efficiency", "fixed")

_add(
    "sm_ppu",
    6.2,
    "kg/kW",
    (4.0, 8.0),
    "S",
    "Donius and Rovey 2011, Sec. A: PPU specific power 161 W/kg "
    "(6.2 kg/kW). Upper edge: Kluever 2012, 7.2 kg/kW for a complete Hall "
    "subsystem",
    "oat",
)
_add(
    "sm_rf",
    12.4,
    "kg/kW",
    (6.2, 18.6),
    "E",
    "RF-chain specific mass, twice sm_ppu; bracket one to three times "
    "sm_ppu. Comparators: Donius and Rovey 2011, 6.2 kg/kW; Kluever 2012, "
    "7.2 kg/kW",
    "oat",
)
_add(
    "sm_eps_inc",
    12.5,
    "kg/kW",
    (10.0, 15.0),
    "E",
    "Incremental power-system specific mass: solar array at about 85 W/kg "
    "(11.8 kg/kW) plus power management and thermal",
    "oat",
    "Electric thrusting in sunlight only; no battery increment.",
)

_add(
    "arr_degr_transfer",
    0.18,
    "-",
    (0.0, 0.24),
    "S",
    "Morante et al. 2020, Sec. 6.2, p. 457: 18-24 % solar-array power "
    "loss in all-electric LEO-GEO transfers",
    "fixed",
    "Applied in proportion to the electric share of the transfer.",
)

_add(
    "power_margin_50",
    30.0,
    "W",
    (15.0, 45.0),
    "D",
    "Power available to propulsion without an added power system, 50 kg bus",
    "fixed",
)
_add(
    "power_margin_150",
    100.0,
    "W",
    (50.0, 150.0),
    "D",
    "Power available to propulsion without an added power system, 150 kg bus",
    "oat",
)
_add(
    "power_margin_300",
    200.0,
    "W",
    (100.0, 300.0),
    "D",
    "Power available to propulsion without an added power system, 300 kg bus",
    "oat",
)

_add(
    "feed_frac",
    0.075,
    "-",
    (0.05, 0.10),
    "S",
    "Feed system as a fraction of propellant mass",
    "fixed",
)
_add(
    "struct_frac",
    0.20,
    "-",
    None,
    "S",
    "Structure margin on the propulsion dry mass",
    "fixed",
)

# Missions

_add("bus_mass_leo", 150.0, "kg", (75.0, 225.0), "D", "LEO mission bus mass", "oat")
_add("bus_mass_geo", 300.0, "kg", (150.0, 450.0), "D", "GEO mission bus mass", "oat")
_add(
    "mission_life_leo",
    5.0,
    "yr",
    (2.5, 7.5),
    "S",
    "LEO mission duration",
    "oat",
    "Scales the station-keeping budget and the schedule windows.",
)
_add("mission_life_geo", 15.0, "yr", None, "S", "GEO mission duration", "fixed")
_add("leo_altitude", 500e3, "m", None, "S", "LEO reference orbit", "fixed")

_add(
    "dv_sk_leo_rate",
    10.0,
    "m/s/yr",
    (5.0, 30.0),
    "S",
    "Drag compensation at 500 km, 5-30 m/s/yr with solar activity",
    "oat",
    "Station-keeping budget = rate x mission_life_leo.",
)
_add(
    "dv_ca_budget_1d",
    0.1,
    "m/s",
    (0.05, 0.15),
    "P",
    "Collision-avoidance delta-v per event at a 1 d warning time: 0.1 m/s, "
    "IRS-P6 case, Proc. 8th European Conference on Space Debris (2021)",
    "oat",
    "Per-event delta-v = dv_ca_budget_1d x (t_warn_ref_ca / t_warn), all modes.",
)
_add(
    "dv_ca_restore",
    0.1,
    "m/s",
    None,
    "P",
    "Return burn per collision-avoidance event, 0.1 m/s: IRS-P6 case, "
    "Proc. 8th European Conference on Space Debris (2021)",
    "fixed",
    "Not time-critical; flown electrically where the architecture can.",
)
_add(
    "t_warn_ref_ca",
    86400.0,
    "s",
    None,
    "E",
    "Reference warning time of dv_ca_budget_1d (1 d)",
    "fixed",
)
_add(
    "n_ca_events",
    15,
    "-",
    (5, 20),
    "S",
    "Collision-avoidance events over the LEO mission (3 per year)",
    "fixed",
)
_add(
    "dv_eol_leo",
    28.0,
    "m/s",
    (22.0, 60.0),
    "V",
    "Impulsive perigee lowering from 500 km circular to 500 x 400 km, "
    "vis-viva (27.92 m/s)",
    "oat",
)
_add(
    "eol_leo_ep_loss_ratio",
    56.0 / 28.0,
    "-",
    (1.25, 3.0),
    "V",
    "Low-thrust spiral from 500 to 400 km circular (55.95 m/s) over the "
    "impulsive perigee lowering (27.92 m/s)",
    "oat",
)

_add(
    "dv_transfer_chem",
    1800.0,
    "m/s",
    (1500.0, 1800.0),
    "P",
    "GTO to GEO, impulsive: 1800 m/s, Vulcan user's guide, standard Cape "
    "Canaveral GTO (27 deg)",
    "oat",
)
_add(
    "eta_planning_transfer",
    0.64,
    "-",
    (0.62, 0.65),
    "P",
    "GTO to GEO planning efficiency, dv_transfer_chem / dv_transfer_ep: "
    "0.64, Loghry et al., IEPC-2017-396. Bracket: Kluever 2012, 0.62 at "
    "30 d of electric transfer rising to about 0.65; Rovey et al. 2019, "
    "0.45-0.65 for Earth orbit",
    "oat",
    "Electric transfer budget = dv_transfer_chem / eta_planning_transfer.",
)

_add(
    "dv_geo_sk_rate",
    50.0,
    "m/s/yr",
    None,
    "S",
    "GEO station keeping, 50 m/s/yr",
    "fixed",
    "Station-keeping budget = rate x mission_life_geo.",
)
_add("dv_eol_geo", 11.0, "m/s", None, "S", "GEO graveyard raise (IADC)", "fixed")

# Feasibility limits and duty factors

_add(
    "sk_duty_max",
    0.5,
    "-",
    None,
    "E",
    "Station-keeping burn time at most 0.5 of the available window",
    "fixed",
)
_add(
    "ca_burn_twarn_frac",
    0.5,
    "-",
    None,
    "E",
    "Collision-avoidance burn at most 0.5 of the warning time",
    "fixed",
)
_add(
    "t_warn",
    86400.0,
    "s",
    (7200.0, 259200.0),
    "S",
    "Collision-avoidance warning time, 1 d",
    "axis",
    "Swept over 2 h, 6 h, 1 d and 3 d.",
)
_add(
    "eol_leo_window",
    5.0,
    "yr",
    None,
    "S",
    "LEO disposal window after the mission (FCC five-year rule; ESA Zero Debris)",
    "fixed",
)
_add(
    "eol_geo_window",
    0.25,
    "yr",
    None,
    "E",
    "GEO graveyard window, three months",
    "fixed",
)
_add(
    "t_transfer_max",
    0.5,
    "yr",
    (0.25, 2.0),
    "D",
    "GEO transfer window, six months (702SP-class electric orbit raising; "
    "Kluever 2012)",
    "axis",
    "Swept over 1 week to 24 months.",
)
_add(
    "w_duty_transfer",
    0.8,
    "-",
    (0.6, 1.0),
    "E",
    "Sunlit fraction during the GTO to GEO transfer",
    "oat",
)
_add(
    "w_duty_leo",
    0.63,
    "-",
    (0.55, 0.75),
    "V",
    "Sunlit orbit fraction at 500 km: 1 - asin(Re / (Re + h)) / pi = 0.622 at beta = 0",
    "oat",
    "Electric thrusting in sunlight only.",
)
_add(
    "w_duty_geo",
    1.0,
    "-",
    None,
    "V",
    "Sunlit fraction at GEO; eclipses only near the equinoxes, at most 72 min per day",
    "fixed",
)

_add(
    "eps_cap_frac",
    0.10,
    "-",
    (0.05, 0.20),
    "D",
    "Power-system limit: added power-system mass at most this fraction of the bus mass",
    "oat",
)

_add(
    "envelope_frac",
    0.5,
    "-",
    (0.35, 1.0),
    "D",
    "Bus mass-fraction limit: propulsion wet mass at most this fraction of "
    "the bus mass",
    "oat",
)

# Closed-form transfer check

_add(
    "geo_m_avg",
    350.0,
    "kg",
    None,
    "E",
    "Average spacecraft mass during the GEO transfer",
    "fixed",
)
_add(
    "geo_chem_stack",
    430.0,
    "kg",
    None,
    "E",
    "Wet stack of the chemical transfer",
    "fixed",
)


def dv_transfer_ep():
    """Electric GTO to GEO budget [m/s] = dv_transfer_chem /
    eta_planning_transfer."""
    return value("dv_transfer_chem") / value("eta_planning_transfer")


def dv_transfer_ep_range():
    """Electric transfer budget [m/s] at the two ends of the
    eta_planning_transfer bracket, (low, high)."""
    lo, hi = REGISTER["eta_planning_transfer"].range
    chem = value("dv_transfer_chem")
    return (chem / hi, chem / lo)


def power_margin_w(bus_class):
    """Power available to propulsion for a bus class [W]."""
    key = f"power_margin_{int(bus_class)}"
    if key not in REGISTER:
        raise KeyError(
            f"No power-margin row for bus class {bus_class} kg (defined: 50/150/300)."
        )
    return value(key)
