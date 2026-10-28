"""
Thruster-model inputs (PARAMS), file locations (PATHS) and run options.
"""

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent

PATHS = {
    "project_root": PROJECT_ROOT,
    "data": PROJECT_ROOT / "data",
    "output": PROJECT_ROOT / "output",
    "figures": PROJECT_ROOT / "output" / "figures",
    "lxcat": PROJECT_ROOT / "data" / "lxcat",
    "xe_lxcat": PROJECT_ROOT / "data" / "lxcat" / "xe_hayashi.txt",
    "ar_lxcat": PROJECT_ROOT / "data" / "lxcat" / "ar_hayashi.txt",
    "cross_sections": PROJECT_ROOT / "data" / "cross_sections",
    "ps2017": PROJECT_ROOT / "data" / "ps2017",
    "sl2022": PROJECT_ROOT / "data" / "sl2022",
}


PARAMS = {
    # Hardware efficiencies
    "eta_RF": 0.90,  # Lafleur 2014; Cannat et al. 2015
    "eta_D": 0.73,  # Sheppard and Little 2022, Sec. IV
    # Mass utilisation efficiency
    # Cannat et al. 2015 Table II value for both propellants (parity).
    "eta_u_xe": 0.61,
    "eta_u_n2o": 0.61,
    # Anchor: Cannat et al. 2015, configuration 2, 0.1 mg/s row
    "T_ref": 0.98e-3,  # N        (Table II)
    "mdot_ref": 0.1e-6,  # kg/s     (Table II)
    "Pin_ref": 30.0,  # W        (Table II)
    "Isp_ref": 993.0,  # s        (Table II)
    "eta_u_xe_ref": 0.61,  #          (Table II)
    "Te_ref": 22.0,  # eV       (Fig. 8, 0.1 mg/s)
    "Rm_ref": 1.70,  #          (Fig. 2(b): B(20 mm) / B(35 mm) = 850 / 500)
    # Mirror ratio for the sensitivity curves
    "Rm_nominal": 5.0,
    # N2O Bohm closure: 'B' common system Bohm velocity, 'A' species-resolved
    "n2o_case": "B",
    # Electron-temperature grid [eV]
    "Te_min": 1.0,
    "Te_max": 30.0,
    "Te_step": 0.1,
    "spot_Te": (2, 3, 5, 7, 10, 12, 15, 20, 22, 25, 30),
}


RUN_OPTIONS = {
    "save_csvs": True,  # run.py --no-save overrides this
}
