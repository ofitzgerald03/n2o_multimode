"""
Matplotlib defaults (applied on import) and the colours, linestyles and
names of the propellants, architectures and loss categories.
"""

from typing import Mapping

import matplotlib as mpl


THESIS_RCPARAMS = {
    "figure.figsize": (8.0, 5.5),
    "figure.dpi": 100,
    "savefig.dpi": 150,
    "savefig.bbox": "tight",
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 12,
    "legend.fontsize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "lines.linewidth": 2.0,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linestyle": "-",
    "legend.framealpha": 0.9,
}


def apply_thesis_style():
    """Apply THESIS_RCPARAMS to matplotlib."""
    mpl.rcParams.update(THESIS_RCPARAMS)


ARCHITECTURE_COLOURS = {
    "single_elec": "tab:blue",
    "multimode": "tab:green",
    "hybrid": "tab:orange",
    "single_chem": "tab:red",
    "mono_multimode": "tab:purple",
}

ARCHITECTURE_LINESTYLES = {
    "single_elec": (0, (1.2, 1.3)),
    "multimode": "-",
    "hybrid": (0, (6.0, 1.6, 1.4, 1.6)),
    "single_chem": "--",
    "mono_multimode": (0, (3.0, 1.2, 1.0, 1.2, 1.0, 1.2)),
}

ARCHITECTURE_NAMES = {
    "single_elec": "Single-electric",
    "single_chem": "Single-chemical",
    "hybrid": "Hybrid",
    "multimode": "Biprop. multimode",
    "mono_multimode": "Monoprop. multimode",
}


apply_thesis_style()


PROPELLANT_STYLE: Mapping[str, dict] = {
    "Xe": {
        "color": "#009E73",  # green
        "linestyle": "-",
        "label": "Xe",
    },
    "Ar": {
        "color": "#7B3294",  # purple
        "linestyle": (0, (6.0, 1.6, 1.4, 1.6)),
        "label": "Ar",
    },
    "N2O": {
        "color": "#0072B2",  # blue
        "linestyle": "--",
        "label": r"N$_2$O",
    },
    "H2O": {
        "color": "#E69F00",  # orange
        "linestyle": (0, (1.2, 1.3)),
        "label": r"H$_2$O",
    },
}

# Second shade of each propellant colour, for literature data drawn
# against the model curve
PROPELLANT_SHADE = {
    "Xe": "#005C43",
    "Ar": "#C29BD1",
    "N2O": "#003F63",
    "H2O": "#8C5F00",
}

# Colours of the energy-loss categories
CATEGORY_COLOURS = {
    "rotation": "#2ca02c",
    "vibration": "#9467bd",
    "elastic": "#0e7f8a",
    "attachment": "#7f7f7f",
    "electronic excitation": "#1f77b4",
    "ionisation": "#d62728",
}
