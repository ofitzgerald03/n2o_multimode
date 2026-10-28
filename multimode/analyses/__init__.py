"""
Analysis tasks run through run.py; importing the package registers them in
the order below.
"""

from . import thruster_results  # noqa: F401 # thruster_anchor, ratio_table, closure_spread
from . import thruster_figures  # noqa: F401 # thruster_figures, sl2022_validation
from . import thruster_checks  # noqa: F401 # check_thruster
from . import spacecraft_results  # noqa: F401 # spacecraft_results, spacecraft_sweeps
from . import spacecraft_figures  # noqa: F401 # response_variant, spacecraft_figures
from . import spacecraft_checks  # noqa: F401 # check_spacecraft

from .registry import TASKS, get, all_task_names  # noqa: F401
