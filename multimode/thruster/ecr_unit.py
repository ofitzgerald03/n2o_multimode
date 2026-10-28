"""
Composed N2O electrodeless (ECR) unit at the anchor as a function of the
mass utilisation efficiency eta_u, for the spacecraft-level model.
"""

# Thrust and absorbed power are held at their anchor values; the total mass
# flow scales as 1/eta_u and the specific impulse as eta_u.

from dataclasses import dataclass

from ..config import PARAMS
from .anchoring import compose_extensive_curves, at_anchor_operating_point
from .pipeline import run_pipeline


@dataclass(frozen=True)
class EcrOperatingPoint:
    """Operating point of the N2O unit at one eta_u (absorbed-power side)."""

    eta_u: float
    thrust: float  # [N]
    isp: float  # [s], T / (mdot g0)
    mdot: float  # [kg/s], total propellant flow
    p_abs: float  # [W], absorbed RF power


def eta_u_parity() -> float:
    """Measured xenon utilisation of the anchor (0.61)."""
    return PARAMS["eta_u_xe_ref"]


_STATE = {}


def _pipeline_state():
    """Framework arrays, scale factors and anchor of one pipeline run,
    cached; none of them depend on eta_u."""
    if "fw_n2o" not in _STATE:
        params = dict(PARAMS)
        params["eta_u_n2o"] = eta_u_parity()
        res = run_pipeline(params, propellants=("Xe", "N2O"))
        case = params.get("n2o_case", "B")
        _STATE.update(
            {
                "fw_xe": res["Xe"]["fw_A"],
                "fw_n2o": res["N2O"][f"fw_{case}"],
                "scale_factors": res["anchor"]["scale_factors"],
                "i_anchor": res["anchor"]["i_anchor"],
                "anchor_ref": res["anchor"]["reference"],
            }
        )
    return _STATE


def ecr_unit_at(eta_u: float) -> EcrOperatingPoint:
    """N2O unit operating point at the anchor for 0 < eta_u <= 1."""
    if not 0.0 < eta_u <= 1.0:
        raise ValueError(f"eta_u = {eta_u} outside (0, 1]")
    s = _pipeline_state()
    composed = compose_extensive_curves(s["fw_n2o"], s["scale_factors"], eta_u)
    op = at_anchor_operating_point(
        composed, s["fw_xe"], s["fw_n2o"], s["anchor_ref"], eta_u, s["i_anchor"]
    )
    return EcrOperatingPoint(
        eta_u=eta_u,
        thrust=op["T"],
        isp=op["Isp_gross"],
        mdot=op["m_dot"],
        p_abs=op["P_in"],
    )
