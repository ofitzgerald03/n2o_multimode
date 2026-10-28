"""
Composition of absolute thrust, mass flow and power from a measured xenon
operating point (Cannat et al. 2015, configuration 2, Table II, 0.1 mg/s).
"""

# For Q in {T, m_dot_ion, P_in}:
#
#     Q_target(Te) = Q_xe,meas * Q_fw,target(Te) / Q_fw,xe(Te_ref)
#
# where Q_fw is the framework prediction per unit n0 A0, which cancels. The
# total mass flow is m_dot_ion / eta_u.

from dataclasses import dataclass
import numpy as np

from ..constants import E_CHARGE, G0


@dataclass(frozen=True)
class AnchorReference:
    """Measured xenon operating point.

    T_meas [N], mdot_meas [kg/s] (total), Pin_meas [W] (absorbed),
    Isp_meas [s], eta_u_meas [-], Te_ref [eV], Rm_ref [-].
    """

    T_meas: float
    mdot_meas: float
    Pin_meas: float
    Isp_meas: float
    eta_u_meas: float
    Te_ref: float
    Rm_ref: float

    def mdot_ion_meas(self) -> float:
        """Ion mass flow [kg/s] = mdot_meas * eta_u_meas."""
        return self.mdot_meas * self.eta_u_meas


def framework_extensive_arrays(q, eta_T, eta_RF, eta_D):
    """Thrust, ion mass flow and absorbed power per unit n0 A0.

    q : per-particle bundle ('v_b_flux', 'm_i_eff', 'u_e', 'phi_sh').
    eta_T : thrust efficiency against Te at the reference mirror ratio.

    T = v_b m_i u_e / 2;  mdot_ion = v_b m_i / 2;
    P_in = (v_b e phi_sh / 2) / (eta_T eta_RF eta_D).
    """
    v_b = q["v_b_flux"]
    m_i = q["m_i_eff"]
    u_e = q["u_e"]
    phi_eV = q["phi_sh"]

    half = 0.5
    T_fw = half * v_b * m_i * u_e
    mdot_ion_fw = half * v_b * m_i
    P_beam_fw = half * v_b * E_CHARGE * phi_eV

    eta_tot = eta_T * eta_RF * eta_D
    with np.errstate(divide="ignore", invalid="ignore"):
        Pin_fw = np.where(eta_tot > 0, P_beam_fw / eta_tot, np.inf)

    return {
        "T_fw": T_fw,
        "mdot_ion_fw": mdot_ion_fw,
        "Pin_fw": Pin_fw,
    }


def compute_scale_factors(fw_arrays_xe, Te_array, anchor_ref):
    """Scale factors Q_meas,xe / Q_fw,xe(Te_ref) for T, mdot_ion and P_in,
    and the index of the grid point closest to Te_ref."""
    i_anchor = int(np.argmin(np.abs(Te_array - anchor_ref.Te_ref)))

    T_fw_xe_anchor = float(fw_arrays_xe["T_fw"][i_anchor])
    mdot_ion_fw_xe_anchor = float(fw_arrays_xe["mdot_ion_fw"][i_anchor])
    Pin_fw_xe_anchor = float(fw_arrays_xe["Pin_fw"][i_anchor])

    if T_fw_xe_anchor <= 0 or mdot_ion_fw_xe_anchor <= 0 or Pin_fw_xe_anchor <= 0:
        raise ValueError(
            "Non-positive Xe framework quantity at the anchor Te; cannot "
            "compute scale factors. Check Te_ref against the Te grid."
        )

    return {
        "scale_T": anchor_ref.T_meas / T_fw_xe_anchor,
        "scale_mdot_ion": anchor_ref.mdot_ion_meas() / mdot_ion_fw_xe_anchor,
        "scale_Pin": anchor_ref.Pin_meas / Pin_fw_xe_anchor,
        "i_anchor": i_anchor,
    }


def compose_extensive_curves(fw_arrays_target, scale_factors, eta_u_target):
    """Absolute-scale curves of a target propellant against Te.

    Returns 'T' [N], 'm_dot' [kg/s] (total), 'm_dot_ion' [kg/s], 'P_in' [W],
    'Isp_gross' = T / (m_dot g0) [s] and 'Isp_intensive' = u_e / g0 [s].
    """
    T_arr = scale_factors["scale_T"] * fw_arrays_target["T_fw"]
    mdot_ion_arr = scale_factors["scale_mdot_ion"] * fw_arrays_target["mdot_ion_fw"]
    mdot_total_arr = mdot_ion_arr / eta_u_target
    Pin_arr = scale_factors["scale_Pin"] * fw_arrays_target["Pin_fw"]

    with np.errstate(divide="ignore", invalid="ignore"):
        Isp_gross = np.where(mdot_total_arr > 0, T_arr / (mdot_total_arr * G0), np.inf)

    # u_e = T_fw / mdot_ion_fw
    with np.errstate(divide="ignore", invalid="ignore"):
        u_e = np.where(
            fw_arrays_target["mdot_ion_fw"] > 0,
            fw_arrays_target["T_fw"] / fw_arrays_target["mdot_ion_fw"],
            0.0,
        )
    Isp_intensive = u_e / G0

    return {
        "T": T_arr,
        "m_dot": mdot_total_arr,
        "m_dot_ion": mdot_ion_arr,
        "P_in": Pin_arr,
        "Isp_gross": Isp_gross,
        "Isp_intensive": Isp_intensive,
    }


def at_anchor_operating_point(
    composed_curves, fw_arrays_xe, fw_arrays_target, anchor_ref, eta_u_target, i_anchor
):
    """Operating point at the anchor and the propellant ratios to xenon.

    R_T, R_mdot_ion and R_Pin are framework ratios target / Xe at Te_ref;
    R_ue = R_T / R_mdot_ion. Returns the scalars 'T', 'm_dot', 'm_dot_ion',
    'P_in', 'Isp_gross', 'Isp_intensive' and 'ratios'.
    """
    T_target = float(composed_curves["T"][i_anchor])
    mdot_total = float(composed_curves["m_dot"][i_anchor])
    mdot_ion = float(composed_curves["m_dot_ion"][i_anchor])
    Pin_target = float(composed_curves["P_in"][i_anchor])
    Isp_gross = float(composed_curves["Isp_gross"][i_anchor])
    Isp_intensive = float(composed_curves["Isp_intensive"][i_anchor])

    T_fw_xe = float(fw_arrays_xe["T_fw"][i_anchor])
    mdot_ion_fw_xe = float(fw_arrays_xe["mdot_ion_fw"][i_anchor])
    Pin_fw_xe = float(fw_arrays_xe["Pin_fw"][i_anchor])

    T_fw_target = float(fw_arrays_target["T_fw"][i_anchor])
    mdot_ion_fw_target = float(fw_arrays_target["mdot_ion_fw"][i_anchor])
    Pin_fw_target = float(fw_arrays_target["Pin_fw"][i_anchor])

    R_T = T_fw_target / T_fw_xe
    R_mdot_ion = mdot_ion_fw_target / mdot_ion_fw_xe
    R_Pin = Pin_fw_target / Pin_fw_xe
    R_ue = R_T / R_mdot_ion

    return {
        "T": T_target,
        "m_dot": mdot_total,
        "m_dot_ion": mdot_ion,
        "P_in": Pin_target,
        "Isp_gross": Isp_gross,
        "Isp_intensive": Isp_intensive,
        "ratios": {
            "R_T": R_T,
            "R_mdot_ion": R_mdot_ion,
            "R_ue": R_ue,
            "R_Pin": R_Pin,
        },
    }
