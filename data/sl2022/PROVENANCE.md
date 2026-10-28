# Sheppard and Little (2022) digitised data

Source: Sheppard and Little (2022), "Performance Analysis of an Electron
Cyclotron Resonance Thruster with Various Propellants", Figs. 7a (argon),
7c (xenon) and 7d (water): electron temperature and exit velocity against
mdot/P, per mass-flow group.

Files: `SL2022_<propellant>_Te.csv` and `SL2022_<propellant>_ue.csv` for
Ar, H2O and Xe. Columns: `mdot_mg_s` (mass-flow group),
`mdot_over_P_ug_per_J`, and `Te_eV` [eV] or `ue_m_s` [m/s].

T_e and u_e are digitised as separate curves; `multimode/thruster/
sl2022_validation.py` pairs them within each mass-flow group.
