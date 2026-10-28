"""
CEA cross-check of the 271 s datasheet specific impulse of the N2O/propene
bipropellant thruster: ideal vacuum Isp at O/F = 7, Pc = 10-20 bar and area
ratio 50-100, with a delivered efficiency of 0.92-0.96. Requires RocketCEA.
"""

from rocketcea.cea_obj_w_units import CEA_Obj

cea = CEA_Obj(oxName="N2O", fuelName="Propylene", pressure_units="bar", isp_units="sec")
OF = 7.0
for pc in (10.0, 15.0, 20.0):
    for eps in (50.0, 100.0):
        isp_vac = cea.get_Isp(Pc=pc, MR=OF, eps=eps)
        print(
            f"Pc {pc:4.0f} bar eps {eps:5.0f}: Isp_vac,ideal "
            f"{isp_vac:6.1f} s -> delivered (x0.92-0.96) "
            f"{0.92 * isp_vac:5.1f}-{0.96 * isp_vac:5.1f} s"
        )
