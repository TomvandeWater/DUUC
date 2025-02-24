import numpy as np
import config
import data.atr_reference as ref
from fuselage import Fuselage
from wing import Wing


class Aircraft:
    def __init__(self, aircraft_type: str, alpha: float, reynolds: float, v_inf: float, mach: float):
        self.aircraft_type = aircraft_type
        self.alpha = alpha
        self.reynolds = reynolds
        self.v_inf = v_inf
        self.mach = mach

        self.wing = Wing(b_wing=ref.b_w, sweep_wing=ref.phi_qc_w, wing_airfoil=ref.wing_airfoil,
                         taper_ratio_wing=ref.tr_w, c_root_wing=ref.c_root_w, alpha=self.alpha, velocity=self.v_inf,
                         mach=self.mach, reynolds_number=self.reynolds, wing_type=self.aircraft_type)

        self.fuselage = Fuselage(fuselage_length=ref.l_cab + ref.l_cockpit + ref.l_cabin,
                                 fuselage_diameter=ref.diameter_fuselage, l_cabin=ref.l_cabin, l_cockpit=ref.l_cockpit,
                                 l_tail=ref.l_tail, velocity=self.v_inf, alpha=self.alpha, ref_area=self.wing.area(),
                                 reynolds_number=self.reynolds, mach=self.mach, cl_wing=self.wing.cl_prime(),
                                 cmac_wing=ref.c_root_w)

        self.empennage = self.make_empennage()

    def make_empennage(self):
        if self.aircraft_type == "DUUC":
            from propulsive_empennage.empennage_assembly_PE import PropulsiveEmpennage
            va_inlet = self.v_inf * (np.pi * (config.duct_diameter / 2))
            return PropulsiveEmpennage(rpm=config.rpm,
                                       alpha=self.alpha,
                                       power_condition=config.power_condition,
                                       va_inlet=va_inlet,
                                       re_duct=self.reynolds,
                                       n_blades=config.n_blades,
                                       prop_diameter=config.duct_diameter * 0.95,
                                       hub_diameter=config.hub_diameter,
                                       prop_airfoil=config.prop_airfoil,
                                       prop_sweep=config.propeller_sweep,
                                       prop_pitch=config.propeller_pitch,
                                       c_root=config.c_root,
                                       c_tip=config.c_tip,
                                       duct_diameter=config.duct_diameter,
                                       duct_chord=config.duct_chord,
                                       duct_profile=config.duct_airfoil,
                                       cant_angle=config.cant_angle,
                                       pylon_length=config.pylon_length,
                                       pylon_chord=config.pylon_chord,
                                       pylon_profile=config.pylon_airfoil,
                                       d_exit=config.d_exit,
                                       nacelle_length=config.nacelle_length,
                                       nacelle_diameter=config.nacelle_diameter,
                                       support_length=config.support_length,
                                       support_chord=config.support_chord,
                                       support_profile=config.support_airfoil,
                                       hcv_span=config.control_vane_length,
                                       hcv_chord=config.control_vane_chord,
                                       control_profile=config.control_vanes_airfoil,
                                       vcv_span=config.control_vane_length,
                                       vcv_chord=config.control_vane_chord,
                                       propulsor_type=config.propulsor_type,
                                       v_inf=self.v_inf,
                                       mach=self.mach,
                                       ref_area=ref.s_w)
        if self.aircraft_type == "conventional":
            from conventional_empennage.empennage_assembly_conv import ConventionalEmpennage
            return ConventionalEmpennage(ht_span=ref.b_h,
                                         ht_chord=ref.c_root_h,
                                         ht_profile=ref.airfoil_ht,
                                         ht_taper=ref.tr_h,
                                         ht_sweep=ref.phi_qc_h,
                                         ht_croot=ref.c_root_h,
                                         vt_span=ref.b_v,
                                         vt_chord=ref.c_root_v,
                                         vt_profile=ref.airfoil_vt,
                                         vt_taper=ref.tr_v,
                                         vt_sweep=ref.phi_qc_v,
                                         vt_croot=ref.c_root_v,
                                         tail_type="t-tail",
                                         rpm=config.rpm,
                                         alpha=self.alpha,
                                         power_condition=config.power_condition,
                                         n_blades=config.n_blades,
                                         prop_diameter=config.duct_diameter * 0.95,
                                         hub_diameter=config.hub_diameter,
                                         prop_airfoil=config.prop_airfoil,
                                         prop_sweep=config.propeller_sweep,
                                         prop_pitch=config.propeller_pitch,
                                         c_root=config.c_root,
                                         c_tip=config.c_tip,
                                         nacelle_length=ref.l_nacelle,
                                         nacelle_diameter=ref.d_nacelle,
                                         v_inf=self.v_inf,
                                         mach=self.mach,
                                         reynolds=self.reynolds)
        else:
            print("Wrong aircraft type defined!!")
            return None

    def cd(self):
        cd_ac = 0
        return cd_ac

    def cl(self):
        cl_ac = 0
        return cl_ac

    def thrust(self):
        t_ac = 0
        return t_ac


""" Test section"""

if __name__ == "__main__":
    DUUC = Aircraft(aircraft_type="DUUC", alpha=0,
                    reynolds=8422274, v_inf=128,
                    mach=0.576)

    print(f"DUUC drag: {DUUC.empennage.cd_prime()}")
