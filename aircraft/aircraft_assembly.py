import numpy as np
import config
import data.atr_reference as ref
from aircraft.fuselage import Fuselage
from aircraft.wing import Wing
from analysis_modules.cg_calculation import CenterOfGravity
from analysis_modules.htail_sizing import slopes
from aircraft.landing_gear import LandingGear
from aircraft.operations import Operations


class Aircraft:
    def __init__(self, aircraft_type: str, alpha: float, reynolds_duct: float, v_inf: float, mach: float,
                 reynolds_wing: float, reynolds_htail: float, power: str):
        self.aircraft_type = aircraft_type
        self.alpha = alpha
        self.reynolds_duct = reynolds_duct
        self.reynolds_wing = reynolds_wing
        self.reynolds_htail = reynolds_htail
        self.v_inf = v_inf
        self.mach = mach
        self.pc = power

        self.wing = Wing(b_wing=ref.b_w, sweep_wing=ref.phi_qc_w, wing_airfoil=ref.wing_airfoil,
                         taper_ratio_wing=ref.tr_w, c_root_wing=ref.c_root_w, alpha=self.alpha, velocity=self.v_inf,
                         mach=self.mach, reynolds_number=self.reynolds_wing, wing_type=self.aircraft_type)

        self.fuselage = Fuselage(fuselage_length=ref.l_cab + ref.l_cockpit + ref.l_cabin,
                                 fuselage_diameter=ref.diameter_fuselage, l_cabin=ref.l_cabin, l_cockpit=ref.l_cockpit,
                                 l_tail=ref.l_tail, velocity=self.v_inf, alpha=self.alpha, ref_area=self.wing.area(),
                                 reynolds_number=self.reynolds_wing, mach=self.mach, cl_wing=self.wing.cl_prime(),
                                 cmac_wing=ref.c_root_w)

        self.lg = LandingGear(aircraft_type=self.aircraft_type)

        self.operation = Operations(MTOM=ref.MTOW, pax=config.n_pax)

        self.empennage = self.make_empennage()

    def make_empennage(self):
        ar_wing = self.wing.aspect_ratio()
        cl_wing = self.wing.cl_prime()
        cla_wing = self.wing.cl_al()

        if self.aircraft_type == "DUUC":
            from aircraft.propulsive_empennage.empennage_assembly_PE import PropulsiveEmpennage
            va_inlet = self.v_inf * (np.pi * (config.duct_diameter / 2))
            print(f"va inlet: {va_inlet}")
            return PropulsiveEmpennage(rpm=config.rpm,
                                       alpha=self.alpha,
                                       power_condition=self.pc,
                                       va_inlet=va_inlet,
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
                                       ref_area=ref.s_w,
                                       ref_chord=ref.c_mac_w,
                                       ar_wing=ar_wing,
                                       cl_wing=cl_wing,
                                       cla_wing=cla_wing)
        if self.aircraft_type == "conventional":
            from aircraft.conventional_empennage.empennage_assembly_conv import ConventionalEmpennage
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
                                         power_condition=self.pc,
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
                                         reynolds=self.reynolds_htail,
                                         ar_wing=ar_wing,
                                         cl_wing=cl_wing,
                                         cla_wing=cla_wing)
        else:
            print("Wrong aircraft type defined!!")
            return None

    def x_cog(self):
        if self.aircraft_type == "conventional":
            cog = CenterOfGravity(w_lg_nose=self.lg.weight()[0],
                                  w_lg_main=self.lg.weight()[1],
                                  w_eng_wing=self.empennage.weight_prop() / 2,
                                  w_wing=self.wing.weight(),
                                  w_fuse=self.fuselage.weight(),
                                  w_nac_wing=self.empennage.weight_nac() / 2,
                                  w_vt=self.empennage.vt_tail.weight(),
                                  w_ht=self.empennage.ht_tail.weight(),
                                  w_duct=0,
                                  w_sys=self.operation.weight_sys(),
                                  l_fuselage=(ref.l_tail + ref.l_cockpit + ref.l_cabin),
                                  aircraft_type=self.aircraft_type,
                                  c_mac_wing=ref.c_mac_w)
            return cog.x_cg()[0], cog.x_cg()[1]
        if self.aircraft_type == "DUUC":
            cog = CenterOfGravity(w_lg_nose=self.lg.weight()[0],
                                  w_lg_main=self.lg.weight()[1],
                                  w_eng_wing=0,
                                  w_wing=self.wing.weight(),
                                  w_fuse=self.fuselage.weight(),
                                  w_nac_wing=0, w_vt=0, w_ht=0,
                                  w_duct=self.empennage.weight() * 2,
                                  w_sys=self.operation.weight_sys(),
                                  l_fuselage=(ref.l_tail + ref.l_cockpit + ref.l_cabin),
                                  aircraft_type=self.aircraft_type,
                                  c_mac_wing=ref.c_mac_w)
            return cog.x_cg()[0], cog.x_cg()[1]
        else:
            return 0, 0

    def cd0_vector(self):
        if self.aircraft_type == "conventional":
            cd0_nacelle = self.empennage.cd0_vector()[0] * 2
            cd0_ht = self.empennage.cd0_vector()[1]
            cd0_vt = self.empennage.cd0_vector()[2]
            cd0_wing = self.wing.cd0()
            cd0_fuselage = self.fuselage.cd0()
            return [cd0_nacelle, cd0_ht, cd0_vt, cd0_wing, cd0_fuselage]
        if self.aircraft_type == "DUUC":
            cd0_duct = self.empennage.cd0_vector()[0] * 2
            cd0_pylon = self.empennage.cd0_vector()[1] * 2
            cd0_nacelle = self.empennage.cd0_vector()[2] * 2
            cd0_support = self.empennage.cd0_vector()[3] * 2
            cd0_control = self.empennage.cd0_vector()[4] * 8
            cd0_wing = self.wing.cd0()
            cd0_fuselage = self.fuselage.cd0()
            return [cd0_duct, cd0_pylon, cd0_nacelle, cd0_support, cd0_control, cd0_wing, cd0_fuselage]

    def cd0_empennage(self):
        if self.aircraft_type == "conventional":
            cd0_ht = self.empennage.cd0_vector()[1]
            cd0_vt = self.empennage.cd0_vector()[2]
            return [cd0_ht, cd0_vt]
        if self.aircraft_type == "DUUC":
            cd0_duct = self.empennage.cd0_vector()[0] * 2
            cd0_pylon = self.empennage.cd0_vector()[1] * 2
            cd0_support = self.empennage.cd0_vector()[3] * 2
            cd0_control = self.empennage.cd0_vector()[4] * 8
            return [cd0_duct, cd0_pylon, cd0_support, cd0_control,]

    def cl_vector(self):
        """ not normalized values """
        if self.aircraft_type == "conventional":
            cl_fuselage = self.fuselage.cl()
            cl_wing = self.wing.cl()
            cl_ht = self.empennage.ht_tail.cl()
            return [cl_ht, cl_wing, cl_fuselage]

        if self.aircraft_type == "DUUC":
            cl_fuselage = self.fuselage.cl()
            cl_wing = self.wing.cl()
            cl_duct = self.empennage.duct.cl() * 2
            cl_pylon = self.empennage.pylon.cl() * 2
            cl_support = self.empennage.support.cl() * 2
            cl_control = self.empennage.rudder.cl() * 4
            return [cl_duct, cl_pylon, cl_support, cl_control, cl_wing, cl_fuselage]

    def cl_ac(self):
        if self.aircraft_type == "conventional":
            cl_fuselage = self.fuselage.cl_prime()
            cl_wing = self.wing.cl_prime()
            cl_ht = self.empennage.ht_tail.cl_prime()

            cl_tot = cl_fuselage + cl_wing + cl_ht
            return cl_tot

        if self.aircraft_type == "DUUC":
            cl_fuselage = self.fuselage.cl_prime()
            cl_wing = self.wing.cl_prime()
            cl_empennage = self.empennage.cl_prime()

            cl_tot = cl_fuselage + cl_wing + cl_empennage
            return cl_tot

    def slope_htail(self):
        cog = self.x_cog()[0]
        x_lemac = self.x_cog()[1]

        print(f"cog: {cog}, xlemac: {x_lemac}")
        print(f"diff: {cog - x_lemac}")

        cl_h = self.empennage.cl_prime()
        cl = self.cl_ac()

        cl_h = 4.586
        cl = 1.44

        a1 = slopes("control", self.aircraft_type, ref.z_h, ref.phi_qc_w, self.wing.aspect_ratio(),
                    self.fuselage.length(), x_lemac, 2.234, 0.9, cl_h, cl, ref.tr_w, ref.b_w, 0,
                    self.mach, ref.ar_h, 0, 0)[0]

        b1 = slopes("control", self.aircraft_type, ref.z_h, ref.phi_qc_w, self.wing.aspect_ratio(),
                    self.fuselage.length(), x_lemac, 2.234, 0.9, cl_h, cl, ref.tr_w, ref.b_w, 0,
                    self.mach, ref.ar_h, 0, 0)[1]

        a2 = slopes("stability", self.aircraft_type, ref.z_h, ref.phi_qc_w, self.wing.aspect_ratio(),
                    self.fuselage.length(), x_lemac, 2.234, 0.9, cl_h, cl, ref.tr_w, ref.b_w, 0,
                    self.mach, ref.ar_h, 0, 0)[0]

        return a1, b1, a2

    def thrust(self):
        t_ac = 0
        return t_ac


""" Test section"""
"""
if __name__ == "__main__":
    DUUC = Aircraft(aircraft_type="DUUC", alpha=0,
                    reynolds=8422274, v_inf=128,
                    mach=0.576)

    print(f"DUUC drag: {DUUC.empennage.cd_prime()}")"""
