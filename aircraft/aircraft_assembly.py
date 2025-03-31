import numpy as np
import config
import data.atr_reference as ref
from aircraft.fuselage import Fuselage
from aircraft.wing import Wing
from analysis_modules.cg_calculation import CenterOfGravity
from analysis_modules.htail_sizing import slopes
from aircraft.landing_gear import LandingGear
from aircraft.operations import Operations
from analysis_modules.factors import area_ratio


class Aircraft:
    def __init__(self, aircraft_type: str, alpha: float, v_inf: float, mach: float, power: str, bem_input,
                 altitude: float, delta_e: float, delta_r: float, wing_span: float, wing_sweep: float, wing_airfoil: str,
                 wing_tr: float, wing_cr: float, l_coc: float, l_cab: float, l_tail: float, fus_d: float, pax: float,
                 d_duct: float, c_duct: float, duct_airfoil: str, b_pylon: float, c_pylon: float, pylon_airfoil: str,
                 cant_angle: float, b_support: float, c_support: float, support_airfoil: str, cv_airfoil: str,
                 n_blades: float, l_nacelle: float, d_nacelle: float, l_cv: float, c_cv: float, propulsor_type: str,
                 prop_airfoil: str, d_hub: float, rpm: float, prop_c_root: float, prop_c_tip: float, prop_diameter: float):
        self.aircraft_type = aircraft_type
        self.alpha = alpha
        self.v_inf = v_inf
        self.mach = mach
        self.pc = power
        self.bem_input = bem_input
        self.altitude = altitude
        self.delta_e = delta_e
        self.delta_r = delta_r
        self.wing_span = wing_span
        self.wing_sweep = wing_sweep
        self.wing_airfoil = wing_airfoil
        self.wing_tr = wing_tr
        self.wing_cr = wing_cr
        self.l_coc = l_coc
        self.l_cab = l_cab
        self.l_tail = l_tail
        self.l_fuse = self.l_coc + self.l_coc + self.l_cab
        self.fus_d = fus_d
        self.pax = pax
        self.d_duct = d_duct
        self.c_duct = c_duct
        self.duct_airfoil = duct_airfoil
        self.l_pylon = b_pylon
        self.c_pylon = c_pylon
        self.pylon_airfoil = pylon_airfoil
        self.cant_angle = cant_angle
        self.b_support = b_support
        self.c_support = c_support
        self.support_airfoil = support_airfoil
        self.control_vanes_airfoil = cv_airfoil
        self.n_blades = n_blades
        self.l_nacelle = l_nacelle
        self.d_nacelle = d_nacelle
        self.control_vane_length = l_cv
        self.control_vane_chord = c_cv
        self.propulsor_type = propulsor_type
        self.hub_diameter = d_hub
        self.propeller_airfoil = prop_airfoil
        self.RPM = rpm
        self.c_root = prop_c_root
        self.c_tip = prop_c_tip
        self.prop_diameter = prop_diameter

        self.wing = Wing(b_wing=self.wing_span, sweep_wing=self.wing_sweep, wing_airfoil=self.wing_airfoil,
                         taper_ratio_wing=self.wing_tr, c_root_wing=self.wing_cr, alpha=self.alpha, velocity=self.v_inf,
                         mach=self.mach, wing_type=self.aircraft_type, altitude=self.altitude)

        self.fuselage = Fuselage(fuselage_length=self.l_fuse,
                                 fuselage_diameter=self.fus_d, l_cabin=self.l_cab, l_cockpit=self.l_coc,
                                 l_tail=self.l_tail, velocity=self.v_inf, alpha=self.alpha, ref_area=self.wing.area(),
                                 mach=self.mach, cl_wing=self.wing.cl_prime(),
                                 cmac_wing=self.wing_cr, altitude=self.altitude)

        self.lg = LandingGear(aircraft_type=self.aircraft_type)

        self.operation = Operations(MTOM=ref.MTOW, pax=self.pax)

        self.empennage = self.make_empennage()

    def make_empennage(self):
        ar_wing = self.wing.aspect_ratio()
        cl_wing = self.wing.cl_prime()
        cla_wing = self.wing.cl_al()

        if self.aircraft_type == "DUUC":
            from aircraft.propulsive_empennage.empennage_assembly_PE import PropulsiveEmpennage
            va_inlet = self.v_inf * (np.pi * (self.d_duct / 2))
            return PropulsiveEmpennage(rpm=self.RPM,
                                       alpha=self.alpha,
                                       power_condition=self.pc,
                                       va_inlet=va_inlet,
                                       n_blades=self.n_blades,
                                       prop_diameter=self.prop_diameter,
                                       hub_diameter=self.hub_diameter,
                                       prop_airfoil=self.propeller_airfoil,
                                       prop_sweep=0,
                                       prop_pitch=0,
                                       c_root=self.c_root,
                                       c_tip=self.c_tip,
                                       duct_diameter=self.d_duct,
                                       duct_chord=self.c_duct,
                                       duct_profile=self.duct_airfoil,
                                       cant_angle=self.cant_angle,
                                       pylon_length=self.l_pylon,
                                       pylon_chord=self.c_pylon,
                                       pylon_profile=self.pylon_airfoil,
                                       d_exit=area_ratio("0016", self.c_duct, self.d_duct/2, 1)[0],
                                       nacelle_length=self.l_nacelle,
                                       nacelle_diameter=self.d_nacelle,
                                       support_length=self.b_support,
                                       support_chord=self.c_support,
                                       support_profile=self.support_airfoil,
                                       hcv_span=self.control_vane_length,
                                       hcv_chord=self.control_vane_chord,
                                       control_profile=self.control_vanes_airfoil,
                                       vcv_span=self.control_vane_length,
                                       vcv_chord=self.control_vane_chord,
                                       propulsor_type=self.propulsor_type,
                                       v_inf=self.v_inf,
                                       mach=self.mach,
                                       ref_area=ref.s_w,
                                       ref_chord=ref.c_mac_w,
                                       ar_wing=ar_wing,
                                       cl_wing=cl_wing,
                                       cla_wing=cla_wing,
                                       bem_input=self.bem_input,
                                       altitude=self.altitude,
                                       delta_e=self.delta_e,
                                       delta_r=self.delta_r)
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
                                         ar_wing=ar_wing,
                                         cl_wing=cl_wing,
                                         cla_wing=cla_wing,
                                         bem_input=self.bem_input,
                                         altitude=self.altitude)
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
                                  l_fuselage=(self.l_tail + self.l_coc + self.l_cab),
                                  aircraft_type=self.aircraft_type,
                                  c_mac_wing=ref.c_mac_w)
            return cog.x_cg()[0], cog.x_cg()[1], cog.cg_fuselage_group()[0], cog.cg_wing_group()[0]
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
            return cog.x_cg()[0], cog.x_cg()[1], cog.cg_fuselage_group()[0], cog.cg_wing_group()[0]
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
