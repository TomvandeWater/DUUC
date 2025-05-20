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
from analysis_modules.ISA import air_density_isa
from analysis_modules.aerodynamic import tail_volume
import matplotlib.pyplot as plt


class Aircraft:
    def __init__(self, conditions, reference, geometry_duct, geometry_pylon, geometry_support, geometry_control,
                 geometry_nacelle, geometry_propeller, geometry_ht, geometry_vt, aircraft_type: str, power: str,
                 bem_input, delta_e: float, delta_r: float, wing_span: float, wing_sweep: float, wing_airfoil: str,
                 wing_tr: float, wing_cr: float, l_coc: float, l_cab: float, l_tail: float, fus_d: float, pax: float,
                 propulsor_type: str, x_wing: float, x_duct: float, comp_pe, rpm: float, z_duct: float,
                 cv_mode: str):
        self.rpm = rpm
        self.aircraft_type = aircraft_type
        self.comp_pe = comp_pe
        self.pc = power
        self.bem_input = bem_input
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
        self.propulsor_type = propulsor_type
        self.x_wing = x_wing
        self.x_duct = x_duct
        self.z_duct = z_duct
        self.conditions = conditions
        self.reference = reference
        self.geometry_duct = geometry_duct
        self.geometry_pylon = geometry_pylon
        self.geometry_support = geometry_support
        self.geometry_control = geometry_control
        self.geometry_nacelle = geometry_nacelle
        self.geometry_propeller = geometry_propeller
        self.geometry_ht = geometry_ht
        self.geometry_vt = geometry_vt
        self.cv_mode = cv_mode

        self.wing = Wing(b_wing=self.wing_span, sweep_wing=self.wing_sweep, wing_airfoil=self.wing_airfoil,
                         taper_ratio_wing=self.wing_tr, c_root_wing=self.wing_cr, alpha=self.conditions[1],
                         velocity=self.conditions[0], mach=self.conditions[3], wing_type=self.aircraft_type,
                         altitude=self.conditions[2])

        self.fuselage = Fuselage(fuselage_length=self.l_fuse, fuselage_diameter=self.fus_d, l_cabin=self.l_cab,
                                 l_cockpit=self.l_coc, l_tail=self.l_tail, velocity=self.conditions[0],
                                 alpha=self.conditions[1], ref_area=self.reference[0], mach=self.conditions[3],
                                 cl_wing=self.wing.cl_prime(), cmac_wing=self.wing_cr, altitude=self.conditions[2])

        self.lg = LandingGear(aircraft_type=self.aircraft_type)

        self.operation = Operations(MTOM=ref.MTOW, pax=self.pax)

        self.empennage = self.make_empennage()

    def make_empennage(self):
        ar_wing = self.wing.aspect_ratio()
        cl_wing = self.wing.cl_prime()
        cla_wing = self.wing.cl_al()

        if self.aircraft_type == "DUUC":
            from aircraft.propulsive_empennage.empennage_assembly_PE import PropulsiveEmpennage
            va_inlet = self.conditions[0] * (np.pi * (self.geometry_duct[0] / 2))
            return PropulsiveEmpennage(rpm=self.rpm,
                                       power_condition=self.pc,
                                       va_inlet=va_inlet,
                                       d_exit=area_ratio("0016", self.geometry_duct[1], self.geometry_duct[0]/2, 1)[0],
                                       propulsor_type=self.propulsor_type,
                                       ar_wing=ar_wing,
                                       cl_wing=cl_wing,
                                       cla_wing=cla_wing,
                                       bem_input=self.bem_input,
                                       delta_e=self.delta_e,
                                       delta_r=self.delta_r,
                                       conditions=self.conditions,
                                       reference=self.reference,
                                       geometry_duct=self.geometry_duct,
                                       geometry_pylon=self.geometry_pylon,
                                       geometry_control=self.geometry_control,
                                       geometry_nacelle=self.geometry_nacelle,
                                       geometry_support=self.geometry_support,
                                       geometry_propeller=self.geometry_propeller,
                                       comp_pe=self.comp_pe,
                                       cv_mode=self.cv_mode)
        if self.aircraft_type == "conventional":
            from aircraft.conventional_empennage.empennage_assembly_conv import ConventionalEmpennage
            return ConventionalEmpennage(rpm=config.rpm,
                                         power_condition=self.pc,
                                         ar_wing=ar_wing,
                                         cl_wing=cl_wing,
                                         cla_wing=cla_wing,
                                         bem_input=self.bem_input,
                                         conditions=self.conditions,
                                         reference=self.reference,
                                         ht_geometry=self.geometry_ht,
                                         vt_geometry=self.geometry_vt,
                                         nacelle_geometry=self.geometry_nacelle,
                                         geometry_prop=self.geometry_propeller)
        else:
            print("Wrong aircraft type defined!!")
            return None

    def h_lever_PE(self):
        h_lever = self.x_duct - self.x_cog()[0]
        return h_lever

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
                                  c_mac_wing=ref.c_mac_w,
                                  x_wing=0,
                                  x_duct=0,
                                  w_fuel=0,
                                  w_pax=0,
                                  z_PE=0)
            return cog.x_cg()[0], cog.x_cg()[1], cog.cg_fuselage_group()[0], cog.cg_wing_group()[0], cog.z_cg()
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
                                  c_mac_wing=ref.c_mac_w,
                                  x_wing=self.x_wing,
                                  x_duct=self.x_duct,
                                  w_fuel=0,
                                  w_pax=0,
                                  z_PE=self.z_duct + 0.5 * self.geometry_duct[0])
            return cog.x_cg()[0], cog.x_cg()[1], cog.cg_fuselage_group()[0], cog.cg_wing_group()[0], cog.z_cg()
        else:
            return 0, 0

    def x_cog_ferry(self):
        length_fuselage = self.l_tail + self.l_coc + self.l_cab
        if self.aircraft_type == "conventional":
            cog_start = CenterOfGravity(w_lg_nose=self.lg.weight()[0], w_lg_main=self.lg.weight()[1],
                                        w_eng_wing=self.empennage.weight_prop() / 2, w_wing=self.wing.weight(),
                                        w_fuse=self.fuselage.weight(), w_nac_wing=self.empennage.weight_nac() / 2,
                                        w_vt=self.empennage.vt_tail.weight(), w_ht=self.empennage.ht_tail.weight(),
                                        w_duct=0, w_sys=self.operation.weight_sys(), l_fuselage=length_fuselage,
                                        aircraft_type=self.aircraft_type, c_mac_wing=ref.c_mac_w,
                                        x_wing=0, x_duct=0, w_fuel=config.w_fuel_ferry_start,  w_pax=0, z_PE=0)
            vector_start = [cog_start.x_cg()[0], cog_start.x_cg()[1], cog_start.cg_fuselage_group()[0],
                            cog_start.cg_wing_group()[0]]

            cog_end = CenterOfGravity(w_lg_nose=self.lg.weight()[0], w_lg_main=self.lg.weight()[1],
                                      w_eng_wing=self.empennage.weight_prop() / 2, w_wing=self.wing.weight(),
                                      w_fuse=self.fuselage.weight(), w_nac_wing=self.empennage.weight_nac() / 2,
                                      w_vt=self.empennage.vt_tail.weight(), w_ht=self.empennage.ht_tail.weight(),
                                      w_duct=0, w_sys=self.operation.weight_sys(), l_fuselage=length_fuselage,
                                      aircraft_type=self.aircraft_type, c_mac_wing=ref.c_mac_w,
                                      x_wing=0, x_duct=0, w_fuel=config.w_fuel_ferry_end,  w_pax=0, z_PE=0)
            vector_end = [cog_end.x_cg()[0], cog_end.x_cg()[1], cog_end.cg_fuselage_group()[0],
                          cog_end.cg_wing_group()[0]]
            return vector_start, vector_end
        if self.aircraft_type == "DUUC":
            cog_start = CenterOfGravity(w_lg_nose=self.lg.weight()[0], w_lg_main=self.lg.weight()[1], w_eng_wing=0,
                                        w_wing=self.wing.weight(), w_fuse=self.fuselage.weight(), w_nac_wing=0, w_vt=0,
                                        w_ht=0, w_duct=self.empennage.weight() * 2, w_sys=self.operation.weight_sys(),
                                        l_fuselage=length_fuselage, aircraft_type=self.aircraft_type,
                                        c_mac_wing=ref.c_mac_w, x_wing=self.x_wing, x_duct=self.x_duct,
                                        w_fuel=config.w_fuel_ferry_start, w_pax=0,
                                        z_PE=self.z_duct + 0.5 * self.geometry_duct[0])
            vector_start = [cog_start.x_cg()[0], cog_start.x_cg()[1], cog_start.cg_fuselage_group()[0],
                            cog_start.cg_wing_group()[0]]

            cog_end = CenterOfGravity(w_lg_nose=self.lg.weight()[0], w_lg_main=self.lg.weight()[1], w_eng_wing=0,
                                      w_wing=self.wing.weight(), w_fuse=self.fuselage.weight(), w_nac_wing=0, w_vt=0,
                                      w_ht=0, w_duct=self.empennage.weight() * 2, w_sys=self.operation.weight_sys(),
                                      l_fuselage=length_fuselage, aircraft_type=self.aircraft_type,
                                      c_mac_wing=ref.c_mac_w, x_wing=self.x_wing, x_duct=self.x_duct,
                                      w_fuel=config.w_fuel_ferry_end, w_pax=0,
                                      z_PE=self.z_duct + 0.5 * self.geometry_duct[0])
            vector_end = [cog_end.x_cg()[0], cog_end.x_cg()[1], cog_end.cg_fuselage_group()[0],
                          cog_end.cg_wing_group()[0]]
            return vector_start, vector_end
        else:
            return 0, 0

    def x_cog_payload(self):
        length_fuselage = self.l_tail + self.l_coc + self.l_cab
        if self.aircraft_type == "conventional":
            cog_start = CenterOfGravity(w_lg_nose=self.lg.weight()[0], w_lg_main=self.lg.weight()[1],
                                        w_eng_wing=self.empennage.weight_prop() / 2, w_wing=self.wing.weight(),
                                        w_fuse=self.fuselage.weight(), w_nac_wing=self.empennage.weight_nac() / 2,
                                        w_vt=self.empennage.vt_tail.weight(), w_ht=self.empennage.ht_tail.weight(),
                                        w_duct=0, w_sys=self.operation.weight_sys(), l_fuselage=length_fuselage,
                                        aircraft_type=self.aircraft_type, c_mac_wing=ref.c_mac_w,
                                        x_wing=0, x_duct=0, w_fuel=config.w_fuel_full_start,
                                        w_pax=self.operation.weight_pax(), z_PE=0)
            vector_start = [cog_start.x_cg()[0], cog_start.x_cg()[1], cog_start.cg_fuselage_group()[0],
                            cog_start.cg_wing_group()[0]]

            cog_end = CenterOfGravity(w_lg_nose=self.lg.weight()[0], w_lg_main=self.lg.weight()[1],
                                      w_eng_wing=self.empennage.weight_prop() / 2, w_wing=self.wing.weight(),
                                      w_fuse=self.fuselage.weight(), w_nac_wing=self.empennage.weight_nac() / 2,
                                      w_vt=self.empennage.vt_tail.weight(), w_ht=self.empennage.ht_tail.weight(),
                                      w_duct=0, w_sys=self.operation.weight_sys(), l_fuselage=length_fuselage,
                                      aircraft_type=self.aircraft_type, c_mac_wing=ref.c_mac_w,
                                      x_wing=0, x_duct=0, w_fuel=config.w_fuel_full_end,
                                      w_pax=self.operation.weight_pax(), z_PE=0)
            vector_end = [cog_end.x_cg()[0], cog_end.x_cg()[1], cog_end.cg_fuselage_group()[0],
                          cog_end.cg_wing_group()[0]]
            return vector_start, vector_end
        if self.aircraft_type == "DUUC":
            cog_start = CenterOfGravity(w_lg_nose=self.lg.weight()[0], w_lg_main=self.lg.weight()[1], w_eng_wing=0,
                                        w_wing=self.wing.weight(), w_fuse=self.fuselage.weight(), w_nac_wing=0, w_vt=0,
                                        w_ht=0, w_duct=self.empennage.weight() * 2, w_sys=self.operation.weight_sys(),
                                        l_fuselage=length_fuselage, aircraft_type=self.aircraft_type,
                                        c_mac_wing=ref.c_mac_w, x_wing=self.x_wing, x_duct=self.x_duct,
                                        w_fuel=config.w_fuel_full_start, w_pax=self.operation.weight_pax(),
                                        z_PE=self.z_duct + 0.5 * self.geometry_duct[0])
            vector_start = [cog_start.x_cg()[0], cog_start.x_cg()[1], cog_start.cg_fuselage_group()[0],
                            cog_start.cg_wing_group()[0]]

            cog_end = CenterOfGravity(w_lg_nose=self.lg.weight()[0], w_lg_main=self.lg.weight()[1], w_eng_wing=0,
                                      w_wing=self.wing.weight(), w_fuse=self.fuselage.weight(), w_nac_wing=0, w_vt=0,
                                      w_ht=0, w_duct=self.empennage.weight() * 2, w_sys=self.operation.weight_sys(),
                                      l_fuselage=length_fuselage, aircraft_type=self.aircraft_type,
                                      c_mac_wing=ref.c_mac_w, x_wing=self.x_wing, x_duct=self.x_duct,
                                      w_fuel=config.w_fuel_full_end, w_pax=self.operation.weight_pax(),
                                      z_PE=self.z_duct + 0.5 * self.geometry_duct[0])
            vector_end = [cog_end.x_cg()[0], cog_end.x_cg()[1], cog_end.cg_fuselage_group()[0],
                          cog_end.cg_wing_group()[0]]
            return vector_start, vector_end
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
            return [cd0_duct, cd0_pylon, cd0_support, cd0_control]

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
            cl_empennage = self.empennage.cl_norm_vector()

            cl_tot = cl_fuselage + cl_wing + cl_empennage
            return cl_tot

    def cm(self):
        x_cog = self.x_cog()[0]
        z_cog = self.x_cog()[4]
        x_ac = self.x_wing + 0.25 * ref.c_mac_w
        q = 0.5 * self.conditions[0] ** 2 * air_density_isa(self.conditions[2])[0] * ref.s_w
        if self.aircraft_type == "conventional":
            l_h = ref.lever_h
            z_e = self.fus_d

            cm_sums = self.empennage.cm_emp() + self.fuselage.cm() + self.wing.cm()
            cm_l_wing = self.wing.cn() * q * (x_cog - x_ac)
            cm_l_emp = self.empennage.cn() * q * l_h
            cm_t_eng = self.empennage.thrust() * (z_e - z_cog)

            cm_ac = cm_sums + (cm_l_wing + cm_l_emp - cm_t_eng) / q
            #print(f"Cm aircraft (ATR): {cm_ac}")
            return cm_ac
        elif self.aircraft_type == "DUUC":
            l_h = self.x_duct
            z_e = self.z_duct + 0.5 * self.geometry_duct[0]
            cm_sums = self.empennage.cm_emp() + self.fuselage.cm() + self.wing.cm()

            cm_l_wing = self.wing.cn() * q * (x_cog - x_ac)
            cm_l_emp = self.empennage.cn() * q * (l_h - x_cog)
            cm_t_eng = self.empennage.thrust() * (z_e - z_cog)

            cm_ac = cm_sums + (cm_l_wing + cm_l_emp - cm_t_eng) / q
            #print(f"Cm aircraft (DUUC): {cm_ac}")
            return cm_ac
        else:
            raise ValueError("Unknown aircraft type specified")

    def tc(self):
        if self.aircraft_type == "conventional":
            ct_ac = 0
        elif self.aircraft_type == "DUUC":
            ct_ac = 1.2
        else:
            raise ValueError("Aircraft type not specified properly")
        return ct_ac

    """ ----------------------------------------------- Derrivatives ----------------------------------------------- """
    # ---------------------------- coefficient from Propulsive Empennage ---------
    def vt_volume_coefficient_PE(self):
        vt_coeff_pe = tail_volume(2 * self.empennage.duct.proj_area(), self.h_lever_PE(), self.reference[0], ref.b_w)
        return vt_coeff_pe

    def cm_a_PE(self):
        cm_a_pe = self.dcn_da_pe() * 0.95 * (self.h_lever_PE() / ref.c_mac_w)
        return cm_a_pe

    def cn_beta_PE(self):
        cn_beta_pe = (self.vt_volume_coefficient_PE() * self.empennage.duct.cl_a() * 0.90) / (1 + ref.side_wash)
        return cn_beta_pe

    def cm_de_PE(self):
        cm_de_pe = (- 0.9 * tail_volume(4 * self.empennage.elevator.area(), self.h_lever_PE(), ref.b_w, ref.c_mac_w)
                    * self.empennage.elevator.cl_a())
        return cm_de_pe

    def dcn_da_pe(self):
        from aircraft.propulsive_empennage.empennage_assembly_PE import PropulsiveEmpennage
        ar_wing = self.wing.aspect_ratio()
        cl_wing = self.wing.cl_prime()
        cla_wing = self.wing.cl_al()
        va_inlet = self.conditions[0] * (np.pi * (self.geometry_duct[0] / 2))
        conditions0 = [128, 10, 7000, 0.44, 0, 0, 0]
        cn10: PropulsiveEmpennage = PropulsiveEmpennage(rpm=self.rpm,
                                                        power_condition=self.pc,
                                                        va_inlet=va_inlet,
                                                        d_exit=area_ratio("0016", self.geometry_duct[1],
                                                                         self.geometry_duct[0] / 2, 1)[0],
                                                        propulsor_type=self.propulsor_type,
                                                        ar_wing=ar_wing,
                                                        cl_wing=cl_wing,
                                                        cla_wing=cla_wing,
                                                        bem_input=self.bem_input,
                                                        delta_e=0,
                                                        delta_r=0,
                                                        conditions=conditions0,
                                                        reference=self.reference,
                                                        geometry_duct=self.geometry_duct,
                                                        geometry_pylon=self.geometry_pylon,
                                                        geometry_control=self.geometry_control,
                                                        geometry_nacelle=self.geometry_nacelle,
                                                        geometry_support=self.geometry_support,
                                                        geometry_propeller=self.geometry_propeller,
                                                        comp_pe=self.comp_pe, cv_mode=self.cv_mode)
        cn_a10 = cn10.cn()
        conditions1 = [128, 0, 7000, 0.44, 0, 0, 0]
        cn0: PropulsiveEmpennage = PropulsiveEmpennage(rpm=self.rpm,
                                                       power_condition=self.pc,
                                                       va_inlet=va_inlet,
                                                       d_exit=area_ratio("0016", self.geometry_duct[1],
                                                                         self.geometry_duct[0] / 2, 1)[0],
                                                       propulsor_type=self.propulsor_type,
                                                       ar_wing=ar_wing,
                                                       cl_wing=cl_wing,
                                                       cla_wing=cla_wing,
                                                       bem_input=self.bem_input,
                                                       delta_e=0,
                                                       delta_r=0,
                                                       conditions=conditions1,
                                                       reference=self.reference,
                                                       geometry_duct=self.geometry_duct,
                                                       geometry_pylon=self.geometry_pylon,
                                                       geometry_control=self.geometry_control,
                                                       geometry_nacelle=self.geometry_nacelle,
                                                       geometry_support=self.geometry_support,
                                                       geometry_propeller=self.geometry_propeller,
                                                       comp_pe=self.comp_pe, cv_mode=self.cv_mode)
        cn_a0 = cn0.cn()

        dcn = (cn_a10 - cn_a0) / 10
        dcn_rad = dcn * 180 / np.pi
        return dcn_rad

    # ---------------------------- coefficients due to alpha ----------------------------
    def cl_a(self):
        cl_a_ac = self.wing.cl_al() + self.empennage.cl_a()
        return cl_a_ac

    def cm_a(self):
        if self.aircraft_type == "conventional":
            cm_a_ac = -1.6677
        elif self.aircraft_type == "DUUC":
            cm_a_fus = 0.05
            cm_a_wing = 0.17
            cm_a_pe = self.dcn_da_pe() * 0.95 * (self.h_lever_PE() / ref.c_mac_w)
            cm_a_ac = cm_a_fus + cm_a_wing - cm_a_pe
        else:
            raise ValueError()
        return cm_a_ac

    # ---------------------------- coefficient due to elevator deflection ----------------------------
    def cl_de(self):
        cl_de_ac = self.empennage.cl_de()
        return cl_de_ac

    def cm_de(self):
        if self.aircraft_type == "conventional":
            cm_de_ac = self.empennage.cm_de()
        elif self.aircraft_type == "DUUC":
            cm_de_ac = self.cm_de_PE()
        else:
            raise ValueError()
        return cm_de_ac

    # ---------------------------- coefficients due to rudder deflection ----------------------------
    def cn_dr(self):
        if self.aircraft_type == "conventional":
            cn_dr_ac = -0.3273
            return cn_dr_ac
        elif self.aircraft_type == "DUUC":
            l_v = self.x_duct - self.x_cog()[0]
            cn_dr_ac = - self.cy_dr() * (l_v / ref.b_w)
            return cn_dr_ac
        else:
            raise ValueError

    def cy_dr(self):
        cy_dr_ac = self.empennage.cy_dr()
        return cy_dr_ac

    # ---------------------------- coefficients due to sideslip ----------------------------
    def cy_beta(self):
        cy_beta_ac = self.empennage.cy_beta()
        return cy_beta_ac

    def cn_beta(self):
        """ wing contribution neglected as the effect is small due to small sweep """
        if self.aircraft_type == "conventional":
            """ this includes fuselage contribution as the value is from the full aircraft """
            cn_beta_conv = self.empennage.cn_beta()
            return cn_beta_conv
        elif self.aircraft_type == "DUUC":
            cn_beta_emp = self.cn_beta_PE()
            cn_beta_fus = - 0.08
            cn_beta_duuc = cn_beta_emp + cn_beta_fus
            return cn_beta_duuc
        else:
            raise ValueError()

    """ ----------------------------------------------- FORCES ----------------------------------------------------- """
    def lift(self):
        lift_wing = self.wing.lift()
        lift_fuselage = self.fuselage.lift()
        if self.aircraft_type == "conventional":
            lift_empennage = self.empennage.lift()
        elif self.aircraft_type == "DUUC":
            lift_empennage = self.empennage.lift()
        else:
            raise ValueError("Aircraft type not specified properly")
        lift_ac = lift_wing + lift_empennage + lift_fuselage
        return lift_ac

    def drag(self):
        drag_wing = self.wing.drag()
        drag_fuselage = self.fuselage.drag()
        if self.aircraft_type == "conventional":
            drag_empennage = self.empennage.drag()
        elif self.aircraft_type == "DUUC":
            drag_empennage = self.empennage.drag()
        else:
            raise ValueError("Aircraft type not specified properly")
        drag_ac = drag_wing + drag_empennage + drag_fuselage
        return drag_ac

    def thrust(self):
        if self.aircraft_type == "conventional":
            thrust_eng = self.bem_input[0]
            thrust_duct = 0
        elif self.aircraft_type == "DUUC":
            thrust_eng = self.bem_input[0]
            thrust_duct = self.empennage.duct.thrust()[0] * 2
        else:
            raise ValueError("Aircraft type not specified properly")
        thrust_ac = thrust_eng + thrust_duct
        return thrust_ac

    def weight(self):
        weight_fuselage = self.fuselage.weight()
        weight_wing = self.wing.weight()
        weight_landing = self.lg.weight()[2]
        weight_systems = self.operation.weight_sys() + self.operation.weight_pax() + self.operation.weight_pax()

        if self.aircraft_type == "conventional":
            weight_empennage = self.empennage.weight()
        elif self.aircraft_type == "DUUC":
            weight_empennage = self.empennage.weight() * 2
        else:
            raise ValueError("Aircraft type not specified properly")

        weight_aircraft = weight_empennage + weight_systems + weight_wing + weight_fuselage + weight_landing
        return weight_aircraft

    def slope_htail(self):
        cog = self.x_cog()[0]
        x_lemac = self.x_cog()[1]

        print(f"cog: {cog}, xlemac: {x_lemac}")
        print(f"diff: {cog - x_lemac}")

        cl_h = self.empennage.cl_norm_vector()
        cl = self.cl_ac()

        cl_h = 4.586
        cl = 1.44

        a1 = slopes("control", self.aircraft_type, ref.z_h, ref.phi_qc_w, self.wing.aspect_ratio(),
                    self.fuselage.length(), x_lemac, 2.234, 0.9, cl_h, cl, ref.tr_w, ref.b_w, 0,
                    self.conditions[3], ref.ar_h, 0, 0)[0]

        b1 = slopes("control", self.aircraft_type, ref.z_h, ref.phi_qc_w, self.wing.aspect_ratio(),
                    self.fuselage.length(), x_lemac, 2.234, 0.9, cl_h, cl, ref.tr_w, ref.b_w, 0,
                    self.conditions[3], ref.ar_h, 0, 0)[1]

        a2 = slopes("stability", self.aircraft_type, ref.z_h, ref.phi_qc_w, self.wing.aspect_ratio(),
                    self.fuselage.length(), x_lemac, 2.234, 0.9, cl_h, cl, ref.tr_w, ref.b_w, 0,
                    self.conditions[3], ref.ar_h, 0, 0)[0]

        return a1, b1, a2
