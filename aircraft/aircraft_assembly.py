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
    def __init__(self, conditions, reference, geometry_duct, geometry_pylon, geometry_support, geometry_control,
                 geometry_nacelle, geometry_propeller, geometry_ht, geometry_vt, aircraft_type: str, power: str,
                 bem_input, delta_e: float, delta_r: float, wing_span: float, wing_sweep: float, wing_airfoil: str,
                 wing_tr: float, wing_cr: float, l_coc: float, l_cab: float, l_tail: float, fus_d: float, pax: float,
                 propulsor_type: str, x_wing: float, x_duct: float, comp_pe, rpm: float):
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
                                       comp_pe=self.comp_pe)
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
                                  x_duct=0)
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
                                  c_mac_wing=ref.c_mac_w,
                                  x_wing=self.x_wing,
                                  x_duct=self.x_duct)
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

    def ct(self):
        if self.aircraft_type == "conventional":
            ct_ac = 0
        elif self.aircraft_type == "DUUC":
            ct_ac = 1.2
        else:
            raise ValueError("Aircraft type not specified properly")
        return ct_ac

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
            thrust_eng = 0
            thrust_duct = 0
        elif self.aircraft_type == "DUUC":
            thrust_eng = 0
            thrust_duct = 0
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

    def thrust(self):
        t_ac = 0
        return t_ac
