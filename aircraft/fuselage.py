import numpy as np
import data.atr_reference as ref
import matplotlib.pyplot as plt
import config
from analysis_modules.factors import skin_friction, mach_correction
from analysis_modules.ISA import air_density_isa
from analysis_modules.aerodynamic import reynolds
import flow_conditions


class Fuselage:
    def __init__(self, fuselage_length: float, fuselage_diameter: float, l_cabin: float, l_cockpit: float,
                 l_tail: float, velocity: float, alpha: float, ref_area: float, mach: float,
                 cl_wing: float, cmac_wing: float, altitude:float):
        super().__init__()
        self.fuselage_length = fuselage_length
        self.fuselage_diameter = fuselage_diameter
        self.l_cabin = l_cabin
        self.l_cockpit = l_cockpit
        self.l_tail = l_tail
        self.velocity = velocity
        self.alpha = alpha
        self.ref_area = ref_area
        self.mach = mach
        self.cl_wing = cl_wing
        self.c_wing = cmac_wing
        self.altitude = altitude

    """ ------------------------- Determine inflow properties ----------------------------------------------------- """
    def inflow_velocity(self):
        """ flow is undisturbed in the freestream, returned in m/s"""
        vel_fus = self.velocity
        return vel_fus

    def inflow_angle(self):
        """ flow is undisturbed in the free stream, returned in radians"""
        ang_fus = self.alpha
        return ang_fus

    def rey_fuselage(self):
        re_fus = reynolds(air_density_isa(self.altitude), self.inflow_velocity(), self.fuselage_length)
        return re_fus

    """ ---------------------------- Determine geometric properties ---------------------------------------------- """
    @staticmethod
    def l_cabin():
        k_cabin = 1
        length_cabin = k_cabin * (config.n_pax / 4)

        return np.round(length_cabin, 0)

    def length(self):
        """ 4 m is cockpit length, 1.4 * D is tail length"""
        l_fus = self.l_cabin + 1.4 * self.fuselage_diameter + 4
        return l_fus

    def area_proj(self):
        area_proj_fus = self.length() * self.fuselage_diameter
        return area_proj_fus

    def area_front(self):
        area_front_fus = (np.pi / 4) * self.fuselage_diameter ** 2
        return area_front_fus

    def slenderness(self):
        slenderness_fus = self.length() / self.fuselage_diameter
        return slenderness_fus

    def aspect_ratio(self):
        ar_fus = (4 * self.length() ** 2) / (np.pi * self.area_front())

        ar_fus = (np.pi / 4) * self.fuselage_length/self.fuselage_diameter
        return ar_fus

    def area_wetted(self):
        """ Based on Sadreay wetted area fuselage"""
        radius = self.fuselage_diameter / 2
        area_wet_fus = 2 * np.pi * radius ** 2 + 2 * (np.pi * radius ** 2 * self.fuselage_length) / radius

        """ based on TextNita"""
        l_bug = 1.4 * self.fuselage_diameter
        l_heck = 3 * self.fuselage_diameter
        l_zyl = self.l_cabin
        d = self.fuselage_diameter

        area_wet = (np.sqrt(l_bug ** 2 + (d / 2) ** 2) * (d * np.pi) / 2
                    + (self.fuselage_diameter * np.pi * l_zyl) + np.sqrt(l_heck ** 2 + (d / 2) ** 2)
                    * (d * np.pi) / 2)

        return area_wet

    """ --------------------------------- Determine coefficients ------------------------------------------------- """
    def cd0(self):
        """ Based on Sadraey fuselage model"""
        cf = skin_friction(self.re, "t")
        fm = mach_correction(self.mach)
        f_fus = (1 + 60 / (self.fuselage_length / self.fuselage_diameter) ** 3 + 0.0025
                 * self.fuselage_length / self.fuselage_diameter)
        norm_area = self.area_wetted() / self.ref_area

        cd0_fus = cf * fm * f_fus * norm_area
        return cd0_fus

    def e_fuselage(self):
        """ Based on Raymer's approximation"""
        e_fus = 1 / (1.05 + 0.38 * self.cd0() * np.pi * self.aspect_ratio())
        return e_fus

    def cdi(self):
        cdi_fus = self.cl() ** 2 / (np.pi * self.aspect_ratio() * self.e_fuselage())
        return cdi_fus

    def cd(self):
        cd_fus = self.cdi() + self.cd0()
        return cd_fus

    def cl(self):
        cl_fus = (2 / 3) * (self.area_front() / self.ref_area) * self.cl_wing
        return cl_fus

    @staticmethod
    def cd_interference():
        """ interference drag because of the wing fuselage interaction"""
        dCl_dc = -0.07 * ref.alpha_install_wing

        cd_int_fuse = 0.015 * (dCl_dc ** 2)
        return cd_int_fuse

    def cm0(self):
        """ Defined by Hoerners 1965"""
        cm_fuselage = - self.length() / (4 * self.fuselage_diameter)
        return cm_fuselage

    def cmi(self):
        """ Defined by Hoerners 1965"""
        lr = self.length() / self.fuselage_diameter

        cm_fus = - 1.2 / lr
        return cm_fus

    def cm(self):
        cm_fuselage = self.cm0() + self.cmi() * self.alpha
        return cm_fuselage

    """ -------------------------------- Determine prime outputs -------------------------------------------------- """
    def cd_prime(self):
        norm_speed = self.inflow_velocity() ** 2 / self.velocity ** 2
        norm_area = self.area_proj() / self.ref_area

        alpha = np.radians(self.inflow_angle())

        cd_cd = self.cd() * np.cos(alpha) * norm_area * norm_speed

        cd_cl = self.cl() * np.sin(alpha) * norm_area * norm_speed

        cd_prime_fus = cd_cd + cd_cl
        return cd_prime_fus

    def cl_prime(self):
        norm_speed = self.inflow_velocity() ** 2 / self.velocity ** 2
        norm_area = self.area_proj() / self.ref_area

        alpha = np.radians(self.inflow_angle())

        cl_cl = self.cl() * np.cos(alpha) * norm_speed * norm_area

        cl_cd = self.cd() * np.sin(alpha) * norm_speed * norm_area

        cl_prime_fus = cl_cl + cl_cd
        return cl_prime_fus

    def cm_prime(self):
        norm_speed = self.inflow_velocity() ** 2 / self.velocity
        norm_area = self.area_proj() / self.ref_area

        cm = self.cm() * norm_area * norm_speed * self.c_wing
        return cm

    """  --------------------------------------- Determine the forces of the fuselage ---------------------------- """
    def lift(self):
        lift_fuselage = self.cl_prime() * self.inflow_velocity() ** 2 * self.area_proj() * 0.5 * flow_conditions.rho
        return lift_fuselage

    def drag(self):
        drag_fuselage = self.cd_prime() * self.inflow_velocity() ** 2 * self.area_proj() * 0.5 * flow_conditions.rho
        return drag_fuselage

    """  --------------------------------------- Determine the weight of the fuselage ---------------------------- """
    def weight(self):
        vd = ref.v_dive
        hf = self.fuselage_diameter
        wf = self.fuselage_diameter

        lh = ref.lever_h  # lever arm between fuselage and tail

        w_fuselage = 0.23 * np.sqrt(vd * lh / (wf + hf)) * self.area_wetted() ** 1.2
        return w_fuselage


""" Test section """

if __name__ == "__main__":
    fuselage = Fuselage(fuselage_length=ref.l_tail+ref.l_cabin+ref.l_cockpit,
                        fuselage_diameter=2.77,
                        l_cabin=ref.l_cabin,
                        l_cockpit=ref.l_cockpit,
                        l_tail=ref.l_tail,
                        velocity=128,
                        alpha=0,
                        ref_area=ref.s_w,
                        mach=0.4,
                        cl_wing=0.813,
                        cmac_wing=2.626,
                        altitude=7000)

    print(f"inflow vel: {fuselage.inflow_velocity()}")
    print(f"inflow ang: {fuselage.inflow_angle()}, oswald: {fuselage.e_fuselage()}")
    print(f"area: {fuselage.area_proj()}, wetted area: {fuselage.area_wetted()}")
    print(f"aspect ratio: {fuselage.aspect_ratio()}")
    print(f"cd0: {fuselage.cd0()}, cdi: {fuselage.cdi()}, cd: {fuselage.cd()}, cdprime: {fuselage.cd_prime()}")
    print(f" cl: {fuselage.cl()}, cl_prime: {fuselage.cl_prime()}")
    print(f"cm0: {fuselage.cm0()}, cma: {fuselage.cmi()}, cm: {fuselage.cm()}")
    print(f"weight: {fuselage.weight()}")
    print(f"Slenderness: {fuselage.slenderness()}")
