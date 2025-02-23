import numpy as np
import data.atr_reference as ref
import matplotlib.pyplot as plt
import config
from analysis_modules.factors import skin_friction, mach_correction


class Fuselage:
    def __init__(self, fuselage_length: float, fuselage_diameter: float, l_cabin: float, l_cockpit: float, l_tail: float,
                 velocity: float, alpha: float, ref_area: float, reynolds_number: float, mach: float, cl_wing: float):
        super().__init__()
        self.fuselage_length = fuselage_length
        self.fuselage_diameter = fuselage_diameter
        self.l_cabin = l_cabin
        self.l_cockpit = l_cockpit
        self.l_tail = l_tail
        self.velocity = velocity
        self.alpha = alpha
        self.ref_area = ref_area
        self.re = reynolds_number
        self.mach = mach
        self.cl_wing = cl_wing
    """ Determine inflow velocity and angle"""
    def inflow_velocity(self):
        """ flow is undisturbed in the freestream, returned in m/s"""
        vel_fus = self.inflow_velocity()
        return vel_fus

    def inflow_angle(self):
        """ flow is undisturbed in the free stream, returned in radians"""
        ang_fus = self.alpha
        return ang_fus

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

    def aspect_ratio(self):
        ar_fus = (4 * self.length() ** 2) / (np.pi * self.area_front())
        return ar_fus

    def area_wetted(self):
        """ Based on Sadreay wetted area fuselage"""
        radius = self.fuselage_diameter / 2
        area_wet_fus = 2 * np.pi * radius ** 2 + 2 * (np.pi * radius ** 2 * self.fuselage_length) / radius
        return area_wet_fus

    """ Determine coefficients """
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

    @staticmethod
    def cd_interference():
        """ interference drag because of the wing fuselage interaction"""
        dCl_dc = -0.07 * ref.alpha_install_wing

        cd_int_fuse = 0.015 * (dCl_dc ** 2)
        return cd_int_fuse

    """ Determine the weight of the fuselage"""
    def weight(self):
        w_fuselage = config.n_pax * config.w_pax + 15000 # still add W structure
        return  w_fuselage


