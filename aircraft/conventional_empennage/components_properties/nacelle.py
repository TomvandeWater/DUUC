import numpy as np
from analysis_modules.aerodynamic import drag_interference
from analysis_modules.factors import mach_correction, skin_friction


class Nacelle:
    """ Nacelle class for conventional configuration"""
    def __init__(self, nacelle_length: float, nacelle_diameter: float, v_inf: float, alpha: float,
                 area_ref: float, prop_airfoil: str, n_blades: float, mach: float, reynolds: float):
        super().__init__()
        self.nacelle_length = nacelle_length
        self.nacelle_diameter = nacelle_diameter
        self.v_inf = v_inf
        self.alpha = alpha
        self.area_ref = area_ref
        self.prop_airfoil = prop_airfoil
        self.n_blades = n_blades
        self.mach = mach
        self.re_ref = reynolds

    def inflow_velocity(self):
        inflow_nac = self.v_inf
        return inflow_nac

    def area_wetted(self):
        """ only the side area is assumed no closing sides"""
        area_cylinder = np.pi * self.nacelle_diameter * self.nacelle_length

        """ half a sphere is used to close of the rear end"""
        area_rear = 0.5 * np.pi * self.nacelle_diameter ** 2

        area_wet_nac = area_rear + area_cylinder
        return area_wet_nac

    def area(self):
        area_nac = self.nacelle_length * self.nacelle_diameter
        return area_nac

    def t_c(self):
        num_list = [int(digit) for digit in self.prop_airfoil]
        thickness = num_list[2] * 10 + num_list[3]  # NACA thickness of profile
        thickness = thickness / 100  # returns value in percentage of normalized chord
        return thickness

    def cd_int(self):
        norm_area = self.area() / self.area_ref
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2

        cd_int_nac = (drag_interference(self.t_c(), 'plane') * norm_speed * norm_area
                      * self.n_blades)
        return cd_int_nac

    def cd0(self):
        cf = skin_friction(self.re_ref, "t")
        fm = mach_correction(self.mach)
        l_d = self.nacelle_length / self.nacelle_diameter
        f_nac = 1 + 60 / l_d ** 3 + 0.0025 * l_d
        norm_area = self.area_wetted() / self.area_ref

        cd0_nacelle = cf * fm * f_nac * norm_area

        return cd0_nacelle

    def cd_prime(self):
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2

        cd_nac = self.cd0() * norm_speed
        return cd_nac

    @staticmethod
    def cl_prime():
        """assume the nacelle does not produce lift"""
        cl_prime_nac = 0
        return cl_prime_nac


