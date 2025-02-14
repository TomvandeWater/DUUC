import numpy as np
import constants
import flow_conditions
from analysis_modules.factors import skin_friction, mach_correction

class Nacelle:

    def __init__(self, nacelle_length: float, nacelle_diameter: float,
                 propulsor_type: str, re_induct: float, power_condition: str, u_mom: float,
                 alpha: float, ref_area: float):
        super().__init__()
        self.nacelle_length = nacelle_length
        self.nacelle_diameter = nacelle_diameter
        self.propulsor_type = propulsor_type
        self.re_induct = re_induct
        self.pc = power_condition
        self.u_mom = u_mom
        self.alpha = alpha
        self.ref_area = ref_area

    """ The inflow speed for the nacelle is determined by the propeller and
    hence the speed after the propeller is used """
    def inflow_velocity(self):
        if self.pc == "off":
            u_nacelle = flow_conditions.u_inf
            return u_nacelle
        else:
            u_nacelle = self.u_mom
            return u_nacelle

    def inflow_angle(self):
        inflow_nacelle = self.alpha
        return inflow_nacelle

    """ For the area of the nacelle, geometric properties are used to determine
     the outside area of the nacelle. """

    def wet_area(self):
        """ only the side area is assumed no closing sides"""
        area_cylinder = np.pi * self.nacelle_diameter * self.nacelle_length

        """ half a sphere is used to close of the rear end"""
        area_rear = 0.5 * np.pi * self.nacelle_diameter ** 2

        area_nacelle = area_rear + area_cylinder
        return area_nacelle

    """ The coefficients for the nacells are determined as follows """
    def cd0(self):
        cf = skin_friction(self.re_induct, "t")
        fm = mach_correction(flow_conditions.Mach)
        l_d = self.nacelle_length / self.nacelle_diameter
        f_nac = 1 + 60 / l_d ** 3 + 0.0025 * l_d
        norm_area = self.wet_area() / self.ref_area

        cd0_nacelle = cf * fm * f_nac * norm_area

        return cd0_nacelle

    def cd_prime(self):
        norm_speed = self.inflow_velocity() / flow_conditions.u_inf

        cd_nacelle = self.cd0() * norm_speed ** 2
        return cd_nacelle

    @staticmethod
    def cl_prime():
        """ assume no lift is produced by the nacelle"""
        cl_nacelle = 0
        return cl_nacelle

    """ The weight is depending on the propulsor type """

    def weight(self):
        if self.propulsor_type == 'conventional':
            weight = 100000
            return weight

        else:
            weight = 50000
            return weight
