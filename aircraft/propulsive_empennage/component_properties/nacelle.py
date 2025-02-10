import numpy as np
import constants
import flow_conditions
from analysis_modules.factors import skin_friction, mach_correction

class Nacelle:

    def __init__(self, nacelle_length: float, nacelle_diameter: float,
                 propulsor_type: str, re_induct: float):
        super().__init__()
        self.nacelle_length = nacelle_length
        self.nacelle_diameter = nacelle_diameter
        self.propulsor_type = propulsor_type
        self.re_induct = re_induct

    """ The inflow speed for the nacelle is determined by the propeller and
    hence the speed after the propeller is used """

    @staticmethod
    def velocity():
        u_nacelle = flow_conditions.u_inf * 1.25
        return u_nacelle

    """ For the area of the nacelle, geometric properties are used to determine
     the outside area of the nacelle. """

    def area(self):
        """ only the side area is assumed no closing sides"""
        area_cylinder = np.pi * self.nacelle_diameter * self.nacelle_length

        """ half a sphere is used to close of the rear end"""
        area_rear = 0.5 * np.pi * self.nacelle_diameter ** 2

        area_nacelle = area_rear + area_cylinder
        return area_nacelle

    """ The weight is depending on the propulsor type """

    def weight(self):
        if self.propulsor_type == 'conventional':
            weight = 100000
            return weight

        else:
            weight = 50000
            return weight

    """ The nacelle forces are determined """

    def fx(self):
        cf = skin_friction(self.re_induct, "t")
        fm = mach_correction(flow_conditions.Mach)
        l_d = self.nacelle_length / self.nacelle_diameter
        area_ratio = self.area() / (0.25 * np.pi * self.nacelle_diameter ** 2)

        cd_nacelle = (cf * fm * area_ratio
                      * (1 + (60/(l_d ** 3)) + 0.0225 * l_d))

        s_flat = self.nacelle_diameter * self.nacelle_length

        fx_nacelle = cd_nacelle * (0.5 * flow_conditions.rho * self.velocity()
                                   * s_flat)
        return fx_nacelle

    @staticmethod
    def fy(self):
        """ it is assumed that the nacelle does not produce a force in y or
        z direction due to axisymmetric flow conditions"""

        fy_nacelle = 0
        return fy_nacelle

    @staticmethod
    def fz(self):
        """ it is assumed that the nacelle does not produce a force in y or
        z direction due to axisymmetric flow conditions"""

        fz_nacelle = 0
        return fz_nacelle

    def moment(self):
        M = 1000 * self.nacelle_length

        return M
