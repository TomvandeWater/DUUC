import numpy as np
from analysis_modules.aerodynamic import lift, drag, moment, drag_interference
import constants
import flow_conditions
from data.read_data import airfoil_polar
import data.atr_reference as ref


class Pylon:
    """ The reference point for connecting the pylon to the fuselage is on
    y=0 and in at half of the chord length of the airfoil. The pylon is
    completely outside the duct. """
    def __init__(self, pylon_length: float, pylon_chord: float, pylon_profile: str,
                 power_condition: str, cant_angle: float, alpha: float, ref_area: float,
                 v_inf: float):
        super().__init__()
        self.pylon_length = pylon_length
        self.pylon_chord = pylon_chord
        self.pylon_airfoil = pylon_profile
        self.cant_angle = cant_angle
        self.alpha = alpha
        self.pc = power_condition
        self.ref_area = ref_area
        self.v_inf = v_inf

    """ The inflow speed for the pylon is affected by the outside of the duct
    and also affected by the downwash of the wing. The downwash of the wing and the duct affect
    the inlfow angle of on the pylon"""
    def inflow_velocity(self):
        e = (2 * ref.cl_wing) / (np.pi * ref.ar_w)

        if self.pc == "off":
            u_pylon = self.v_inf * (1 - (e ** 2) / 2)
            return u_pylon
        else:
            u_pylon = self.v_inf * (1 - (e ** 2) / 2)
            return u_pylon

    def inflow_angle(self):
        e = (2 * ref.cl_wing) / (np.pi * ref.ar_w)

        inflow_angle = self.alpha - e
        return inflow_angle

    """ For the area of the pylon, it is assumed to be a rectangle and the 
    exact curvature at the connection of the pylon is not incorporated. """

    def area(self):
        area_pylon = self.pylon_length * self.pylon_chord
        return area_pylon

    def t_c(self):
        num_list = [int(digit) for digit in self.pylon_airfoil]
        thickness = num_list[2] * 10 + num_list[3]  # naca thickness of profile
        thickness = thickness / 100  # returns value in percentage of normalized chord
        return thickness

    """ Coefficients for force calculations """

    def cl_da(self):
        cl0 = airfoil_polar(f"pylon{self.pylon_airfoil}.txt", float(0.0))
        cl0_val = float(cl0[0])
        cl5 = airfoil_polar(f"pylon{self.pylon_airfoil}.txt", float(10.0))
        cl5_val = float(cl5[0])
        cl_da_pylon = cl5_val - cl0_val / 5

        cl_da_pylon = np.pi * 2
        return cl_da_pylon

    def cl(self):
        cl_pylon = self.cl_da() * self.inflow_angle() * np.cos(np.radians(self.cant_angle))
        return cl_pylon

    def cd(self):
        cd = airfoil_polar(f"support{self.pylon_airfoil}.txt", float(self.inflow_angle()))
        cd_val = float(cd[1] + cd[2])
        return cd_val

    def cd_interference(self):
        norm_area = (self.t_c() * self.pylon_chord ** 2) / self.ref_area
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2

        cd_pylon_duct = (drag_interference(self.t_c(), "t-junction")
                         * norm_area * norm_speed)  # interference pylon duct

        cd_pylon_fuse = (drag_interference(self.t_c(), "plane")
                         * norm_area * norm_speed)  # interference pylon fuselage

        cd_int_pylon = cd_pylon_duct + cd_pylon_fuse
        return cd_int_pylon

    def cl_prime(self):
        norm_area = self.area() / self.ref_area
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2

        cl_cl = (self.cl() * np.cos(np.radians(self.alpha - self.inflow_angle())) * norm_speed
                 * norm_area)

        cl_cd = (self.cd() * np.sin(np.radians(self.alpha - self.inflow_angle())) * norm_speed
                 * norm_area)

        cl_pylon = cl_cl + cl_cd
        return cl_pylon

    def cd_prime(self):
        norm_area = self.area() / self.ref_area
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2

        cd_cd = (self.cd() * np.cos(np.radians(self.alpha - self.inflow_angle())) * norm_speed
                 * norm_area)

        cd_cl = (self.cl() * np.sin(np.radians(self.alpha - self.inflow_angle())) * norm_speed
                 * norm_area)

        cd_pylon = cd_cd + cd_cl
        return cd_pylon

    """ Weight definition of the pylon"""
    def weight(self):
        """ Based on weight approximation treating it as a horizontal
        stabilizer (Vos, 2018) multiplied by correction factor found in
        research Stavreva 2020"""

        pylon_weight = ((self.area() * (3.81 * self.area() ** 0.2
                                        * flow_conditions.V_d - 0.287))
                        * constants.g) * constants.K_pylon
        return pylon_weight

    def cog(self):
        cog_x = 0.5 * self.pylon_chord
        cog_y = - 0.5 * self.pylon_length * np.cos(np.radians(self.cant_angle))
        cog_z = 0.5 * self.pylon_length * np.sin(np.radians(self.cant_angle))
        return cog_x, cog_y, cog_z


""" test section"""
"""
pylon: Pylon = Pylon(pylon_length=config.pylon_length, pylon_chord=config.pylon_chord, pylon_profile=config.pylon_airfoil,
                     power_condition="off", cant_angle=30, alpha=0, ref_area=15,
                     v_inf=183)

pylon.cd_prime() """