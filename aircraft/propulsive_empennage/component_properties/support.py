import numpy as np
from analysis_modules.aerodynamic import drag_interference
from data.read_data import airfoil_polar


class SupportStrut:
    def __init__(self, support_length: float, support_chord: float, support_profile: str,
                 cant_angle: float, power_condition: str, v_prop: float, u_mom: float,
                 alpha: float, tc_prop: float, cn_prop: float, ref_area: float, v_inf: float):
        super().__init__()
        self.support_length = support_length
        self.support_chord = support_chord
        self.support_profile = support_profile
        self.cant_angle = cant_angle
        self.pc = power_condition
        self.v_prop = v_prop
        self.u_mom = u_mom
        self.alpha = alpha
        self.tc_prop = tc_prop
        self.cn_prop = cn_prop
        self.ref_area = ref_area
        self.v_inf = v_inf

    """ Inflow velocity on the strut is affected by the propeller"""
    def inflow_velocity(self):
        if self.pc == "off":
            u_support = self.v_prop
            return u_support
        else:
            u_support = self.u_mom
            return u_support

    def inflow_angle(self):
        thrust_c = self.tc_prop
        c_norm = self.cn_prop
        de_da_prop = ((thrust_c / (4 + 8/7 * thrust_c)) + (c_norm * (1 + 1.3 * thrust_c) ** 0.5)
                      / (4 + 2 * thrust_c))

        if self.pc == "off":
            inflow_support = self.alpha
            return inflow_support
        else:
            inflow_support = self.alpha * (1 - de_da_prop)
            return inflow_support

    """" The area is based on a rectangle, note that the support goes through the nacelle in this 
    model. The area is hence overestimated. """
    def area(self):
        s_support = self.support_chord * self.support_length
        return s_support

    def t_c(self):
        num_list = [int(digit) for digit in self.support_profile]
        thickness = num_list[2] * 10 + num_list[3]  # naca thickness of profile
        thickness = thickness / 100  # returns value in percentage of normalized chord
        return thickness

    """ coefficients for support forces """
    def cl_da(self):
        cl0 = airfoil_polar(f"support{self.support_profile}.txt", float(0.0))
        cl0_val = float(cl0[0])
        cl5 = airfoil_polar(f"support{self.support_profile}.txt", float(10.0))
        cl5_val = float(cl5[0])
        cl_da_support = cl5_val - cl0_val / 5

        cl_da_support = np.pi * 2
        return cl_da_support

    def cl(self):
        cl_support = self.cl_da() * self.inflow_angle() * np.cos(self.cant_angle)
        return cl_support

    def cd(self):
        cd = airfoil_polar(f"support{self.support_profile}.txt", float(self.inflow_angle()))
        cd_val = float(cd[1] + cd[2])
        return cd_val

    def cd_interference(self):
        norm_area = (self.t_c() * self.support_chord ** 2) / self.ref_area
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2

        cd_support_nacelle = (2 * drag_interference(self.t_c(), "plane")
                              * norm_area * norm_speed)  # multiplied by 2 for 2 nac-supp inter

        cd_support_duct = (2 * drag_interference(self.t_c(), "t-junction")
                           * norm_area * norm_speed)  # multiplied by 2 for 2 supp-duct inter

        cd_int_support = cd_support_nacelle + cd_support_duct
        return cd_int_support

    def cl_prime(self):
        norm_area = self.area() / self.ref_area
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2

        cl_cl = (self.cl() * np.cos(np.radians(self.alpha - self.inflow_angle())) * norm_speed
                 * norm_area)

        cl_cd = (self.cd() * np.sin(np.radians(self.alpha - self.inflow_angle())) * norm_speed
                 * norm_area)

        cl_support = cl_cl + cl_cd
        return cl_support

    def cd_prime(self):
        norm_area = self.area() / self.ref_area
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2

        cd_cd = (self.cd() * np.cos(np.radians(self.alpha - self.inflow_angle())) * norm_speed
                 * norm_area)

        cd_cl = (self.cl() * np.sin(np.radians(self.alpha - self.inflow_angle())) * norm_speed
                 * norm_area)

        cd_support = cd_cd + cd_cl
        return cd_support


""" Test section """
"""
support:  SupportStrut = SupportStrut(support_length=1, support_chord=0.5, support_profile="0012",
                                      cant_angle=30, power_condition="off", v_prop=150, u_mom=50,
                                      alpha=0, tc_prop=0.8, cn_prop=0.12, ref_area=15)
"""
