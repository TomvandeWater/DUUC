import numpy as np
from analysis_modules.factors import mach_correction
from analysis_modules.factors import skin_friction
from data.read_data import airfoil_polar


class Duct:
    def __init__(self, duct_diameter: float, duct_chord: float, duct_profile: str,
                 alpha: float, re_duct: float, power_condition: str, u_mom: float,
                 tc_prop: float, v_inf: float, mach: float):
        super().__init__()
        self.duct_diameter = duct_diameter
        self.duct_chord = duct_chord
        self.duct_profile = duct_profile
        self.alpha = alpha
        self.re_duct = re_duct
        self.pc = power_condition
        self.u_mom = u_mom
        self.tc_prop = tc_prop
        self.v_inf = v_inf
        self.mach = mach

    """ Define velocities and angles"""
    def inflow_velocity(self):
        if self.pc == "off":
            u_duct = self.v_inf
            return u_duct
        else:
            u_duct = self.u_mom
            return u_duct

    def inflow_angle(self):
        inflow_duct = self.alpha
        return inflow_duct

    """ For the area calculation, the project area and wetted area are differentiated, 
    thickness ratio is defined based on NACA airfoil"""
    def wetted_area(self):
        num_list = [int(digit) for digit in self.duct_profile]
        thickness = num_list[2] * 10 + num_list[3]  # naca thickness of profile
        t_max = thickness / 100

        s_wet = (2 * np.pi * self.duct_chord * self.duct_diameter + 2 * np.pi
                 * self.duct_diameter * 0.5 * t_max * self.duct_chord)
        return s_wet

    def proj_area(self):
        proj_area = self.duct_chord * self.duct_diameter
        return proj_area

    def aspect_ratio(self):
        ar_duct = self.duct_diameter / self.duct_chord
        return ar_duct

    def t_c(self):
        num_list = [int(digit) for digit in self.duct_profile]
        thickness = num_list[2] * 10 + num_list[3]  # naca thickness of profile
        thickness = thickness / 100  # returns value in percentage of normalized chord
        return thickness

    """" coefficient """
    def zeta(self):
        delta = 1 / self.aspect_ratio()
        zeta_duct = 1 / (1 + delta * np.pi / 2 + np.arctan(1.2 * delta) * delta)
        return zeta_duct

    def cl_a(self):
        cl_a_naca = 2 * np.pi
        return cl_a_naca

    def cl_da(self):
        cl_da_duct = np.pi / 2 * self.zeta() * self.cl_a()
        return cl_da_duct

    def cl(self):
        k_prop = 0.2 * np.sqrt(self.tc_prop)
        if self.pc == "off":
            cl_duct = self.cl_da() * self.inflow_angle()
            return cl_duct
        else:
            cl_duct = (1 + k_prop) * self.cl_da() * self.inflow_angle()
            return cl_duct

    def cd0(self):
        cf = skin_friction(self.re_duct, 'T')
        fm = mach_correction(self.mach)
        ftc = 1 + 2.7 * self.t_c() + 100 * self.t_c() ** 4
        coeff = airfoil_polar(f"support{self.duct_profile}.txt", float(0.0))
        cdmin = float(coeff[1] + coeff[2])

        cd0_duct = fm * ftc * cf * self.wetted_area()/self.proj_area() * (cdmin / 0.004) ** 4
        return cd0_duct

    def cdi(self):
        cdi_duct = self.cl() ** 2 / (2 * np.pi * self.aspect_ratio())
        return cdi_duct

    def cd_prime(self):
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2

        cd_duct = (self.cd0() + self.cdi()) * norm_speed
        return cd_duct

    def cl_prime(self):
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2

        cl_duct = self.cl() * norm_speed
        return cl_duct
