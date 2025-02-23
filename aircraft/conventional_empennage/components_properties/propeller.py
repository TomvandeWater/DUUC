import numpy as np
from analysis_modules.BEM_simplified import BEM
from analysis_modules.factors import skin_friction
import flow_conditions


class Propeller:
    """ propeller class for the conventional configuration"""
    def __init__(self, rpm: float, alpha: float, power_condition: str, n_blades: float,
                 prop_diameter: float, hub_diameter: float, prop_airfoil: str, prop_sweep: float,
                 prop_pitch: float, c_root: float, c_tip: float, v_inf: float,
                 area_ref: float, reynolds: float):
        super().__init__()
        self.rpm = rpm
        self.alpha = alpha
        self.pc = power_condition
        self.n_blades = n_blades
        self.prop_diameter = prop_diameter
        self.hub_diameter = hub_diameter
        self.prop_airfoil = prop_airfoil
        self.prop_sweep = prop_sweep
        self.prop_pitch = prop_pitch
        self.c_root = c_root
        self.c_tip = c_tip
        self.v_inf = v_inf
        self.area_ref = area_ref
        self.reynolds = reynolds

    """ inflow velocity and angles """
    def inflow_velocity(self):
        v_prop = self.v_inf
        return v_prop

    def inflow_angle(self):
        a_prop = self.alpha
        return a_prop

    """ areas and geometric properties"""
    def area(self):
        area_propel = np.pi * (self.prop_diameter / 2) ** 2
        return area_propel

    @staticmethod
    def t_c():
        """ defined on comparable propeller blade"""
        t_c_prop_root = 0.10
        t_c_prop_tip = 0.08
        t_c_av = (t_c_prop_tip + t_c_prop_root) / 2
        return t_c_prop_root, t_c_prop_tip, t_c_av

    def wet_area(self):
        wet_prop_bl = (2 * (1 + 0.5 * self.t_c()[2]) * self.prop_diameter / 2
                       * (self.c_root + self.c_tip) / 2)

        wet_prop = self.n_blades * wet_prop_bl
        return wet_prop

    """ define coefficients """
    def cd0(self):
        if self.pc == "off":
            """ in power off conditions the drag is estimated by only the skin friction drag"""
            cd0_prop = skin_friction(self.reynolds, "t") * self.n_blades
            return cd0_prop
        else:
            """ in power on conditions the skin friction drag is assumed to be neglible"""
            cd0_prop = 0
            return cd0_prop

    def cd_prime(self):
        a = np.radians(self.alpha)
        inflow_a = np.radians(self.inflow_angle())
        norm_area = self.wet_area() / self.area_ref
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2

        if self.pc == "off":
            cd_prime = (self.cd0() * np.cos(a - inflow_a) * self.n_blades * norm_area
                        * norm_speed)
            return cd_prime
        else:

            cd_prime = self.cn() * np.sin(a - inflow_a) * norm_area * norm_speed
            return cd_prime

    def cl_prime(self):
        a = np.radians(self.alpha)
        inflow_a = np.radians(self.inflow_angle())
        norm_area = self.wet_area() / self.area_ref
        norm_speed = self.inflow_velocity() ** 2 / self.v_inf ** 2

        if self.pc == "off":
            cl_prime = (self.cd0() * np.sin(a - inflow_a) * self.n_blades * norm_area
                        * norm_speed)
            return cl_prime
        else:
            cl_1 = self.tc() * np.sin(a)

            cl_2 = self.cn() * np.cos(a - inflow_a) * norm_area * norm_speed

            cl_prime = cl_1 + cl_2
            return cl_prime

    def thrust(self):
        """
        prop = BEM(radius=self.prop_diameter / 2,
                   num_blades=self.n_blades,
                   rpm=self.rpm,
                   prop_airfoil=self.prop_airfoil)
        thrust_prop, torque_prop, normal_f = prop.calculate_thrust(self.v_inf)"""
        thrust_prop = 8000
        return thrust_prop

    def tc(self):
        tc_prop = self.thrust() / (0.5 * flow_conditions.rho * self.inflow_velocity() ** 2 * self.area())
        return tc_prop

    def cn(self):
        """
        prop = BEM(radius=self.prop_diameter / 2,
                   num_blades=self.n_blades,
                   rpm=self.rpm,
                   prop_airfoil=self.prop_airfoil)
        thrust_prop, torque_prop, normal_f = prop.calculate_thrust(self.v_inf) """
        normal_f = 3200
        cn_prop = normal_f / (0.5 * flow_conditions.rho * self.inflow_velocity() ** 2
                              * self.area())
        return cn_prop

    def u1(self):
        u1_prop = np.sqrt((2 / flow_conditions.rho) * (self.thrust() / self.area()) +
                          self.inflow_velocity() ** 2)
        return u1_prop
