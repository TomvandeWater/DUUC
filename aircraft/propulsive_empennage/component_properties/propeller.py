import numpy as np
import constants
import flow_conditions
from analysis_modules.BEM import BEM
from analysis_modules.aerodynamic import drag_interference
from analysis_modules.factors import skin_friction


class Propeller:
    def __init__(self, n_blades: float, prop_diameter: float, hub_diameter: float,
                 prop_airfoil: str, prop_sweep: float, prop_pitch: float, rpm: float,
                 power_condition: str, va_inlet: float, alpha: float, re_inflow: float,
                 s_ref: float, c_root: float, c_tip: float):
        super().__init__()
        self.n_blades = n_blades
        self.prop_diameter = prop_diameter
        self.hub_diameter = hub_diameter
        self.prop_airfoil = prop_airfoil
        self.prop_sweep = prop_sweep
        self.prop_pitch = prop_pitch
        self.rpm = rpm
        self.pc = power_condition
        self.va_inlet = va_inlet
        self.alpha = alpha
        self.re_inflow = re_inflow
        self.area_duct = s_ref
        self.c_root = c_root
        self.c_tip = c_tip

    def inflow_velocity(self):
        v_ind = self.va_inlet / (np.pi / 4 * self.prop_diameter ** 2)
        u_prop = (flow_conditions.u_inf + v_ind) / 2
        return u_prop

    def inflow_angle(self):
        inflow_prop = self.alpha
        return inflow_prop

    def area(self):
        area_prop = np.pi / 4 * self.prop_diameter ** 2
        return area_prop

    def t_c(self):
        """ define this based on airfoil profile"""
        t_c_prop_root = 0.12
        t_c_prop_tip = 0.12
        return t_c_prop_root, t_c_prop_tip

    """ Coefficient calculations """
    def cd0(self):
        if self.pc == "off":
            """ in power off conditions the drag is estimated by only the skin friction drag"""
            cd0_prop = skin_friction(self.re_inflow, "T")
            return cd0_prop
        else:

            return cd0_prop

    def cd(self):
        if self.pc == "off":
            cd = self.cd0() * self.n_blades

    def cd_interference(self):
        if self.pc == "off":
            norm_area_nac = ((self.t_c()[0] * self.c_root) * self.c_root) / self.area_duct
            norm_area_tip = ((self.t_c()[1] * self.c_tip) * self.c_tip) / self.area_duct
            norm_speed = self.inflow_velocity() / flow_conditions.u_inf

            cd_int_nac = (self.n_blades * drag_interference(self.t_c()[0], "plane")
                          * norm_speed ** 2 * norm_area_nac)

            cd_int_tip = (self.n_blades * drag_interference(self.t_c()[1], "plane")
                          * norm_speed ** 2 * norm_area_tip)

            cd_int_prop = cd_int_nac + cd_int_tip
            return cd_int_prop
        else:
            """ Assume that interference drag is a neglible contribution when power is on"""
            cd_int_prop = 0
            return cd_int_prop

    def u1(self):
        u1_prop = np.sqrt((2/flow_conditions.rho) * (T_prop / self.area()) +
                          self.inflow_velocity() ** 2)
        return u1_prop

    def thrust(self):
        prop_thrust = BEM(self.prop_pitch, flow_conditions.u_inf, self.prop_diameter,
                          self.n_blades, self.prop_airfoil, self.rpm)
        return prop_thrust

    """ Weight definition"""
    def weight(self):
        prop_weight = 10 * self.prop_diameter
        return prop_weight






duuc: Propeller = Propeller(n_blades=3, prop_diameter=3.6,
                            hub_diameter=0.2, prop_airfoil='ARAD8',
                            prop_sweep=0, prop_pitch=2, rpm=6000)

duuc.thrust()
