import numpy as np
from analysis_modules.aerodynamic import lift, drag, moment, drag_interference
import flow_conditions
from data.read_data import airfoil_polar
from analysis_modules.factors import k_control, skin_friction, mach_correction, oswald


class ControlVane:
    """ Reference for each control vane is the attachment to the nacelle the values calculated
    are for 1 control vane. """
    def __init__(self, cv_span: float, cv_chord: float, cv_profile: str, power_condition: str,
                 va_inlet: float, d_exit: float, u1: float, ref_area: float, deflection: float,
                 re_inflow: float):
        super().__init__()
        self.cv_span = cv_span
        self.cv_chord = cv_chord
        self.cv_profile = cv_profile
        self.pc = power_condition
        self.va_inlet = va_inlet
        self.d_exit = d_exit
        self.u1 = u1
        self.ref_area = ref_area
        self.deflection_angle = deflection
        self.ref_inflow = re_inflow

    """ The inflow velocity is affected by the propeller """
    def inflow_velocity(self):
        area_vane = np.pi / 4 * self.d_exit ** 2
        v_vane = self.va_inlet / area_vane

        if self.pc == "off":
            v_cv = v_vane
            return v_cv
        else:
            v_cv = self.u1
            return v_cv

    def inflow_angle(self):
        alpha_power_off = 0
        deflection = self.deflection_angle
        inflow = np.radians(alpha_power_off)
        return inflow

    """ The area is based on a rectangle and is calculated for 1 control vane, assumed zero degrees
     sweep of the control vane. Aspect ratio is span^2 / area """
    def area(self):
        area_control_vane = self.cv_span * self.cv_chord
        return area_control_vane

    def aspect_ratio(self):
        aspect_ratio_control = self.cv_span ** 2 / self.area()
        return aspect_ratio_control

    def t_c(self):
        num_list = [int(digit) for digit in self.cv_profile]
        thickness = num_list[2] * 10 + num_list[3]  # NACA thickness of profile
        thickness = thickness / 100  # returns value in percentage of normalized chord
        return thickness

    def wet_area(self):
        wet_vane = 2 * (1 + 0.5 * self.t_c()) * self.cv_span * self.cv_chord
        return wet_vane

    """ determine coefficients based on theory"""
    def cd_interference(self):
        # assumed constant chord and thickness for the control vanes
        norm_area = (self.t_c() * self.cv_chord ** 2) / self.ref_area
        norm_speed = self.inflow_velocity() / flow_conditions.u_inf

        cd_cv_nac = (drag_interference(self.t_c(), "plane")
                     * norm_speed ** 2 * norm_area)

        cd_cv_duct = (drag_interference(self.t_c(), "t-junction")
                      * norm_speed ** 2 * norm_area)

        cd_interference_vcv = cd_cv_nac + cd_cv_duct
        return cd_interference_vcv

    def cl_d_deflection(self):
        cl_dd_vane = ((2 * np.pi * self.aspect_ratio() * k_control(0, self.aspect_ratio())) /
                      (2 + np.sqrt((self.aspect_ratio() ** 2 * (1 - flow_conditions.Mach ** 2))
                                   / k_control(0, self.aspect_ratio()) + 4)))
        return cl_dd_vane

    def cl(self):
        cl_vane = self.cl_d_deflection() * self.deflection_angle
        return cl_vane

    def cd0(self):
        cf = skin_friction(self.ref_inflow, "t")
        fm = mach_correction(flow_conditions.Mach)
        norm_area = self.wet_area() / self.ref_area
        ftc = 1 + 2.7 * self.t_c() + 100 * self.t_c() ** 4

        coeff = airfoil_polar(f"support{self.cv_profile}.txt", float(0.0))
        cdmin = float(coeff[1] + coeff[2])
        cd0_vane = cf * fm * ftc * norm_area * (cdmin / 0.004) ** 0.4
        return cd0_vane

    def cdi(self):
        e = oswald(self.aspect_ratio(), 0)
        cdi_vane = self.cl() ** 2 / (np.pi * self.aspect_ratio() * e)
        return cdi_vane

    def cd(self):
        cd_vane = self.cdi() + self.cd0()
        return cd_vane

    def cl_prime(self):
        norm_speed = self.inflow_velocity() / flow_conditions.u_inf
        norm_area = self.area() / self.ref_area

        alpha = self.inflow_angle()

        cl_cl = self.cl() * np.cos(alpha) * norm_speed ** 2 * norm_area
        cl_cd = self.cd() * np.sin(alpha) * norm_speed ** 2 * norm_area

        cl_vane = cl_cl + cl_cd

        return cl_vane

    def cd_prime(self):
        norm_speed = self.inflow_velocity() / flow_conditions.u_inf
        norm_area = self.area() / self.ref_area

        alpha = self.inflow_angle()

        cd_cl = self.cl() * np.sin(alpha) * norm_speed ** 2 * norm_area
        cd_cd = self.cd() * np.cos(alpha) * norm_speed ** 2 * norm_area

        cd_vane = cd_cl + cd_cd

        return cd_vane

    """ For the force calculations, a 2D airfoil section is analysed for the
    given angle of attack. This is then translated to a 3D force. """