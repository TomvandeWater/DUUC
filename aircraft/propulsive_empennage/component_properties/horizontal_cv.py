import numpy as np
from analysis_modules.aerodynamic import lift, drag, moment
import constants
import flow_conditions
from data.read_data import airfoil_polar


class HorizontalControlVane:
    """ Reference for each horizontal control vane is the attachment to the nacelle
    the values calculated are for 1 control vane. """
    def __init__(self, cv_span: float, cv_chord: float, cv_profile: str, power_condition: str,
                 duct_diameter: float, d_exit: float, u1: float):
        super().__init__()
        self.cv_span = cv_span
        self.cv_chord = cv_chord
        self.cv_profile = cv_profile
        self.pc = power_condition
        self.d_fan = duct_diameter
        self.d_exit = d_exit
        self.u1 = u1

    """ The inflow velocity is affected by the propeller """
    def inflow_velocity(self):
        area_inlet = np.pi / 4 * self.d_fan ** 2
        area_vane = np.pi / 4 * self.d_exit ** 2
        v_vane = (flow_conditions.u_inf * area_inlet) / area_vane

        if self.pc == "off":
            v_hcv = v_vane
            return v_hcv
        else:
            v_hcv = self.u1
            return v_hcv

    @staticmethod
    def inflow_angle():
        alpha_power_off = 0
        inflow_hcv = alpha_power_off + flow_conditions.delta_e
        return inflow_hcv

    """ The area is based on a rectangle and is calculated for 1 control vane """
    def area(self):
        area_control_hvane = self.cv_span * self.cv_chord
        return area_control_hvane

    def t_c(self):
        num_list = [int(digit) for digit in self.cv_profile]
        thickness = num_list[2] * 10 + num_list[3]  # NACA thickness of profile
        thickness = thickness / 100  # returns value in percentage of normalized chord
        return thickness

    """ For the force calculations, a 2D airfoil section is analysed for the
    given angle of attack. This is then translated to a 3D force. """
    def coefficients(self):
        coeff = airfoil_polar(f"hcv{self.cv_profile}.txt", self.inflow_angle())
        cl = float(coeff[0])
        cd = float(coeff[1] + coeff[2])
        cm = float(coeff[3])
        return cl, cd, cm

    def airfoil_lift(self):
        l_airfoil = lift(self.coefficients()[0], flow_conditions.rho, self.inflow_velocity(), self.area())
        return l_airfoil

    def airfoil_drag(self):
        d_airfoil = drag(self.coefficients()[1], flow_conditions.rho, self.inflow_velocity(),
                         self.area())
        return d_airfoil

    def fx(self):
        """ Correct for angle of attack """
        fx_hcv = (self.airfoil_lift() * np.sin(np.radians(self.inflow_angle()))
                  + self.airfoil_drag() * np.cos(np.radians(self.inflow_angle())))

        if fx_hcv >= 0:
            print(f"Drag = {np.round(fx_hcv, 1)} [N]")
        else:
            print(f"Thrust = {np.round(fx_hcv, 1)} [N]")
        return fx_hcv

    @staticmethod
    def fy(self):
        """ assume control vane does not produce a side force in normal conditions """
        fy_hcv = 0
        print(f"Side force = {np.round(fy_hcv,1)} [N]")
        return fy_hcv

    def fz(self):
        """ Correct for angle of attack """
        fz_hcv = (self.airfoil_lift() * np.cos(np.radians(self.inflow_angle())) -
                  self.airfoil_drag() * np.sin(np.radians(self.inflow_angle())))

        print(f"Lift = {np.round(fz_hcv,1)} [N]")
        return fz_hcv

    def moment(self):
        pitching_moment_hcv = moment(self.coefficients()[2], flow_conditions.rho,
                                     self.inflow_velocity(), self.area(),
                                     self.cv_chord)
        print(f"Moment = {np.round(pitching_moment_hcv, 1)} [Nm]")
        return pitching_moment_hcv
