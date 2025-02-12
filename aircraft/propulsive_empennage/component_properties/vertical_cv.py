import numpy as np
from analysis_modules.aerodynamic import lift, drag, moment
import flow_conditions
from data.read_data import airfoil_polar


class VerticalControlVane:
    """ Reference for each horizontal control vane is the attachment to the nacelle
    the values calculated are for 1 control vane. """
    def __init__(self, cv_span: float, cv_chord: float, cv_profile: str):
        super().__init__()
        self.cv_span = cv_span
        self.cv_chord = cv_chord
        self.cv_profile = cv_profile

    """ The inflow velocity is affected by the propeller """
    @staticmethod
    def inflow_velocity():
        u_control_vane = flow_conditions.u_inf * 1.25
        return u_control_vane

    """ The area is based on a rectangle and is calculated for 1 control vane """
    def area(self):
        area_control_vane = self.cv_span * self.cv_chord
        return area_control_vane

    """ For the force calculations, a 2D airfoil section is analysed for the
    given angle of attack. This is then translated to a 3D force. """
    @staticmethod
    def alpha():
        angle_vcontrol_vane = flow_conditions.delta_r
        return angle_vcontrol_vane

    def coefficients(self):
        coeff = airfoil_polar(f"vcv{self.cv_profile}.txt", self.alpha())
        cl = float(coeff[0])
        cd = float(coeff[1] + coeff[2])
        cm = float(coeff[3])
        return cl, cd, cm

    def airfoil_lift(self):
        l_airfoil = lift(self.coefficients()[0], flow_conditions.rho, self.inflow_velocity(),
                         self.area())
        return l_airfoil

    def airfoil_drag(self):
        d_airfoil = drag(self.coefficients()[1], flow_conditions.rho, self.inflow_velocity(),
                         self.area())
        return d_airfoil

    def fx(self):
        """ Correct for angle of attack """
        fx_vcv = (self.airfoil_lift() * np.sin(np.radians(self.alpha()))
                  + self.airfoil_drag() * np.cos(np.radians(self.alpha())))

        if fx_vcv >= 0:
            print(f"Drag = {np.round(fx_vcv, 1)} [N]")
        else:
            print(f"Thrust = {np.round(fx_vcv, 1)} [N]")
        return fx_vcv

    def fy(self):
        """ Correct for angle of attack """
        fy_vcv = (self.airfoil_lift() * np.cos(np.radians(self.alpha())) -
                  self.airfoil_drag() * np.sin(np.radians(self.alpha())))

        print(f"Side force = {np.round(fy_vcv,1)} [N]")
        return fy_vcv

    @staticmethod
    def fz():
        """ assume control vane does not produce a lift force in normal conditions """
        fz_vcv = 0
        print(f"Lift = {np.round(fz_vcv,1)} [N]")
        return fz_vcv

    def moment(self):
        pitching_moment_vcv = moment(self.coefficients()[2], flow_conditions.rho,
                                     self.inflow_velocity(), self.area(),
                                     self.cv_chord)
        print(f"Moment = {np.round(pitching_moment_vcv, 1)} [Nm]")
        return pitching_moment_vcv
