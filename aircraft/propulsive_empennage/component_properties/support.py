import numpy as np
from analysis_modules.aerodynamic import lift, drag, moment
import flow_conditions
from data.read_data import airfoil_polar


class SupportStrut:
    def __init__(self, support_length: float, support_chord: float, support_profile: str,
                 cant_angle: float):
        self.support_length = support_length
        self.support_chord = support_chord
        self.support_profile = support_profile
        self.cant_angle = cant_angle

    """ Inflow velocity on the strut is affected by the propeller"""
    @staticmethod
    def inflow_velocity():
        u_support = flow_conditions.u_inf * 1.25
        return u_support

    """" The area is based on a rectangle, note that the support goes through the nacelle in this 
    model. The area is hence overestimated. """
    def area(self):
        s_support = self.support_chord * self.support_length
        return s_support

    @staticmethod
    def alpha():
        """ assume undisturbed flow for now"""
        alpha = 0
        return alpha

    """ Forces are determined based on a 2D airfoil analysis"""
    def coefficients(self):
        coeff = airfoil_polar(f"support{self.support_profile}.txt", float(self.alpha()))
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
        fx_pylon = (self.airfoil_lift() * np.sin(np.radians(self.alpha()))
                    + self.airfoil_drag() * np.cos(np.radians(self.alpha())))

        if fx_pylon >= 0:
            print(f"Drag = {np.round(fx_pylon,1)} [N]")
        else:
            print(f"Thrust = {np.round(fx_pylon, 1)} [N]")
        return fx_pylon

    def fy(self):
        """ Correct for angle of attack, then break down for cant angle"""
        fy = (self.airfoil_lift() * np.cos(np.radians(self.alpha())) -
              self.airfoil_drag() * np.sin(np.radians(self.alpha())))

        fy_pylon = fy * np.sin(np.radians(self.cant_angle))
        print(f"Side force = {np.round(fy_pylon,1)} [N]")
        return fy_pylon

    def fz(self):
        """ Correct for angle of attack, then break down for cant angle"""
        fz = (self.airfoil_lift() * np.cos(np.radians(self.alpha())) -
              self.airfoil_drag() * np.sin(np.radians(self.alpha())))

        fz_pylon = fz * np.cos(np.radians(self.cant_angle))
        print(f"Lift = {np.round(fz_pylon,1)} [N]")
        return fz_pylon

    def moment(self):
        """ still has to be adjusted for the translation for the pylon and reference frame"""
        pitching_moment = moment(self.coefficients()[2],
                                 flow_conditions.rho,
                                 self.inflow_velocity(), self.area(),
                                 self.support_chord)
        print(f"Moment = {np.round(pitching_moment, 1)} [Nm]")
        return pitching_moment

