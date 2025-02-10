import numpy as np
from analysis_modules.aerodynamic import lift, drag, moment
import constants
import flow_conditions
from data.read_data import airfoil_polar


class Pylon:
    """ The reference point for connecting the pylon to the fuselage is on
    y=0 and in at half of the chord length of the airfoil. The pylon is
    completely outside the duct. """
    def __init__(self, pylon_length: float, pylon_chord: float,
                 pylon_profile: str,
                 cant_angle: float, alpha: float):
        super().__init__()
        self.pylon_length = pylon_length
        self.pylon_chord = pylon_chord
        self.pylon_airfoil = pylon_profile
        self.cant_angle = cant_angle
        self.alpha = alpha

    """ The inflow speed for the pylon is affected by the outside of the duct
    and potentially also affected by the downwash of the wing."""

    @staticmethod
    def inflow_velocity():
        u_pylon = flow_conditions.u_inf * 1.25
        return u_pylon

    """ For the area of the pylon, it is assumed to be a rectangle and the 
    exact curvature at the connection of the pylon is not incorporated. """

    def area(self):
        area_pylon = self.pylon_length * self.pylon_chord
        return area_pylon

    """ For the force calculations, a 2D airfoil section is analysed for the
    given angle of attack. This is then translated to a 3D force. """

    def coefficients(self):
        coeff = airfoil_polar("pylon0012.txt", self.alpha)
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
        fx_pylon = (self.airfoil_lift() * np.sin(np.radians(self.alpha))
                    + self.airfoil_drag() * np.cos(np.radians(self.alpha)))

        if fx_pylon >= 0:
            print(f"Drag = {np.round(fx_pylon,1)} [N]")
        else:
            print(f"Thrust = {np.round(fx_pylon, 1)} [N]")
        return fx_pylon

    def fy(self):
        """ Correct for angle of attack, then break down for cant angle"""
        fz = (self.airfoil_lift() * np.cos(np.radians(self.alpha)) -
              self.airfoil_drag() * np.sin(np.radians(self.alpha)))

        fy_pylon = fz * np.sin(np.radians(self.cant_angle))
        print(f"Side force = {np.round(fy_pylon,1)} [N]")
        return fy_pylon

    def fz(self):
        """ Correct for angle of attack, then break down for cant angle"""
        fz = (self.airfoil_lift() * np.cos(np.radians(self.alpha)) -
              self.airfoil_drag() * np.sin(np.radians(self.alpha)))

        fz_pylon = fz * np.cos(np.radians(self.cant_angle))
        print(f"Lift = {np.round(fz_pylon,1)} [N]")
        return fz_pylon

    def moment(self):
        pitching_moment = moment(self.coefficients()[2],
                                 flow_conditions.rho,
                                 self.inflow_velocity(), self.area(),
                                 self.pylon_chord)
        print(f"Moment = {np.round(pitching_moment, 1)} [Nm]")
        return pitching_moment

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
