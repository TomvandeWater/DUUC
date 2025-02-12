import numpy as np
import constants
import flow_conditions


class Duct:
    def __init__(self, duct_diameter: float, duct_chord: float, duct_profile: str,
                 alpha: float):
        super().__init__()
        self.duct_diameter = duct_diameter
        self.duct_chord = duct_chord
        self.duct_profile = duct_profile
        self.alpha = alpha

    @staticmethod
    def inflow_velocity():
        u_duct = flow_conditions.u_inf
        return u_duct

    """ For the area calculation, the project area and wetted area are differentiated"""
    def wetted_area(self):
        s_wet = np.pi * self.proj_area()
        return s_wet

    def proj_area(self):
        proj_area = self.duct_chord * self.duct_diameter
        return proj_area

    """" Force calculation """
