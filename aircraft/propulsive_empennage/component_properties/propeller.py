import numpy as np
import constants
import flow_conditions
from analysis_modules.BEM import BEM


class Propeller:
    def __init__(self, n_blades: float, prop_diameter: float, hub_diameter: float,
                 prop_airfoil: str, prop_sweep: float, prop_pitch: float, rpm: float):
        super().__init__()
        self.n_blades = n_blades
        self.prop_diameter = prop_diameter
        self.hub_diameter = hub_diameter
        self.prop_airfoil = prop_airfoil
        self.prop_sweep = prop_sweep
        self.prop_pitch = prop_pitch
        self.rpm = rpm

    def thrust(self):
        prop_thrust = BEM(self.prop_pitch, flow_conditions.u_inf, self.prop_diameter,
                          self.n_blades, self.prop_airfoil, self.rpm)
        return prop_thrust

    def weight(self):
        prop_weight = 10 * self.prop_diameter
        return prop_weight


duuc: Propeller = Propeller(n_blades=3, prop_diameter=3.6,
                            hub_diameter=0.2, prop_airfoil='ARAD8',
                            prop_sweep=0, prop_pitch=2, rpm=6000)

duuc.thrust()
