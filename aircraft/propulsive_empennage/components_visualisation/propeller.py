from parapy.core import *
from parapy.geom import *
from aircraft.propulsive_empennage.components_visualisation.propeller_blade import PropellerBlade
import numpy as np
import constants


class Propeller(GeomBase):
    """Input section"""
    n_blade = Input(3)
    propeller_length = Input(3.6)

    """Attributes"""
    @Attribute
    def blade_angle(self):
        a_blade = 360 / self.n_blade
        return a_blade

    @Attribute
    def prop_forces(self):
        prop_l = 1
        prop_d = 1
        prop_t = 1
        prop_m = 1
        return prop_l, prop_d, prop_t, prop_m

    @Attribute
    def prop_weight(self):
        prop_weight = 15 * constants.g
        return prop_weight

    @Attribute
    def prop_cog(self):
        cog_x = 1
        cog_y = 2
        cog_z = 3
        return cog_x, cog_y, cog_z

    """Parts"""
    @Part
    def propeller_blades(self):
        return PropellerBlade(pass_down="propeller_length",
                              quantify=self.n_blade,
                              position=rotate(self.position, 'x',
                                              child.index * self.blade_angle,
                                              deg=True),
                              label="Blades")


if __name__ == '__main__':
    from parapy.gui import display
    obj = (Propeller(label="Propeller"))
    display(obj)
