from parapy.core import *
from parapy.geom import *
import numpy as np
import flow_conditions
from analysis_modules.factors import skin_friction, mach_correction


class CenterBody(GeomBase):
    """Input section"""
    cb_diameter = Input(0.25)  # center body diameter [m]
    cb_length = Input(2.5)  # center body length [m]
    duct_diameter = Input(3.6)  # duct diameter [m]
    duct_chord = Input(2.5)  # duct chord in [m]

    """Parts"""
    @HiddenPart
    def center_body_cyl(self):
        return Cylinder(radius=self.cb_diameter/2, height=self.cb_length,
                        position=rotate(self.position, 'y', 90, deg=True))

    @HiddenPart
    def rear_shape(self):
        return Sphere(radius=self.cb_diameter / 2, angle=np.pi,
                      position=rotate(self.position, 'z', 90, deg=True))

    @Part
    def center_body(self):
        return FusedSolid(self.center_body_cyl, self.rear_shape,
                          label='Center Body', color='White')

    """Attributes"""
    @Attribute
    def cb_forces(self):
        """ Lift calculation """
        cb_l = 0

        """ Drag calculation """
        fld = (1 + (60/(self.cb_length/self.cb_diameter)**3) + 0.0025
               * self.cb_length / self.cb_diameter)

        s_wet_cb = np.pi * (self.cb_length * self.cb_diameter ** 3)
        s_duct = (2 * np.pi * self.duct_diameter * self.duct_chord * np.pi
                  * 0.12 * self.duct_chord)

        cb_d0 = (skin_friction(flow_conditions.Re, "t")
                 * mach_correction(flow_conditions.Mach) * fld
                 * (s_wet_cb / s_duct))

        cb_d = cb_d0

        """ Thrust calculation  """
        cb_t = 0  # center body is assumed to not produce thrust

        """ Moment calculation"""
        cb_m = 0
        return cb_l, cb_d, cb_t, cb_m

    @Attribute
    def cb_weight(self):
        cb_weight = 15
        return cb_weight

    @Attribute
    def cb_cog(self):
        cog_x = 1
        cog_y = 2
        cog_z = 3
        return cog_x, cog_y, cog_z


if __name__ == '__main__':
    from parapy.gui import display
    obj = CenterBody(label="Center Body")
    display(obj)
