from parapy.core import *
from parapy.geom import *
import numpy as np


class CenterBody(GeomBase):
    """Input section"""
    cb_diameter = Input(0.25)  # center body diameter [m]
    cb_length = Input(2.5)  # center body length [m]

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
        cb_l = 1
        cb_d = 1
        cb_t = 1
        cb_m = 1
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
