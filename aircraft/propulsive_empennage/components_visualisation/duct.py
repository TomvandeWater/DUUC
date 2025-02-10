from parapy.core import *
from parapy.geom import *
from geometry_modules.airfoil_read import Airfoil
import constants


class Duct(GeomBase):
    """Input section"""
    duct_diameter = Input(3.6)  # [m]
    duct_airfoil = Input("Naca0006")
    duct_chord = Input(0.5)
    duct_t_c = Input(1)
    y_duct = Input(-10)
    z_duct = Input(10)

    """Attributes"""
    @Attribute
    def duct_forces(self):
        duct_l = 1
        duct_d = 1
        duct_t = 1
        duct_m = 1
        return duct_l, duct_d, duct_t, duct_m

    @Attribute
    def duct_weight(self):
        duct_weight = 15 * constants.g
        return duct_weight

    @Attribute
    def duct_cog(self):
        cog_x = 1
        cog_y = 2
        cog_z = 3
        return cog_x, cog_y, cog_z

    """Parts"""
    @HiddenPart
    def duct_profile(self):
        return Airfoil(airfoil_name=self.duct_airfoil,
                       chord=self.duct_chord,
                       t_c=self.duct_t_c,
                       position=translate(rotate(self.position, 'y',
                                                 180, deg=True),
                                          'y', self.y_duct,
                                          'z', -self.z_duct),
                       color='red')

    @HiddenPart
    def diameter_circle(self):
        return Circle(radius=self.duct_diameter/2,
                      position=translate(rotate(self.duct_profile.position,
                                                'y', 90, deg=True),
                                         'x', -self.duct_diameter/2))

    @Part
    def duct(self):
        return RevolvedSolid(built_from=self.duct_profile,
                             center=Point(0, self.y_duct, self.z_duct - self.duct_diameter/2),
                             direction=Vector(1, 0, 0),
                             color='White')


if __name__ == '__main__':
    from parapy.gui import display
    obj = (Duct(label="Duct"))
    display(obj)
