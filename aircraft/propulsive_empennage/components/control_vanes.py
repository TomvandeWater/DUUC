from parapy.core import *
from parapy.geom import *
from geometry_modules.airfoil_read import Airfoil
import constants
import flow_conditions


class ControlVanes(GeomBase):
    """Input section"""
    cv_airfoil = Input("Naca0006")
    cv_length = Input(3)
    cv_chord = Input(0.5)
    cv_t_c = Input(2)

    """Attributes"""
    @Attribute
    def cv_forces(self):
        cv_l = 1
        cv_d = 1
        cv_t = 1
        cv_m = 1
        return cv_l, cv_d, cv_t, cv_m

    @Attribute
    def cv_weight(self):
        cv_weight = 15 * constants.g
        return cv_weight

    @Attribute
    def cv_cog(self):
        cog_x = 1
        cog_y = 2
        cog_z = 3
        return cog_x, cog_y, cog_z

    """Parts"""
    @HiddenPart
    def hcv_r_profile(self):
        return Airfoil(airfoil_name=self.cv_airfoil,
                       chord=self.cv_chord,
                       t_c=self.cv_t_c,
                       position=translate(rotate(self.position,
                                                 'y', 180,
                                                 deg=True),
                                          'x', -0.5 * self.cv_chord,
                                          'y', 0.5 * self.cv_length))

    @HiddenPart
    def hcv_t_profile(self):
        return Airfoil(airfoil_name=self.cv_airfoil,
                       chord=self.cv_chord,
                       t_c=self.cv_t_c,
                       position=translate(rotate(self.position,
                                                 'y', 180,
                                                 deg=True),
                                          'x', -0.5 * self.cv_chord,
                                          'y', -0.5 * self.cv_length))

    @HiddenPart
    def vcv_r_profile(self):
        return Airfoil(airfoil_name=self.cv_airfoil,
                       chord=self.cv_chord,
                       t_c=self.cv_t_c,
                       position=translate(rotate(rotate(self.position,
                                                 'y', 180,
                                                        deg=True),
                                                 'x', 90,
                                                 deg=True),
                                          'x', -0.5 * self.cv_chord,
                                          'y', 0.5 * self.cv_length))

    @HiddenPart
    def vcv_t_profile(self):
        return Airfoil(airfoil_name=self.cv_airfoil,
                       chord=self.cv_chord,
                       t_c=self.cv_t_c,
                       position=translate(rotate(rotate(self.position,
                                                 'y', 180,
                                                        deg=True),
                                                 'x', 90,
                                                 deg=True),
                                          'x', -0.5 * self.cv_chord,
                                          'y', -0.5 * self.cv_length))

    @Part
    def horizontal_control_vanes(self):
        return LoftedSolid(profiles=[self.hcv_r_profile,
                                     self.hcv_t_profile],
                           color='Blue',
                           mesh_deflection=0.0001,
                           label='Horizontal Control Vanes')

    @Part
    def vertical_control_vanes(self):
        return LoftedSolid(profiles=[self.vcv_r_profile,
                                     self.vcv_t_profile],
                           color='Blue',
                           mesh_deflection=0.0001,
                           label='Vertical Control Vanes')


if __name__ == '__main__':
    from parapy.gui import display
    obj = (ControlVanes(label="Control Vanes"))
    display(obj)
