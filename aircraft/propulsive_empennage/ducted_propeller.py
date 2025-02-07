from parapy.core import *
from parapy.geom import *
from aircraft.propulsive_empennage.components.pylon import Pylon
from aircraft.propulsive_empennage.components.duct import Duct
from aircraft.propulsive_empennage.components.control_vanes import ControlVanes
from aircraft.propulsive_empennage.components.propeller import Propeller
from aircraft.propulsive_empennage.components.center_body import CenterBody
import numpy as np


class DuctedPropeller(GeomBase):
    """Input section"""
    cant_angle = Input(30)  # cant angle pylon [deg]

    # these should be passed down from top level
    pylon_length = Input(5)
    pylon_chord = Input(1)
    duct_diameter = Input(3.6)
    duct_chord = Input(2.5)
    n_blade = Input(3)
    hub_diameter = Input(0.75)
    propeller_length = Input(1.8)

    """Attributes for positioning"""
    @Attribute
    def center_duct(self):
        y_center = (self.pylon_length * np.cos(np.radians(self.cant_angle)) -
                    (self.duct_diameter/2) *
                    np.cos(np.radians(self.cant_angle)))

        z_center = (self.pylon_length * np.sin(np.radians(self.cant_angle)) -
                    self.duct_diameter/2 * np.sin(np.radians(self.cant_angle)))
        return y_center, z_center

    """Parts"""
    @Part
    def pylon(self):
        return Pylon(pass_down="pylon_length, pylon_chord, cant_angle, "
                               "duct_diameter",
                     position=rotate(self.position, 'x',
                                     -self.cant_angle,
                                     deg=True), label='Pylon')

    @Part
    def control_vanes(self):
        return ControlVanes(cv_length=self.duct_diameter,
                            position=translate(self.position,
                                               'y', -self.center_duct[0],
                                               'z', self.center_duct[1],
                                               'x', -self.pylon_chord),
                            label='Control Vanes')

    @Part
    def propeller(self):
        return Propeller(pass_down="n_blade, propeller_length",
                         position=translate(self.control_vanes.position,
                                            'x', 2),
                         label='Propeller',
                         color=(128, 128, 128))

    @Part
    def hub_shape(self):
        return Sphere(radius=self.hub_diameter / 2, angle=np.pi,
                      position=rotate(self.propeller.position,
                                      'z', -90, deg=True),
                      label='Spinner Hub',
                      color='White')

    @Part
    def center_body(self):
        return CenterBody(pass_down='duct_diameter, duct_chord',
                          cb_diameter=self.hub_diameter, cb_length=2.5,
                          position=translate(self.propeller.position,
                                             'x', -2.5),
                          label='Center Body')

    @Part
    def duct(self):
        return Duct(pass_down="duct_diameter, duct_chord",
                    y_duct=-self.center_duct[0],
                    z_duct=self.center_duct[1]+self.duct_diameter/2,
                    label="Duct",
                    position=translate(self.position, 'x', 1.5))

    """Attributes"""
    @Attribute
    def ducted_propeller_forces(self):
        ductp_l = 1
        ductp_d = 2
        ductp_t = 3
        ductp_m = 4
        return ductp_l, ductp_d, ductp_t, ductp_m

    @Attribute
    def ducted_propeller_cog(self):
        x_cog = 1
        y_cog = 2
        z_cog = 3
        return x_cog, y_cog, z_cog

    @Attribute
    def ducted_propeller_weight(self):
        ductp_weight = 10
        return ductp_weight


if __name__ == '__main__':
    from parapy.gui import display
    obj = (DuctedPropeller(label="Ducted Propeller"))
    display(obj)
