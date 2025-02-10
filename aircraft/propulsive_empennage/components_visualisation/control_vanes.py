from parapy.core import *
from parapy.geom import *
from geometry_modules.airfoil_read import Airfoil
import constants
import flow_conditions
from data.read_data import search_dat_file
import numpy as np


class ControlVanes(GeomBase):
    """Input section"""
    cv_airfoil = Input("Naca0006")
    cv_length = Input(3)
    cv_chord = Input(0.5)
    cv_t_c = Input(2)

    """Attributes"""
    @Attribute
    def hcv_forces(self):
        """ Read coefficients from XFoil run based on alpha NACA0016 """
        coefficients = search_dat_file("naca0016_inv.dat",
                                       (flow_conditions.alpha
                                        + flow_conditions.delta_e))

        s_hcv = self.cv_length * self.cv_chord

        v_cv = 2 * flow_conditions.u_inf

        """ Determine lift forces"""
        hcv_l = (coefficients[0] * (0.5 * flow_conditions.rho * v_cv ** 2
                                    * s_hcv)
                 * np.cos(np.radians(flow_conditions.alpha
                                     + flow_conditions.delta_e)))

        """ Determine drag forces"""
        hcv_d = (coefficients[1] * (0.5 * flow_conditions.rho * v_cv ** 2
                                    * s_hcv)
                 * np.sin(np.radians(flow_conditions.alpha
                                     + flow_conditions.delta_e)))

        """ Determine thrust force"""
        hcv_t = 0  # assume horizontal control vanes do not produce thrust

        """ Determine pitching moment"""
        hcv_m = (coefficients[2] * (0.5 * flow_conditions.rho * v_cv ** 2
                                    * s_hcv * self.cv_chord))

        return (np.round(hcv_l, 1), np.round(hcv_d, 1),
                np.round(hcv_t, 1), np.round(hcv_m, 1))

    @Attribute
    def vcv_forces(self):
        """ Read coefficients from XFoil run based on alpha NACA0016 """
        coefficients = search_dat_file("naca0016_inv.dat",
                                       flow_conditions.delta_r)

        s_vcv = self.cv_chord * self.cv_length

        v_cv = 2 * flow_conditions.u_inf

        """ Determine lift forces """
        vcv_l = 0  # assume rudder does not produce lift

        """ Determine drag forces """
        vcv_d = (coefficients[1] * (0.5 * flow_conditions.rho * v_cv ** 2
                                    * s_vcv)
                 * np.sin(np.radians(flow_conditions.delta_r)))

        """ Determine thrust forces """
        vcv_t = 0  # assume rudder does not produce thrust

        """ Determine pitching moment"""
        vcv_m = (coefficients[2] * (0.5 * flow_conditions.rho * v_cv ** 2
                                    * s_vcv * self.cv_chord))

        return (np.round(vcv_l, 1), np.round(vcv_d, 1),
                np.round(vcv_t, 1), np.round(vcv_m, 1))

    @Attribute
    def cv_weight(self):
        hcv_weight = ((self.cv_chord * self.cv_length)
                      * constants.K_weight_h_elevator * constants.g)
        vcv_weight = ((self.cv_chord * self.cv_length)
                      * constants.K_weight_v_elevator * constants.g)
        return hcv_weight + vcv_weight

    @Attribute
    def cv_cog(self):
        """ Symmetrical airfoil so z=0, symmetrically placed with respect
        to the y-axis so hence y=0 """
        cog_x = 0.25 * self.cv_chord  # assume cog to be at 25% of chord
        cog_y = 0
        cog_z = 0
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
