import numpy as np
from parapy.core import *
from parapy.geom import *
from geometry_modules.airfoil_read import Airfoil
from analysis_modules.aerodynamic import reynolds
from analysis_modules.ISA import air_density_isa
import constants
import flow_conditions
from data.read_data import search_dat_file


class Pylon(GeomBase):
    """Input section"""
    pylon_airfoil = Input("Naca0006")
    pylon_length = Input(5)
    pylon_chord = Input(0.5)
    pylon_t_c = Input(2)
    cant_angle = Input(30)
    duct_diameter = Input(3.6)

    """ Reference point """
    @Part
    def bbox(self):
        return Box(0.1, 0.1, 0.1, color='red')

    """Attributes"""
    @Attribute
    def pylon_conditions(self):
        # pylon outside the duct
        v_out = flow_conditions.u_inf
        pylon_re_out = reynolds(air_density_isa(flow_conditions.altitude)[0],
                                air_density_isa(flow_conditions.altitude)[1],
                                v_out, self.pylon_chord)

        # pylon inside the duct
        v_in = 2 * flow_conditions.u_inf
        pylon_re_in = reynolds(air_density_isa(flow_conditions.altitude)[0],
                               air_density_isa(flow_conditions.altitude)[1],
                               v_in, self.pylon_chord)

        return pylon_re_out, pylon_re_in, v_out, v_in

    @Attribute
    def area_pylon(self):
        """ determine area of the pylon in and outside the duct"""
        s_pylon_in = self.pylon_chord * self.duct_diameter
        s_pylon_out = ((self.pylon_length - self.duct_diameter)
                       * self.pylon_chord)

        return s_pylon_in + s_pylon_out, s_pylon_in, s_pylon_out

    @Attribute
    def pylon_forces(self):
        """ Read coefficients corresponding to alpha from polar [-5  to 15]"""
        coefficient = search_dat_file("naca0012_inv.dat",
                                      flow_conditions.alpha)

        """ Determine lift force on the pylon"""
        pylon_l_in = coefficient[0] * (0.5 * flow_conditions.rho
                                       * pow(self.pylon_conditions[3], 2)
                                       * self.area_pylon[1])
        pylon_l_out = coefficient[0] * (0.5 * flow_conditions.rho
                                        * pow(self.pylon_conditions[2], 2)
                                        * self.area_pylon[2])

        pylon_l = ((pylon_l_in + pylon_l_out)
                   * np.cos(np.radians(self.cant_angle))
                   * np.cos(np.radians(flow_conditions.alpha)))

        """ Determine drag force on the pylon"""
        pylon_d_in = coefficient[1] * (0.5 * flow_conditions.rho
                                       * pow(self.pylon_conditions[3], 2)
                                       * self.area_pylon[1])
        pylon_d_out = coefficient[1] * (0.5 * flow_conditions.rho
                                        * pow(self.pylon_conditions[2], 2)
                                        * self.area_pylon[2])

        pylon_d = ((pylon_d_in + pylon_d_out)
                   * np.sin(np.radians(flow_conditions.alpha)))

        """ Determine thrust force on the pylon"""
        pylon_t = 0  # Assume not thrust produced on the pylon

        """ Determine moment around the pylon ac"""
        pylon_m_in = coefficient[2] * (0.5 * flow_conditions.rho
                                       * pow(self.pylon_conditions[3], 2)
                                       * self.area_pylon[1]
                                       * self.pylon_chord)
        pylon_m_out = coefficient[2] * (0.5 * flow_conditions.rho
                                        * pow(self.pylon_conditions[2], 2)
                                        * self.area_pylon[2]
                                        * self.pylon_chord)

        pylon_m = pylon_m_in + pylon_m_out

        return (np.round(pylon_l, 1), np.round(pylon_d, 1),
                np.round(pylon_t, 1), np.round(pylon_m, 1))

    def drag(self):
        drag_calc = 200
        return drag_calc

    @Attribute
    def pylon_weight(self):
        """ Based on weight approximation treating it as a horizontal
        stabilizer (Vos, 2018) multiplied by correction factor found in
        research Stavreva 2020"""

        pylon_weight = ((self.area_pylon * (3.81 * pow(self.area_pylon, 0.2)
                                            * flow_conditions.V_d - 0.287))
                        * constants.g) * constants.K_pylon
        return pylon_weight

    """Parts"""
    @HiddenPart
    def pylon_r_profile(self):
        return Airfoil(airfoil_name=self.pylon_airfoil,
                       chord=self.pylon_chord,
                       t_c=self.pylon_t_c,
                       position=translate(rotate(self.position,
                                                 'y', 180,
                                                 deg=True),
                                          'x', -0.5 * self.pylon_chord))

    @HiddenPart
    def pylon_t_profile(self):
        return Airfoil(airfoil_name=self.pylon_airfoil,
                       chord=self.pylon_chord,
                       t_c=self.pylon_t_c,
                       position=translate(rotate(self.position,
                                                 'y', 180,
                                                 deg=True),
                                          'x', -0.5 * self.pylon_chord,
                                          'y', -self.pylon_length))

    @Part
    def pylon(self):
        return LoftedSolid(profiles=[self.pylon_r_profile,
                                     self.pylon_t_profile],
                           color='White',
                           mesh_deflection=0.0001)


if __name__ == '__main__':
    from parapy.gui import display
    obj = (Pylon(label="Pylon"))
    display(obj)
