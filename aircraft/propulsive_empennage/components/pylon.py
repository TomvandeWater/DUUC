from parapy.core import *
from parapy.geom import *
from geometry_modules.airfoil_read import Airfoil
from analysis_modules.aerodynamic import reynolds
from analysis_modules.ISA import air_density_isa
import constants
import flow_conditions


class Pylon(GeomBase):
    """Input section"""
    pylon_airfoil = Input("Naca0006")
    pylon_length = Input(3)
    pylon_chord = Input(0.5)
    pylon_t_c = Input(2)

    """ reference point """
    @Part
    def bbox(self):
        return Box(0.1, 0.1, 0.1, color='red')

    """Attributes"""
    @Attribute
    def pylon_conditions(self):
        # pylon reynolds number outside the duct
        pylon_re_out = reynolds(air_density_isa(constants.altitude),
                                constants.u_inf,
                                self.pylon_chord)

        # pylon reynolds number inside the duct
        u_prop = 2 * constants.u_inf
        pylon_re_in = reynolds(air_density_isa(constants.altitude),
                               u_prop,
                               self.pylon_chord)

        return pylon_re_out, pylon_re_in

    @Attribute
    def pylon_forces(self):
        pylon_l = 1
        pylon_d = 1
        pylon_t = 1
        pylon_m = 1
        return pylon_l, pylon_d, pylon_t, pylon_m

    @Attribute
    def area_pylon(self):
        s_pylon = self.pylon_length * self.pylon_chord
        return s_pylon

    @Attribute
    def pylon_weight(self):
        """ Based on weight approximation treating it as a horizontal
        stabilizer (Vos, 2018) """
        pylon_weight = ((self.area_pylon * (3.81 * self.area_pylon ^ 0.2
                                            * flow_conditions.V_d - 0.287))
                        * constants.g)
        return pylon_weight

    @Attribute
    def pylon_cog(self):
        cog_x = 1
        cog_y = 2
        cog_z = 3
        return cog_x, cog_y, cog_z

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
