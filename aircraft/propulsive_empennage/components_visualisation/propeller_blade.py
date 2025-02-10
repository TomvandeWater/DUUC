from parapy.core import *
from parapy.geom import *
from geometry_modules.airfoil_read_blade import Airfoil
from math import radians, tan


class PropellerBlade(LoftedSolid):
    profile = Input("Naca0006")
    propeller_length = Input(1)
    t_c = Input(1)
    propeller_chord = Input(0.2)
    sweep = Input(0)
    twist = Input(20)

    @Attribute
    def profiles(self):
        return [self.root_airfoil, self.tip_airfoil]

    """             Root airfoil displacment
    displaced by half of it chord to have it positioning at the 
    center of the airfoil                                                 """
    @HiddenPart
    def root_airfoil(self):
        return Airfoil(airfoil_name=self.profile,
                       chord=self.propeller_chord,
                       t_c=self.t_c,
                       position=translate(self.position, 'y',
                                          -0.5*self.propeller_chord),
                       mesh_deflection=0.0001)

    @HiddenPart
    def tip_airfoil(self):
        return Airfoil(airfoil_name=self.profile,
                       chord=self.propeller_chord,
                       t_c=self.t_c,
                       position=translate(
                           rotate(self.position, "z",
                                  radians(self.twist)),  # apply twist angle
                           "z", self.propeller_length,
                           "y", self.propeller_length *
                           tan(radians(self.sweep))),  # apply sweep
                       mesh_deflection=0.0001)


if __name__ == '__main__':
    from parapy.gui import display
    obj = (PropellerBlade(label="Propeller blade"))
    display(obj)
