from parapy.geom import *
from parapy.core import *
import os


class Airfoil(FittedCurve):
    chord = Input(1.)
    airfoil_name = Input("Naca0006")
    t_c = Input(2)
    mesh_deflection = Input(0.0001)

    @Attribute
    def points(self):
        file_path = os.path.join(os.path.dirname(__file__), '..', 'data',
                                 self.airfoil_name + ".dat")

        with open(file_path, 'r') as f:
            points = []
            for line in f:
                x, y = line.split(' ', 1)
                points.append(self.position.translate(
                   "x", float(y)*self.chord,
                   "y", float(x)*self.chord*self.t_c))
        return points


if __name__ == '__main__':
    from parapy.gui import display
    obj = Airfoil(label="Airfoil profile")
    display(obj)
