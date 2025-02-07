import constants
import numpy as np
from analysis_modules.ISA import air_density_isa


def reynolds(isa_input, u, char_len, component):
    rho = isa_input[0]
    temperature = isa_input[1]
    if rho == 1.225:
        re = (rho * u * char_len) / constants.mu_air
        print("\n----- Reynolds number calculator:")
        print(f'Kinematic viscosity of the air =', np.round(constants.mu_air),
              2)
        print(f'Reynolds numer is of {component} =', np.round(re, 0))
        return re
    else:
        # Sutherland formula for dynamic viscosity
        mu = constants.mu_air * ((temperature / 288.15) ** 1.5) * (
                    (288.15 + 110) / (temperature + 110))
        print("\n----- Reynolds number calculator:")
        print(f'Kinematic viscosity of the air =', np.round(mu, 2))
        re = (rho * u * char_len) / mu
        print(f'Reynolds number of {component} is =', np.round(re, 0))
        return re


def speed_of_sound(altitude):
    temperature = air_density_isa(altitude)[1]

    a = np.sqrt(constants.gamma_air * constants.R * temperature)
    print("\n----- Speed of Sound calculator:")
    print(f"Speed of sound =", np.round(a, 2), "m/s")
    return a


def drag(cd, rho, speed, area, component):
    """ Calculates the 3D drag based on a 2D coefficient"""
    three_dim_drag = cd * 0.5 * rho * speed ** 2 * area

    print(f"Drag from {component} = {three_dim_drag} [N]")
    return three_dim_drag


def lift(cl, rho, speed, area, component):
    """ Calculates the 3D lift based on a 2D coefficient"""
    three_dim_lift = cl * 0.5 * rho * speed ** 2 * area

    print(f"Lift from {component} = {three_dim_lift} [N]")
    return three_dim_lift


def moment(cm, rho, speed, area, chord, component):
    """ Calculates the pitching moment based on a 2D coefficient"""
    pitching_moment = cm * rho * speed ** 2 * area * chord

    print(f"Pitching moment from {component} = {pitching_moment} [Nm]")
    return pitching_moment

