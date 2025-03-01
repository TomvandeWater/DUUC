import constants
import numpy as np
from analysis_modules.ISA import air_density_isa


def reynolds(isa_input, u, char_len):
    rho = isa_input[0]
    temperature = isa_input[1]
    if rho == 1.225:
        re = (rho * u * char_len) / constants.mu_air
        # print("\n----- Reynolds number calculator -----")
        # print(f'Kinematic viscosity of the air =',
        # np.round(constants.mu_air, 5), '[m^2/s]')
        # print(f'Reynolds numer is of {component} =', np.round(re, 0),
        # '[-]')
        return re
    else:
        # Sutherland formula for dynamic viscosity
        mu = constants.mu_air * ((temperature / 288.15) ** 1.5) * (
                    (288.15 + 110) / (temperature + 110))
        # print("\n----- Reynolds number calculator -----")
        # print(f'Kinematic viscosity of the air =', np.round(mu, 5),
        # '[m^2/s]')
        re = (rho * u * char_len) / mu
        # print(f'Reynolds number of {component} is =', np.round(re, 0),
        # '[-]')
        return re


def speed_of_sound(altitude):
    temperature = air_density_isa(altitude)[1]

    a = np.sqrt(constants.gamma_air * constants.R * temperature)
    # print("\n----- Speed of Sound calculator -----")
    # print(f"Speed of sound =", np.round(a, 2), "m/s")
    return a


def drag(cd, rho, speed, area):
    """ Calculates the 3D drag based on a 2D coefficient"""
    three_dim_drag = cd * 0.5 * rho * (speed ** 2) * area

    return three_dim_drag


def lift(cl, rho, speed, area):
    """ Calculates the 3D lift based on a 2D coefficient"""
    three_dim_lift = cl * 0.5 * rho * (speed ** 2) * area

    return three_dim_lift


def moment(cm, rho, speed, area, chord):
    """ Calculates the pitching moment based on a 2D coefficient"""
    pitching_moment = cm * rho * (speed ** 2) * area * chord

    return pitching_moment


def c_root(area, span, taper_ratio):
    c_root_val = (2 * area) / (span * (1 + taper_ratio))
    return c_root_val


def c_tip(c_root_val, taper_ratio):
    c_tip_val = c_root_val * taper_ratio
    return c_tip_val


def drag_zero_lift(cf, ffc, qc, s_ratio):
    cd_0 = cf * ffc * qc * s_ratio
    return cd_0


def drag_interference(t_c, connection_type: str):
    if connection_type == 'plane':
        """ This is based on the Hoerner book for a strut on a plane wall"""
        cd_int = 0.8 * t_c ** 3 - 0.0003
        return cd_int
    if connection_type == 't-junction':
        cd_int = 17 * t_c ** 2 - 0.05
        return cd_int
    else:
        print("Connection type not defined")
        return None


def propulsive_efficiency(fx, velocity, power_input, k1, k2):
    eta_propulsive = (fx * velocity) / (power_input * k1 * k2)
    return eta_propulsive


def cl_cd_bucket_polar(cd0, cdi, cl):
    cd = cd0 + cdi * cl ** 2
    return cd


def tail_volume(surface, lever_arm, surface_wing, reference_l):
    tail_volume_coefficient = (surface * lever_arm) / (surface_wing * reference_l)
    return tail_volume_coefficient

