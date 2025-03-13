import numpy as np
import flow_conditions
from analysis_modules.aerodynamic import reynolds
from analysis_modules.ISA import air_density_isa


def s_control(aircraft_type, sweep25, l_v, power, eta, approach_velocity, n_engines, y_engine):
    if aircraft_type == 'conventional':
        delta_f_max = 25

        velocity_mc = (1.3 / 1.2) * approach_velocity
        swp = np.radians(sweep25)

        k_theta = (1 - 0.08 * np.cos(swp) ** 2) * np.cos(swp) ** (3/4)
        cl_the = 4.52  # from textNita -> still make this an interpollation based on input
        k_prime = 0.675
        cld_ratio = 0.85

        ne = ((eta * power) / (velocity_mc * n_engines)) * y_engine
        nd = 0.5 * ne

        delta_f = np.radians(delta_f_max)

        s_vertical = ((ne + nd) * 2 * np.pi) / (0.5 * flow_conditions.rho * velocity_mc ** 2 * delta_f * cld_ratio
                                                * cl_the * k_prime * k_theta * l_v)

        return s_vertical
    if aircraft_type == 'DUUC':


        return s_vertical


def s_stability(aircraft_type, s_wing, x_cog, l_f, d_f, b_w, l_v, v_crit, aspect_v, sweep50, mach):
    if aircraft_type == 'conventional':
        cn_beta = 0.0571  # from Roskam
        reynolds_num = reynolds(air_density_isa(flow_conditions.altitude), v_crit, l_f)

        k_n = 0.01 * (0.27 * (x_cog / l_f) - 0.168 * np.log(l_f / d_f) + 0.416) - 0.0005
        k_rj = 0.46 * np.log10(reynolds_num / 10 ** 6) + 1

        swp = np.radians(sweep50)

        cn_beta_f = - 360 / (2 * np.pi) * k_n * k_rj * (l_f ** 2 * d_f) / (s_wing * b_w)
        cy_beta_v = -1 * (2 * np.pi * aspect_v) / (2 + np.sqrt(aspect_v ** 2 * (1 + np.tan(swp) ** 2 - mach ** 2) + 4))

        s_ratio = ((cn_beta - cn_beta_f) / (- cy_beta_v)) * b_w / l_v

        s_vertical = s_ratio * s_wing

        return s_vertical

    if aircraft_type == 'DUUC':


        return s_vertical


def s_vertical_sized(aircraft_type, s_wing, x_cog, l_f, d_f, b_w, l_v, v_crit, aspect_v, sweep50, mach, sweep25, power, eta,
                     approach_velocity, n_engines, y_engine):

    s_stab = s_stability(aircraft_type, s_wing, x_cog, l_f, d_f, b_w, l_v, v_crit, aspect_v, sweep50, mach)
    s_cont = s_control(aircraft_type, sweep25, l_v, power, eta, approach_velocity, n_engines, y_engine)

    s_vert = max(s_stab, s_cont)
    return s_vert
