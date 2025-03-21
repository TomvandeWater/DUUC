import numpy as np
import math


def mach_correction(mach):
    f_m = 1 - 0.08 * mach ** 1.45
    # print("\n----- Mach correction factor:")
    # print(f"Mach correction factor =", f_m)
    return f_m


def skin_friction(Re, flow_characteristic):
    if flow_characteristic == 't':
        c_f = 0.455/(math.pow(math.log10(Re), 2.58))
        # print("\n----- Skin Friction:")
        # print(f"Skin friction =", c_f)
        return c_f
    else:
        c_f = 1.327/(Re ** 0.5)
        # print("\n----- Skin Friction:")
        # print(f"Skin friction =", c_f)
        return c_f


def advance_ratio(airspeed, RPM, fan_diameter):
    j = airspeed / (RPM / 60 * fan_diameter)
    # print("\n----- Advance Ratio:")
    # print(f"Advance ratio =", j, "[-]")
    return j


def k_control(sweep, aspect_ratio):
    if aspect_ratio > 4:
        k = 1 + ((8.2-2.3 * sweep) - (0.22 - 0.153 * sweep) * aspect_ratio) / 100
        return k
    else:
        k = 1 + ((1.87 - 0.000233 * sweep) * aspect_ratio) / 100
        return k


def oswald(aspect_ratio, sweep):
    if sweep == 0:
        e = 1 / (1 + 0.38 / aspect_ratio + 60 / aspect_ratio ** 3)
        return e
    else:
        e = 1.78 * (1 - 0.045 * aspect_ratio ** 0.68) / np.cos(np.radians(sweep)) ** 0.25
        return e


def area_ratio(nacaprofile, chord, radius, x_location):
    num_list = [int(digit) for digit in nacaprofile]
    thickness = num_list[2] * 10 + num_list[3]  # naca 4 series thickness of profile
    t = thickness / 100
    x = x_location / chord

    y = (t / 0.2) * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1015 * x**4)

    diameter = 2 * (radius - (y * chord))

    return diameter

