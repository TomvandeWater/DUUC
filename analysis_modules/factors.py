import numpy as np
import math


def mach_correction(mach):
    f_m = 1 - 0.08 * mach ** 1.45
    print("\n----- Mach correction factor:")
    print(f"Mach correction factor =", f_m)
    return f_m


def skin_friction(Re, flow_characteristic):
    if flow_characteristic == 't':
        c_f = 0.455/(math.pow(math.log10(Re), 2.58))
        print("\n----- Skin Friction:")
        print(f"Skin friction =", c_f)
        return c_f
    else:
        c_f = 1.327/(Re ** 0.5)
        print("\n----- Skin Friction:")
        print(f"Skin friction =", c_f)
        return c_f


def advance_ratio(airspeed, RPM, fan_diameter):
    j = airspeed / (RPM/60*fan_diameter)
    print("\n----- Advance Ratio:")
    print(f"Advance ratio =", j, "[-]")
    return j
