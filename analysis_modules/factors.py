import numpy as np
import math


def mach_correction(mach):
    f_m = 1 - 0.08 * mach ** 1.45
    return f_m


def skin_friction(Re, flow_characteristic):
    if flow_characteristic == 't':
        c_f = 0.455/(math.pow(math.log10(Re), 2.58))
        return c_f
    else:
        c_f = 1.327/(Re ** 0.5)
        return c_f
