import numpy as np


def slopes(type_sizing: str, aircraft_type: str, zh, sweep, aspect_ratio, lf, x_lemac, cmac, eta_h, cl_h, cl, tr,
           span, sweep50, mach, aspect_ratio_h, swp50h):
    a = 0
    b = 0

    if type_sizing == 'control':
        if aircraft_type == 'conventional':
            dcl_flapped = 1.089
            c_r = 1.173
            xcp_cmac = 0.44
            dcm_0 = -0.001
            et = np.radians(-3)

            l_h = lf - x_lemac - 1.5

            cm_e = -1.97 / (0.5 * 1.225 * 64.95 ** 2 * 62.187 * 2.2345) * (0.73 * 4102e3) / 64.95

            cm_0 = -.1
            dcm = dcl_flapped * (0.25 - (xcp_cmac * c_r))
            cm_0_flapped = cm_0 + dcm

            swp = np.radians(sweep)
            aratio = aspect_ratio * np.cos(swp) ** 2 / aspect_ratio + np.cos(swp)

            cm_w = cm_0_flapped * aratio + (dcm_0 / et) * et

            a = cl / (cl_h * eta_h * (l_h / cmac))

            b = (cm_w + cm_e) / (cl_h * eta_h * (l_h / cmac))

        elif aircraft_type == 'DUUC':
            a = 1
            b = 2

    elif type_sizing == 'stability':
        if aircraft_type == 'conventional':
            l_h = lf - x_lemac - 1.5

            ka = 1 / aspect_ratio - 1 / (1 - aspect_ratio ** 1.7)
            kt = (10 - 3 * tr) / 7
            kh = (1 - np.abs(zh / span)) / np.cbrt(2 * l_h / span)

            swp = np.radians(sweep)
            swp50 = np.radians(sweep50)

            root = np.sqrt(aspect_ratio ** 2 * (1 + np.tan(swp50) ** 2 - mach ** 2) + 4)
            root0 = np.sqrt(aspect_ratio ** 2 * (1 + np.tan(swp50) ** 2) + 4)

            cl_a_m = (np.pi * 2 * aspect_ratio) / (2 + root)
            cl_a_m0 = (np.pi * 2 * aspect_ratio) / (2 + root0)

            de_da = 4.44 * (ka * kt * kh * np.sqrt(np.cos(swp))) ** 1.19 * (cl_a_m / cl_a_m0)

            swp50h = np.radians(swp50h)
            rooth = np.sqrt(aspect_ratio_h ** 2 * (1 + np.tan(swp50h) ** 2 - mach ** 2) + 4)

            cl_a_h = (2 * np.pi * aspect_ratio_h) / (2 + rooth)

            a = cl_a_m / (cl_a_h * eta_h * (1 - de_da) * (l_h / cmac))

            b = 0
        elif aircraft_type == 'DUUC':
            a = 1
            b = 0

    return a, b


