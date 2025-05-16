import numpy as np
import flow_conditions
import config


def slopes(type_sizing: str, aircraft_type: str, zh, sweep, aspect_ratio, lf, x_lemac, cmac, eta_h, cl_h, cl_w, tr,
           span, sweep50, mach, aspect_ratio_h, swp50h, cm_h, thrust, v_inf, s_wing, z_e, cl_a_duuc, cl_duuc):
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

            cm_e = z_e / (0.5 * flow_conditions.rho * v_inf ** 2 * s_wing * cmac) * (0.73 * 4102e3) / 64.95

            cm_0 = -.1
            dcm = dcl_flapped * (0.25 - (xcp_cmac * c_r))
            cm_0_flapped = cm_0 + dcm

            swp = np.radians(sweep)
            aratio = aspect_ratio * np.cos(swp) ** 2 / aspect_ratio + np.cos(swp)

            cm_w = cm_0_flapped * aratio + (dcm_0 / et) * et

            a = cl_w / (cl_h * eta_h * (l_h / cmac))

            b = (cm_w + cm_e) / (cl_h * eta_h * (l_h / cmac))

        elif aircraft_type == 'DUUC':
            dcl_flapped = 1.089  # revise how this is calculated for both
            c_r = 1.173  # revise how this is calculated
            xcp_cmac = 0.44
            dcm_0 = -0.001
            et = np.radians(-3)  # revise how this is calculated for both

            l_h = lf - x_lemac - 1.5

            cm_e = (z_e * thrust * eta_h) / (0.5 * flow_conditions.rho * v_inf ** 2 * s_wing * cmac)

            cm_0 = -.1
            dcm = dcl_flapped * (0.25 - (xcp_cmac * c_r)) # also revise this
            cm_0_flapped = cm_0 + dcm

            swp = np.radians(sweep)
            aratio = aspect_ratio * np.cos(swp) ** 2 / aspect_ratio + np.cos(swp)

            cm_w = cm_0_flapped * aratio + (dcm_0 / et) * et

            a = cl_w / (cl_duuc * eta_h * (l_h / cmac))

            b = (cm_w + cm_e + cm_h) / (cl_duuc * eta_h * (l_h / cmac))

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

            cl_a_h = cl_a_duuc

            a = cl_a_m / (cl_a_h * eta_h * (1 - de_da) * (l_h / cmac))
            b = 0

    return a, b


def horizontal_tail_control(cl, cl_h, eta_h, l_h, c_mac, cm_w, cm_e):
    a_control = cl / (cl_h * eta_h * (l_h / c_mac))
    b_control = (cm_w + cm_e) / (cl_h * eta_h * (l_h / c_mac))
    return a_control, b_control


def horizontal_tail_stability(cl_a_w, cl_a_h, eta_h, de_da, l_h, c_mac):
    a_stability = cl_a_w / (cl_a_h * eta_h * (1 - de_da) * (l_h / c_mac))
    return a_stability


def de_da_horizontal(ar, tr, z_h, wing_span, l_h, phi_qc, cl_a_m0, cl_a_m):
    phi_qc_rad = np.radians(phi_qc)
    k_a = (1 / ar) - (1 - (1 - ar ** 1.7))
    k_lamda = (10 - 3 * tr) / 7
    k_h = (1 - abs(z_h / wing_span)) / (np.cbrt(2 * l_h / wing_span))

    de_da = 4.44 * (k_a * k_lamda * k_h * np.cos(phi_qc_rad)) ** 1.19 * (cl_a_m / cl_a_m0)
    return de_da


