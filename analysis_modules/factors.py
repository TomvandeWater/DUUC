import matplotlib.pyplot as plt
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

    radius_out = radius - (y * chord)
    return diameter, radius_out


def cross_sectional_area(x, r_input, x_prop):
    if 0 <= x <= x_prop:
        radius = 0
    elif x_prop < x <= 0.8:
        radius = r_input
    elif 0.8 < x <= 1:
        # Linear decrease from r_input at x=0.8 to 0 at x=1
        radius = r_input * (1 - (x - 0.8) / (1 - 0.8))
    else:
        return 0  # Outside valid x range, return 0

    # Compute area
    area = math.pi * radius ** 2
    return area, radius


def mass_of_section(m_total, l, x1, x2):
    # Validate input
    if not (0 <= x1 <= l and 0 <= x2 <= l and x1 <= x2):
        raise ValueError("x1 and x2 must be within [0, L] and x1 <= x2")

    factor = (2 * m_total) / l
    term1 = x2 - (x2 ** 2) / (2 * l)
    term2 = x1 - (x1 ** 2) / (2 * l)

    m_section = factor * (term1 - term2)
    return m_section

"""
x_c = np.linspace(0, 1, 101)
area = []
area_av = []
r_r1 = []
r_r2 = []

for i in range(len(x_c)):
    area.append(np.pi * area_ratio("0012", 1, 1.8, x_c[i])[1] ** 2)
    area_av.append(np.pi * area_ratio("0012", 1, 1.8, x_c[i])[1] ** 2 - cross_sectional_area(x_c[i], 0.3)[0])
    r_r1.append(0.5 * area_ratio("0012", 1, 1.8, x_c[i])[0]/1.8)
    r_r2.append((0.5 * area_ratio("0012", 1, 1.8, x_c[i])[0] - cross_sectional_area(x_c[i], 0.3)[1])/1.8)


plt.figure('Area vs. chord')
plt.plot(x_c, area, label=r'Duct area')
plt.plot(x_c, area_av, label=r'Duct minus Nacelle area')
plt.axvline(x=0.3, linestyle='--', color='black', label='Propeller location')
plt.xlabel(r'x/c [-]')
plt.ylabel(r'Cross sectional area [$m^2$]')
plt.title(r'Area vs. chord (NACA0012)')
plt.ylim([7, 10.5])
plt.legend()
plt.grid(True)

plt.figure('Radius duct vs. Radius available')
plt.plot(x_c, r_r1, label=r'Duct radius')
plt.plot(x_c, r_r2, label=r'Duct minus Nacelle radius')
plt.axvline(x=0.3, linestyle='--', color='black', label='Propeller location')
plt.axhline(y=1.0, linestyle='--', color='tab:blue', alpha=0.5, label='Duct airfoil centerline')
plt.xlabel(r'x/c [-]')
plt.ylabel(r'$R_{i}$ / $R_{center}$ []')
plt.title(r'Area vs. chord (NACA0012)')
plt.legend()
plt.grid(True)

plt.show() """
