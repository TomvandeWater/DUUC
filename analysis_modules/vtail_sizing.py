import numpy as np
import flow_conditions
from analysis_modules.aerodynamic import reynolds
from analysis_modules.ISA import air_density_isa
import data.atr_reference as ref
import matplotlib.pyplot as plt


def s_control(aircraft_type, sweep25, l_v, power, eta, approach_velocity, n_engines, y_engine, cy_duuc=0, cd_pe=0.00,
              cd_wind=0.00, cy_atr=0):
    velocity_mc = (1.3 / 1.2) * approach_velocity
    delta_f_max = 25
    swp = np.radians(sweep25)
    rho_sea = 1.225

    if aircraft_type == 'conventional':
        k_theta = (1 - 0.08 * np.cos(swp) ** 2) * np.cos(swp) ** (3/4)
        cl_the = 4.52  # from textNita -> still make this an interpollation based on input
        k_prime = 0.675
        cld_ratio = 0.85

        ne = ((eta * power) / (velocity_mc * n_engines)) * y_engine
        nd = 0.25 * ne
        top = ne + nd

        delta_f = np.radians(delta_f_max)
        denominator = (0.5 * rho_sea * velocity_mc ** 2 * delta_f * cld_ratio * cl_the * k_prime * k_theta * l_v)
        denominator = 0.5 * rho_sea * velocity_mc ** 2 * l_v * cy_atr

    elif aircraft_type == 'DUUC':
        ne = ((eta * power) / (velocity_mc * n_engines)) * y_engine
        nd = ((cd_pe * 0.5 * rho_sea * velocity_mc ** 2 * ref.s_w) + cd_wind * 0.5 * rho_sea * velocity_mc ** 2 * 1.52) * y_engine
        nd = ne * 0.25
        top = ne + nd

        denominator = 0.5 * rho_sea * velocity_mc ** 2 * l_v * cy_duuc
    else:
        raise ValueError("Invalid aircraft type specified. Choose 'conventional' or 'DUUC'.")

    s_vertical = (top * 2 * np.pi) / denominator
    print(f"ne: {ne}, nd: {nd}")
    print(f"v approach: {approach_velocity}")
    print(f"s_vertical: {s_vertical}")

    return s_vertical

cyd = 0.975
#a = s_control("conventional", 25, 9.13, 2051*10**3, 0.73, 60.4, 2, 3.608)
a = 13.788741676315057
b = s_control("DUUC", 25, 9.13, 2051*10**3, 0.73, 60.4, 2, 2.8,
          cy_duuc=cyd, cd_pe=0.00568, cd_wind=0.06)





def s_stability(aircraft_type, s_wing, x_cog, l_f, d_f, b_w, l_v, v_crit, aspect_v, sweep50, mach, cy_duuc=0):
    cn_beta = 0.0571  # 1/rad from Roskam
    reynolds_num = (v_crit * l_f) / (5.4603 * 10 ** -5)

    k_n = 0.01 * (0.27 * (x_cog / l_f) - 0.168 * np.log(l_f / d_f) + 0.416) - 0.0005
    k_rj = 0.46 * np.log10(reynolds_num / 10 ** 6) + 1

    cn_beta_f = - 360 / (2 * np.pi) * k_n * k_rj * (l_f ** 2 * d_f) / (s_wing * b_w)

    if aircraft_type == 'conventional':
        cy_beta_v = -1 * (2 * np.pi * aspect_v) / (2 + np.sqrt(aspect_v ** 2 * (1 + 0.31 ** 2 - mach ** 2) + 4))
    elif aircraft_type == 'DUUC':
        cy_beta_v = -cy_duuc
    else:
        raise ValueError("Invalid aircraft type specified. Choose 'conventional' or 'DUUC'.")

    s_ratio = ((cn_beta - cn_beta_f) / (- cy_beta_v)) * (b_w / l_v)
    print(f"reynolds_num: {reynolds_num}")
    print(f"kn: {k_n}, krj: {k_rj}")
    print(f"cn_beta_f: {cn_beta_f}")
    print(f"cy_beta_v: {cy_beta_v}")
    print(f"b_w: {b_w}, lv: {l_v}")

    s_vertical = s_ratio * s_wing
    print(f"s_vertical: {s_vertical}")
    return s_vertical


s_stability("conventional", ref.s_w, 11.5586, 27.13, 2.77, ref.b_w, 9.13, 141, ref.ar_v, 0.31, 0.44)


def s_vertical_sized(aircraft_type, s_wing, x_cog, l_f, d_f, b_w, l_v, v_crit, aspect_v, sweep50, mach, sweep25, power, eta,
                     approach_velocity, n_engines, y_engine):

    s_stab = s_stability(aircraft_type, s_wing, x_cog, l_f, d_f, b_w, l_v, v_crit, aspect_v, sweep50, mach)
    s_cont = s_control(aircraft_type, sweep25, l_v, power, eta, approach_velocity, n_engines, y_engine)

    s_vert = max(s_stab, s_cont)
    return s_vert


array = np.linspace(0.1, 3.0, 301)
s_array = []
s_array_atr = []
s_stab_duuc = []
s_stab_atr = []

for i in range(len(array)):
    s_array.append(s_control("DUUC", 25, 9.13, 2051*10**3, 0.73, 60.4, 2, 2.8,
          cy_duuc=array[i], cd_pe=0.00568, cd_wind=0.06))
    s_array_atr.append(s_control("conventional", 25, 9.13, 2051*10**3, 0.73, 60.4, 2, 3.608,
                                 cy_atr=array[i]))
    s_stab_duuc.append(s_stability("DUUC", ref.s_w, 11.5586, 27.13, 2.77, ref.b_w, 9.13,
                                   141, ref.ar_v, 0.31, 0.44, cy_duuc=array[i]))
    s_stab_atr.append(s_stability("DUUC", ref.s_w, 11.5586, 27.13, 2.77, ref.b_w, 9.13,
                                   141, ref.ar_v, 0.31, 0.44, cy_duuc=array[i]))

plt.figure('CY - vs S_V')
plt.plot(array, s_array, label=r'Prediction line DUUC - control', color="tab:blue")
plt.plot(array, s_array_atr, label=r'Prediction line ATR - control', color="tab:orange")
plt.plot([0, 0.975], [a, a], color="tab:orange", linestyle="dashed")
plt.plot(0.975, a, 'o', markersize=6, color="tab:orange", label="ATR - value")
plt.plot([0.975, 0.975], [0, a], linestyle='dashed', color="tab:orange")
plt.plot([0, cyd], [b, b], color="tab:blue", linestyle="dashed")
plt.plot(cyd, b, 'o', markersize=6, color="tab:blue", label="DUUC - 1st iteration")
plt.plot([cyd, cyd], [0, b], color="tab:blue", linestyle="dashed")
plt.xlim(0, 2)
plt.ylim(0, 100)
plt.xlabel(r'$C_{Y}$ [-]')
plt.ylabel(r'$S_{V}$ [m2]')
plt.title(r'$S_{V}$ vs. $C_{Y}$')
plt.legend()
plt.grid(True)

plt.figure('CYbeta - vs S_V - stability')
plt.plot(array, s_stab_duuc, label=r'Prediction line DUUC - stability', color="tab:blue")
plt.plot([0, 2.228], [15.01, 15.01], color="tab:orange", linestyle="dashed")
plt.plot(2.228, 15.01, 'o', markersize=6, color="tab:orange", label="ATR - value")
plt.plot([2.228, 2.228], [0, 15.01], linestyle='dashed', color="tab:orange")
#plt.plot([0, 1.321], [b, b], color="tab:blue", linestyle="dashed")
#plt.plot(1.321, b, 'o', markersize=6, color="tab:blue", label="DUUC - 1st iteration")
#plt.plot([1.321, 1.321], [0, b], color="tab:blue", linestyle="dashed")
plt.xlim(1, 3)
plt.ylim(0, 100)
plt.xlabel(r'$C_{Y_{\beta}}$ [-]')
plt.ylabel(r'$S_{V}$ [m2]')
plt.title(r'$S_{V}$ vs. $C_{Y}$')
plt.legend()
plt.grid(True)

plt.show()


